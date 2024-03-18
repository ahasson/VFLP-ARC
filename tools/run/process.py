import multiprocessing
import tempfile
from datetime import datetime
import logging
import json
import gzip
import os
import tarfile
import hashlib
from pathlib import Path
import shutil

from run.ligand import process_ligand

AWS_AVAILABLE: bool = False
try:
    import boto3
    import botocore
    from botocore.config import Config

    AWS_AVAILABLE = True
except ImportError:
    logging.warning("AWS: boto3 not available")


def move_file(ctx, move_item):

    if ctx["job_storage_mode"] == "s3" and AWS_AVAILABLE:
        print(
            f"moving |{move_item['local_path']}| to |{move_item['remote_path']}| - bucket:{ctx['main_config']['job_bucket']}"
        )
        try:
            _response = ctx["s3"].upload_file(
                move_item["local_path"],
                ctx["main_config"]["job_bucket"],
                move_item["remote_path"],
            )
        except botocore.exceptions.ClientError as e:
            logging.error(e)
            raise

    elif ctx["job_storage_mode"] == "sharedfs":
        print(f"moveitem_local: {move_item['remote_path']}")

        local_file_path = Path(move_item["remote_path"])
        local_file_parent = local_file_path.parent.absolute()
        local_file_parent.mkdir(parents=True, exist_ok=True)

        print(f"{local_file_path.name} -> parent: {local_file_parent.name}")

        shutil.copyfile(move_item["local_path"], move_item["remote_path"])


def generate_tarfile(dir):
    os.chdir(str(Path(dir).parents[0]))

    with tarfile.open(f"{os.path.basename(dir)}.tar.gz", "x:gz") as tar:
        tar.add(os.path.basename(dir))

    return os.path.join(str(Path(dir).parents[0]), f"{os.path.basename(dir)}.tar.gz")


def get_collection_hash(collection_key):
    string_to_hash = f"{collection_key}"
    return hashlib.sha256(string_to_hash.encode()).hexdigest()


def generate_remote_path(
    ctx, collection, output_type="status", output_format="json.gz"
):

    if ctx["job_storage_mode"] == "s3":
        base_prefix = ctx["main_config"]["object_store_job_output_data_prefix_full"]
    elif ctx["job_storage_mode"] == "sharedfs":
        base_prefix = ctx["workflow_dir"]

    if ctx["main_config"]["job_storage_output_addressing"] == "hash":
        collection_string = f"{collection['metatranche']}_{collection['tranche']}_{collection['collection_name']}"
        hash_string = get_collection_hash(collection_string)

        prefix_components = [
            base_prefix,
            hash_string[0:2],
            hash_string[2:4],
        ]
    else:
        prefix_components = [base_prefix]

    prefix = "/".join(prefix_components)

    if output_format == "json.gz":
        remote_dir = [
            prefix,
            "complete",
            output_type,
            collection["metatranche"],
            collection["tranche"],
        ]

        # Remote path
        return "/".join(remote_dir) + f"/{collection['collection_name']}.json.gz"

    elif output_format == "tar.gz":
        collection_path = [
            prefix,
            "complete",
            output_type,
            collection["metatranche"],
            collection["tranche"],
            collection["collection_name"],
        ]

        return "/".join(collection_path) + ".tar.gz"
    else:
        raise RuntimeError(f"Invalid value of output_format ({output_format})")


def generate_collection_summary_file(
    collection, collection_temp_dir, collection_summary
):

    temp_dir = Path(collection_temp_dir.name)
    status_dir = (
        temp_dir
        / "complete"
        / "status"
        / collection["metatranche"]
        / collection["tranche"]
    )
    status_dir.mkdir(parents=True, exist_ok=True)

    status_file = status_dir / f"{collection['collection_name']}.json.gz"

    with gzip.open(status_file, "wt") as json_gz:
        json.dump(collection_summary, json_gz)

    return status_file


def process_collection(ctx, collection_key, collection, collection_data):

    collection_temp_dir = tempfile.TemporaryDirectory(prefix=ctx["temp_path"])

    tasklist = []

    for ligand_key in collection_data["ligands"]:
        ligand = collection_data["ligands"][ligand_key]

        task = {
            "collection_key": collection_key,
            "metatranche": collection["metatranche"],
            "tranche": collection["tranche"],
            "collection_name": collection["collection_name"],
            "ligand_key": ligand_key,
            "ligand": ligand,
            "config": ctx["main_config"],
            "collection_temp_dir": collection_temp_dir,
            "main_temp_dir": ctx["temp_dir"],
            "versions": ctx["versions"],
            "temp_path": ctx["temp_path"],
            "store_all_intermediate_logs": ctx["main_config"][
                "store_all_intermediate_logs"
            ],
        }

        tasklist.append(task)

    start_time = datetime.now()

    res = []

    if ctx["run_sequential"] == 1:
        for taskitem in tasklist:
            res.append(process_ligand(taskitem))
    else:
        with multiprocessing.Pool(processes=ctx["vcpus_to_use"]) as pool:
            res = pool.map(process_ligand, tasklist)

    end_time = datetime.now()
    difference_time = end_time - start_time

    print(f"time difference is: {difference_time.seconds} seconds")

    # For each completed task, summarize the data into data
    # structures that we can save for later processing
    # and analysis

    unit_failed_count = 0
    unit_success_count = 0
    tautomer_failed_count = 0
    tautomer_success_count = 0

    collection_summary = {"ligands": {}, "seconds": difference_time.seconds}

    # Generate summaries based on collection from each of the results

    for task_result in res:

        ligand_key = task_result["base_ligand"]["key"]

        collection_summary["ligands"][ligand_key] = {
            "timers": task_result["base_ligand"]["timers"],
            "status": task_result["status"],
            "status_sub": task_result["base_ligand"]["status_sub"],
            "tautomers": task_result["ligands"],
            "stereoisomers": task_result["stereoisomers"],
        }
        collection_summary["ligands"][ligand_key]["seconds"] = task_result["seconds"]

        # Remove fields we do not need from the output
        for _tautomer_key, tautomer in collection_summary["ligands"][ligand_key][
            "tautomers"
        ].items():
            tautomer.pop("intermediate_dir", "")
            tautomer.pop("pdb_file", "")

    move_list = []

    # Generate the status file
    status_file = generate_collection_summary_file(
        collection, collection_temp_dir, collection_summary
    )
    status_file_remote = generate_remote_path(
        ctx, collection, output_type="status", output_format="json.gz"
    )

    move_item = {
        "local_path": status_file.as_posix(),
        "remote_path": status_file_remote,
    }

    move_list.append(move_item)

    # Save the output files from the calculations to
    # the location listed in the configuration file
    #
    # In all cases the data is tar.gz for folders

    # Process each target format into a separate file
    for target_format in ctx["main_config"]["target_formats"]:

        local_path_dir = [
            collection_temp_dir.name,
            "complete",
            target_format,
            collection["metatranche"],
            collection["tranche"],
            collection["collection_name"],
        ]

        # In some cases we may not generate any output (for example if
        # all ligands in a collection fail)
        try:
            tar_gz_path = generate_tarfile("/".join(local_path_dir))

            move_item = {
                "local_path": tar_gz_path,
                "remote_path": generate_remote_path(
                    ctx, collection, output_type=target_format, output_format="tar.gz"
                ),
            }

            move_list.append(move_item)
            print(move_item)
        except FileNotFoundError as error:
            logging.error(f"Could not create tarball from {local_path_dir}")

    # We may want to save the intermediate data for debugging

    if ctx["main_config"]["store_all_intermediate_logs"] == "true":
        logging.error("Saving intermediate logs")
        for collection_key, collection in ctx["subjob_config"]["collections"].items():

            local_path_dir = [
                collection_temp_dir.name,
                "intermediate",
                collection["metatranche"],
                collection["tranche"],
                collection["collection_name"],
            ]

            # In some cases we may not generate any output (for example if
            # all ligands in a collection fail)
            try:
                tar_gz_path = generate_tarfile("/".join(local_path_dir))

                move_item = {
                    "local_path": tar_gz_path,
                    "remote_path": generate_remote_path(
                        ctx,
                        collection,
                        output_type="intermediate",
                        output_format="tar.gz",
                    ),
                }

                move_list.append(move_item)
                print(move_item)
            except FileNotFoundError as error:
                logging.error(f"Could not create tarball from {local_path_dir}")

    print(move_list)

    for move_item in move_list:
        move_file(ctx, move_item)
