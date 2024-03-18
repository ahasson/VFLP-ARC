#!/usr/bin/env python3

# Copyright (C) 2019 Christoph Gorgulla
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
#
# This file is part of VirtualFlow.
#
# VirtualFlow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# VirtualFlow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with VirtualFlow.  If not, see <https://www.gnu.org/licenses/>.

# ---------------------------------------------------------------------------
#
# Description: Main runner for the individual workunits/subjobs
#              for ligand preparation
#
# Revision history:
# 2021-09-13  Initial version of VFLP ported to Python
#
# ---------------------------------------------------------------------------
import atexit
import csv
import gzip
import json
import logging
import os
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import time
import typing as t


from run.process import process_collection
from run.tranche import get_mol_attributes


AWS_AVAILABLE: bool = False
try:
    import boto3
    import botocore
    from botocore.config import Config

    AWS_AVAILABLE = True
except ImportError:
    logging.warning("AWS: boto3 not available")


####### Main thread

# We now expect our config file to be a JSON file


def parse_config(filename):

    with open(filename, "rt", encoding="utf8") as read_file:
        config = json.load(read_file)

    return config


def is_cxcalc_used(ctx):
    used = False

    if ctx["main_config"]["tautomerization"] == "true" and (
        ctx["main_config"]["tautomerization_program_1"] == "cxcalc"
        or ctx["main_config"]["tautomerization_program_2"] == "cxcalc"
    ):
        used = True

    if ctx["main_config"]["protonation_state_generation"] == "true" and (
        ctx["main_config"]["protonation_program_1"] == "cxcalc"
        or ctx["main_config"]["protonation_program_2"] == "cxcalc"
    ):
        used = True

    # In addition to these, there are some tranche assignment functions that use
    # cxcalc. Original code does not cover these, but should be addressed

    return used


def is_molconvert_used(ctx):
    used = False
    if ctx["main_config"]["conformation_generation"] == "true" and (
        ctx["main_config"]["conformation_program_1"] == "molconvert"
        or ctx["main_config"]["conformation_program_2"] == "molconvert"
    ):
        used = True

    return used


def is_standardizer_used(ctx):
    used = False
    if ctx["main_config"]["neutralization"] == "true" and (
        ctx["main_config"]["neutralization_program_1"] == "standardizer"
        or ctx["main_config"]["neutralization_program_2"] == "standardizer"
    ):
        used = True

    return used


def is_nailgun_needed(ctx):
    return is_cxcalc_used(ctx) or is_molconvert_used(ctx) or is_standardizer_used(ctx)


# TODO: @Alex -> Chrisopth - this isn't used anywhere
def is_rdkit_needed(ctx):

    attributes = get_mol_attributes()
    if ctx["main_config"]["tranche_assignments"] == "true":
        for tranche_type in ctx["main_config"]["tranche_types"]:
            if tranche_type in attributes:
                if attributes[tranche_type]["prog"] == "rdkit":
                    return True

    elif ctx["main_config"]["stereoisomer_generation"] == "true" and (
        ctx["main_config"]["stereoisomer_generation_program_1"] == "rdkit"
        or ctx["main_config"]["stereoisomer_generation_program_2"] == "rdkit"
    ):
        return True

    return False


def get_helper_versions(ctx):

    ctx["versions"] = {
        "cxcalc": "INVALID",
        "molconvert": "INVALID",
        "standardizer": "INVALID",
        "obabel": "INVALID",
    }
    nailgun_host = ctx["main_config"]["nailgun_host"]
    nailgun_port = ctx["main_config"]["nailgun_port"]
    command_grep: str = r"grep -m 1 version | sed 's/.*version \([0-9. ]*\).*/\\1/'"

    if is_cxcalc_used(ctx):
        command_cxcalc: str = f"chemaxon.marvin.Calculator | {command_grep}"
        ret = subprocess.run(
            f"ng --nailgun-server {nailgun_host} --nailgun-port {nailgun_port} {command_cxcalc}",
            capture_output=True,
            text=True,
            shell=True,
            timeout=15,
            check=False,
        )
        if ret.returncode != 0:
            logging.error("Cannot connect to nailgun server for cxcalc version")
            raise RuntimeError("Cannot connect to nailgun server for cxcalc version")
        ctx["versions"]["cxcalc"] = ret.stdout.strip()
        logging.error(  # TODO: Should this be an error?
            "%s:%s: cxcalc version is %s",
            nailgun_port,
            nailgun_host,
            ctx["versions"]["cxcalc"],
        )

        # TODO: @Alex -> Christoph - this is a duplicate of the above code
        # ret = subprocess.run(
        #    f"ng --nailgun-server {nailgun_host} --nailgun-port {nailgun_port} chemaxon.marvin.Calculator",
        #    capture_output=True,
        #    text=True,
        #    shell=True,
        #    timeout=15,
        #    check=False,
        # )
        # logging.error(
        #    "CXCALC version: %s, %s, %s", ret.returncode, ret.stderr, ret.stdout
        # )

    if is_molconvert_used(ctx):
        command_molconv: str = f"chemaxon.formats.MolConverter | {command_grep}"
        ret = subprocess.run(
            f"ng --nailgun-server {nailgun_host} --nailgun-port {nailgun_port} {command_molconv}",
            capture_output=True,
            text=True,
            shell=True,
            timeout=15,
            check=False,
        )
        if ret.returncode != 0:
            raise RuntimeError("Cannot connect to nailgun for molconvert version")
        ctx["versions"]["molconvert"] = ret.stdout.strip()

    if is_standardizer_used(ctx):
        command_standard: str = (
            "chemaxon.standardizer.StandardizerCLI -h | head -n 1 | awk -F '[ ,]' '{{print $2}}'"
        )
        ret = subprocess.run(
            f"ng --nailgun-server {nailgun_host} --nailgun-port {nailgun_port} {command_standard}",
            capture_output=True,
            text=True,
            shell=True,
            timeout=15,
            check=False,
        )
        if ret.returncode != 0:
            raise RuntimeError("Cannot connect to nailgun for standardizer version")
        ctx["versions"]["standardizer"] = ret.stdout.strip()

    # OpenBabel is required across all configurations
    ret = subprocess.run(
        "obabel -V | awk '{{print $3}}'",
        capture_output=True,
        text=True,
        shell=True,
        timeout=15,
        check=False,
    )
    if ret.returncode != 0:
        raise RuntimeError("Cannot get obabel version")
    ctx["versions"]["obabel"] = ret.stdout.strip()


# This has the potential for a race condition. Better
# solution is if Nailgun could report port that it is using


def get_free_port():
    sock = socket.socket()
    sock.bind(("", 0))
    return sock.getsockname()[1]


def start_ng_server(host, port_num, java_max_heap_size):
    # TODO: Should this be a logging.error?
    logging.error("Nailgun starting on %s:%s", host, port_num)
    cmds = [
        "java",
        f"-Xmx{java_max_heap_size}G",
        "com.martiansoftware.nailgun.NGServer",
        f"{host}:{port_num}",
    ]
    ng_process = subprocess.Popen(cmds, start_new_session=True)

    # Nailgun can take a while to come up
    #
    time.sleep(10)

    return ng_process


def cleanup_ng(ng_process):
    # Stop the Nailgun serve
    ng_process.kill()
    print("Nailgun stopped")


def get_workunit_information() -> t.Tuple[str, str]:

    workunit_id = os.getenv("VFLP_WORKUNIT", "")
    subjob_id = os.getenv("VFLP_WORKUNIT_SUBJOB", "")

    if workunit_id == "" or subjob_id == "":
        raise RuntimeError("Invalid VFLP_WORKUNIT and/or VFLP_WORKUNIT_SUBJOB")

    return workunit_id, subjob_id


def get_subjob_config(ctx, workunit_id, subjob_id):

    ctx["job_storage_mode"] = os.getenv("VFLP_JOB_STORAGE_MODE", "INVALID")

    if ctx["job_storage_mode"] == "s3":
        # Get the initial bootstrap information
        job_object = os.getenv("VFLP_CONFIG_JOB_OBJECT")
        job_bucket = os.getenv("VFLP_CONFIG_JOB_BUCKET")
        download_to_workunit_file = f"{ctx['temp_dir'].name}/{workunit_id}.json.gz"

        # Download workunit from S3
        get_workunit_from_s3(
            ctx,
            workunit_id,
            subjob_id,
            job_bucket,
            job_object,
            download_to_workunit_file,
        )

        # Get the subjob config and main config information (same file)
        try:
            with gzip.open(download_to_workunit_file, "rt") as f:
                subjob_config = json.load(f)
                ctx["main_config"] = subjob_config["config"]

                if subjob_id in subjob_config["subjobs"]:
                    ctx["subjob_config"] = subjob_config["subjobs"][subjob_id]
                else:
                    logging.error("There is no subjob ID with ID: %s", subjob_id)
                    # AWS Batch requires that an array job have at least 2 elements,
                    # sometimes we only need 1 though
                    if subjob_id == "1":
                        exit(0)
                    else:
                        raise RuntimeError(f"There is no subjob ID with ID:{subjob_id}")

        except Exception as err:
            logging.error("Cannot open %s: %s", download_to_workunit_file, str(err))
            raise

        ctx["main_config"]["job_object"] = job_object
        ctx["main_config"]["job_bucket"] = job_bucket

    elif ctx["job_storage_mode"] == "sharedfs":

        config_json = os.getenv("VFLP_CONFIG_JSON", "")
        workunit_json = os.getenv("VFLP_WORKUNIT_JSON", "")

        if config_json == "" or workunit_json == "":
            print(
                "For VFLP_JOB_STORAGE_MODE=sharedfs, VFLP_CONFIG_JSON, VFLP_WORKUNIT_JSON must be set"
            )
            sys.exit(1)

        # load in the main config from the ENV
        try:
            with open(config_json, "rt", encoding="utf8") as f:
                config = json.load(f)
                ctx["main_config"] = config
        except Exception as err:
            logging.error("Cannot open %s: %s", config_json, str(err))
            raise

        # Load workunit information
        try:
            with gzip.open(workunit_json, "rt") as f:
                subjob_config = json.load(f)
                if subjob_id in subjob_config["subjobs"]:
                    ctx["subjob_config"] = subjob_config["subjobs"][subjob_id]
                else:
                    logging.error(
                        "There is no subjob ID with ID:%s in %s", subjob_id, workunit_id
                    )
                    raise RuntimeError(
                        f"There is no subjob ID with ID:{subjob_id} in {workunit_id}"
                    )

        except Exception as err:
            logging.error("Cannot open %s: %s", workunit_json, str(err))
            raise

        # Set paths used later by the sharedfs mode
        ctx["workflow_dir"] = ctx["main_config"]["sharedfs_workflow_path"]
        ctx["collection_dir"] = ctx["main_config"]["sharedfs_collection_path"]

    else:
        raise RuntimeError(
            f"Invalid jobstoragemode of {ctx['job_storage_mode']}. VFLP_JOB_STORAGE_MODE must be 's3' or 'sharedfs' "
        )


# Get only the collection information with the subjob specified


def get_workunit_from_s3(
    ctx, _workunit_id, _subjob_id, job_bucket, job_object, workunit_file
):
    try:
        with open(workunit_file, "wb") as f:
            ctx["s3"].download_fileobj(job_bucket, job_object, f)
    except botocore.exceptions.ClientError as error:
        logging.error(
            "Failed to download from S3 %s/%s to %s, (%s)",
            job_bucket,
            job_object,
            workunit_file,
            error,
        )
        raise


def process(ctx: dict):

    # Figure out what job we are running

    ctx["vcpus_to_use"] = int(os.getenv("VFLP_VCPUS", "1"))
    ctx["run_sequential"] = int(os.getenv("VFLP_RUN_SEQUENTIAL", "0"))

    workunit_id, subjob_id = get_workunit_information()

    # This includes all of the configuration information we need
    # After this point ctx['main_config'] has the configuration options
    # and we have specific subjob information in ctx['subjob_config']

    get_subjob_config(ctx, workunit_id, subjob_id)

    # Setup the Nailgun server [only setup if needed]

    if is_nailgun_needed(ctx):
        ng_port = get_free_port()
        ng_host = os.getenv("VFLP_HOST", "localhost")
        ng_process = start_ng_server(
            ng_host, ng_port, ctx["main_config"]["java_max_heap_size"]
        )
        ctx["main_config"]["nailgun_port"] = str(ng_port)
        ctx["main_config"]["nailgun_host"] = str(ng_host)
        atexit.register(cleanup_ng, ng_process)

    # Get information about the versions of the software
    # being used in the calculations
    get_helper_versions(ctx)

    # Run all of the collections. Do these in series to limit the amount of
    # storage needed until we push back to the filesystem or S3

    for collection_key in ctx["subjob_config"]["collections"]:
        collection = ctx["subjob_config"]["collections"][collection_key]
        collection_data = get_collection_data(ctx, collection_key)
        process_collection(ctx, collection_key, collection, collection_data)


def get_collection_data(ctx, collection_key):

    collection = ctx["subjob_config"]["collections"][collection_key]

    collection_data = {}
    collection_data["ligands"] = {}

    # subjob_config
    # 	collections (collection_key)
    # 		metatranche
    # 		tranche
    # 		collection_name
    #
    # 		s3_bucket
    # 		s3_download_path
    #
    # 		fieldnames [headers for collection file]
    #

    # The collection obj will have information on where it is located
    # Use that to create the directory (if not already created)

    if ctx["job_storage_mode"] == "s3":

        temp_dir = Path(ctx["temp_dir"].name)
        download_dir = temp_dir / collection["metatranche"] / collection["tranche"]
        download_dir.mkdir(parents=True, exist_ok=True)

        collection_file = download_dir / f"{collection['collection_name']}.txt.gz"
        try:
            with collection_file.open(mode="wb") as f:
                ctx["s3"].download_fileobj(
                    collection["s3_bucket"], collection["s3_download_path"], f
                )
        except botocore.exceptions.ClientError as error:
            logging.error(
                "Failed to download from S3 {%s}/%s to %s, (%s)",
                collection["s3_bucket"],
                collection["s3_download_path"],
                str(collection_file),
                error,
            )
            raise

    elif ctx["job_storage_mode"] == "sharedfs":
        collection_file = (
            Path(ctx["collection_dir"])
            / collection["metatranche"]
            / collection["tranche"]
            / f"{collection['collection_name']}.txt.gz"
        )

    # Read in the collection data (regardless of source)

    try:
        with gzip.open(collection_file, "rt") as f:
            reader = csv.DictReader(
                f, fieldnames=collection["fieldnames"], delimiter="\t"
            )
            for row in reader:
                ligand_key = row["ligand-name"]
                collection_data["ligands"][ligand_key] = {
                    "smi": row["smi"],
                    "file_data": row,
                }

    except Exception as err:
        logging.error("Cannot open %s: %s", str(collection_file), str(err))
        raise

    return collection_data


def check_diskfree() -> None:
    """Prints the amount of free disk space.

    Args:
            None
    Returns:
            None
    Raises:
            None
    """
    ret = subprocess.run(["df", "-h"], capture_output=True, text=True, check=False)
    if ret.returncode == 0:
        print(ret.stdout)
    else:
        print("could not run df -h")


def main() -> None:
    """Main function for the script"""

    logging.basicConfig(level=logging.ERROR)

    ctx: dict = {}

    check_diskfree()

    if AWS_AVAILABLE:
        aws_region: str = os.getenv("VFLP_REGION", "us-east-1")

        botoconfig = Config(
            region_name=aws_region, retries={"max_attempts": 15, "mode": "standard"}
        )

        # Get the config information
        ctx["s3"] = boto3.client("s3", config=botoconfig)

    tmp_path = os.getenv("VFLP_TMP_PATH", "/tmp")
    tmp_path = os.path.join(tmp_path, "")

    ctx["temp_path"] = tmp_path
    ctx["temp_dir"] = tempfile.TemporaryDirectory(prefix=tmp_path)

    process(ctx)

    print("end")
    check_diskfree()


if __name__ == "__main__":
    main()
