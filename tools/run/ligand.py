import time
import tempfile
import logging

from run.desalt import desalt
from run.neutral import neutralization
from run.omer import stereoisomer_generation, process_stereoisomer
from run.utils import get_intermediate_dir_ligand


def process_ligand(ctx):

    start_time = time.perf_counter()

    ctx["temp_dir"] = tempfile.TemporaryDirectory(prefix=ctx["temp_path"])
    ctx["intermediate_dir"] = get_intermediate_dir_ligand(ctx)

    completion_event = {
        "status": "failed",
        "base_ligand": ctx["ligand"],
        "ligands": {},
        "stereoisomers": {},
        "seconds": 0,
    }

    base_ligand = completion_event["base_ligand"]
    base_ligand["key"] = ctx["ligand_key"]
    # base_ligand["stereoisomers"]: [] TODO: @Alex -> Christoph, this is never used and in fact not defined properly
    base_ligand["timers"] = []
    base_ligand["status_sub"] = []
    base_ligand["remarks"] = {
        "collection_key": f"Original Collection: {ctx['collection_key']}"
    }

    print(f"* Processing {base_ligand['key']}")

    if base_ligand["smi"] == "":
        logging.error(
            f"    * Warning: Ligand {base_ligand['key']} skipped since SMI is blank"
        )
        return completion_event

    # De-salting
    step_timer_start = time.perf_counter()
    try:
        desalt(ctx, base_ligand)
    except RuntimeError as error:
        logging.warning("    * Warning: The desalting procedure has failed...")
        base_ligand["status_sub"].append(
            ["desalt", {"state": "failed", "text": "desalting failed"}]
        )
        if ctx["config"]["desalting_obligatory"] == "true":
            completion_event["seconds"] = time.perf_counter() - start_time
            return completion_event
        else:
            logging.warning(
                "    * Warning: Ligand will be further processed without desalting"
            )

    base_ligand["timers"].append(["desalt", time.perf_counter() - step_timer_start])

    # Neutralization

    step_timer_start = time.perf_counter()
    try:
        neutralization(ctx, base_ligand)
    except RuntimeError as error:
        logging.error("    * Warning: The neutralization has failed.")
        base_ligand["status_sub"].append(
            ["neutralization", {"state": "failed", "text": f"Failed {str(error)}"}]
        )

        # Can we go on?
        if ctx["config"]["neutralization_obligatory"] == "true":
            logging.error(
                "    * Warning: Ligand will be skipped since a successful neutralization is required according to the controlfile."
            )
            completion_event["seconds"] = time.perf_counter() - start_time
            return completion_event

        logging.warning(
            "    * Warning: Ligand will be further processed without neutralization"
        )
        base_ligand["smi_neutralized"] = base_ligand["smi_desalted"]

    base_ligand["timers"].append(
        ["neutralization", time.perf_counter() - step_timer_start]
    )

    # Stereoisomer generation

    step_timer_start = time.perf_counter()
    try:
        # base_ligand['smi_neutralized']
        stereoisomer_generation(ctx, base_ligand)
    except RuntimeError as error:
        logging.error(
            f"    * Warning: The stereoisomer generation has failed. (error: {str(error)}, smi: {str(base_ligand['smi_neutralized'])}"
        )
        base_ligand["status_sub"].append(
            ["stereoisomer", {"state": "failed", "text": f"Failed {str(error)}"}]
        )

        # Can we go on?
        if ctx["config"]["stereoisomer_obligatory"] == "true":
            logging.warning(
                "    * Warning: Ligand will be skipped since a successful stereoisomer generation is required according to the controlfile."
            )
            completion_event["seconds"] = time.perf_counter() - start_time
            return completion_event
        else:
            base_ligand["stereoisomer_smiles"] = [base_ligand["smi_neutralized"]]

    base_ligand["timers"].append(
        ["stereoisomer", time.perf_counter() - step_timer_start]
    )

    # In some cases this will generate additional steroisomers, which we need to process -- in other cases it will just be
    # the same ligand
    number_of_stereoisomers = len(base_ligand["stereoisomer_smiles"])

    # Loop through every stereoisomer
    for index, stereoisomer_smile_full in enumerate(base_ligand["stereoisomer_smiles"]):
        stereoisomer_timer_start = time.perf_counter()

        print(
            f"      * Processing Stereoisomer ({index+1} of {number_of_stereoisomers}) for {base_ligand['key']}"
        )

        # If the SMILES string has a space, take the first part
        stereoisomer_smile = stereoisomer_smile_full.split()[0]
        stereoisomer_key = f"{base_ligand['key']}_S{index}"
        logging.debug(
            f"processing stereoisomer {index}, stereoisomer_key:{stereoisomer_key} smile:{stereoisomer_smile}"
        )

        stereoisomer = {
            "key": stereoisomer_key,
            "smi": stereoisomer_smile,
            "remarks": base_ligand["remarks"].copy(),
            "status": "failed",
            "status_sub": [],
            "index": index,
            "seconds": 0,
            "timers": [],
        }

        #
        completion_event["stereoisomers"][stereoisomer_key] = stereoisomer

        try:
            process_stereoisomer(
                ctx, base_ligand, stereoisomer, completion_event["ligands"]
            )
        except RuntimeError as error:
            logging.error(f"Failed processing {stereoisomer_key} (error: {str(error)})")

        stereoisomer_timer_end = time.perf_counter()
        stereoisomer["seconds"] = stereoisomer_timer_end - stereoisomer_timer_start

    end_time = time.perf_counter()
    completion_event["seconds"] = end_time - start_time
    completion_event["status"] = "success"

    # captured and processed later
    return completion_event
