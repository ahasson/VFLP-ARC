#######################
# Step 2: Neutralization
import re
import time
import logging

from run.chemaxon import run_chemaxon_neutralization_standardizer
from run.obabel import run_obabel_neutralization


def run_neutralization_instance(ctx, ligand, program):
    step_timer_start = time.perf_counter()
    valid_programs = ("standardizer", "obabel")

    if program == "none":
        raise RuntimeError(f"No neutralization program remaining")
    elif program not in valid_programs:
        raise RuntimeError(f"Neutralization program '{program}' is not valid")

    try:
        if program == "standardizer":
            run_chemaxon_neutralization_standardizer(ctx, ligand)
        elif program == "obabel":
            run_obabel_neutralization(ctx, ligand)
    except RuntimeError as error:
        ligand["timers"].append(
            [f"{program}_neutralize", time.perf_counter() - step_timer_start]
        )
        return error

    ligand["timers"].append(
        [f"{program}_neutralize", time.perf_counter() - step_timer_start]
    )


def run_neutralization_generation(ctx, ligand):
    try:
        run_neutralization_instance(
            ctx, ligand, ctx["config"]["neutralization_program_1"]
        )
    except RuntimeError as error:
        run_neutralization_instance(
            ctx, ligand, ctx["config"]["neutralization_program_2"]
        )


def neutralization(ctx, ligand):

    if ctx["config"]["neutralization"] == "true":
        run = 0

        if ctx["config"]["neutralization_mode"] == "always":
            run = 1
        elif (
            ctx["config"]["neutralization_mode"] == "only_genuine_desalting"
            and ligand["number_of_fragments"] > 1
        ):
            run = 1
        elif (
            ctx["config"]["neutralization_mode"]
            == "only_genuine_desalting_and_if_charged"
            and ligand["number_of_fragments"] > 1
        ):
            match = re.search(r"(?P<charge>(\-\]|\+\]))", ligand["smi_desalted"])
            if match:
                run = 1

        if run:
            try:
                run_neutralization_generation(ctx, ligand)
            except RuntimeError as error:
                raise

            ligand["status_sub"].append(
                ["neutralization", {"state": "success", "text": "genuine"}]
            )
        else:
            # Skipping neutralization
            logging.debug(
                "* This ligand does not need to be neutralized, leaving it untouched."
            )
            ligand["status_sub"].append(
                ["neutralization", {"state": "success", "text": "untouched"}]
            )

            ligand["smi_neutralized"] = ligand["smi_desalted"]

    else:
        # Skipping neutralization
        ligand["smi_neutralized"] = ligand["smi_desalted"]
