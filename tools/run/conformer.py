#######################
# Step 6: 3D conformation generation

import time

from run.chemaxon import chemaxon_conformation
from run.obabel import obabel_conformation


def run_conformation_generation(ctx, tautomer, output_file):

    try:
        run_conformation_instance(
            ctx, tautomer, ctx["config"]["conformation_program_1"], output_file
        )
    except RuntimeError as error:
        run_conformation_instance(
            ctx, tautomer, ctx["config"]["conformation_program_2"], output_file
        )


def run_conformation_instance(ctx, tautomer, program, output_file):
    step_timer_start = time.perf_counter()
    valid_programs = ("molconvert", "obabel")

    if program == "none":
        raise RuntimeError(f"No conformation program remaining")
    elif program not in valid_programs:
        raise RuntimeError(f"Conformation program '{program}' is not valid")

    try:
        if program == "molconvert":
            chemaxon_conformation(ctx, tautomer, output_file)
        elif program == "obabel":
            obabel_conformation(ctx, tautomer, output_file)
    except RuntimeError as error:
        tautomer["timers"].append(
            [f"{program}_conformation", time.perf_counter() - step_timer_start]
        )
        raise error

    tautomer["timers"].append(
        [f"{program}_conformation", time.perf_counter() - step_timer_start]
    )
