#######################
# Step 4: Protonation
import time
from run.chemaxon import cxcalc_protonate
from run.obabel import run_obabel_protonation


def run_protonation_instance(ctx, tautomer, program):
    step_timer_start = time.perf_counter()
    valid_programs = ("cxcalc", "obabel")

    if program == "none":
        raise RuntimeError("No protonation program remaining")
    if program not in valid_programs:
        raise RuntimeError(f"Protonation program '{program}' is not valid")

    try:
        if program == "cxcalc":
            cxcalc_protonate(ctx, tautomer)
        elif program == "obabel":
            run_obabel_protonation(ctx, tautomer)
    except RuntimeError as error:
        tautomer["timers"].append(
            [f"{program}_protonate", time.perf_counter() - step_timer_start]
        )
        return error

    tautomer["timers"].append(
        [f"{program}_protonate", time.perf_counter() - step_timer_start]
    )


def run_protonation_generation(ctx, tautomer):
    try:
        run_protonation_instance(ctx, tautomer, ctx["config"]["protonation_program_1"])
    except RuntimeError as _error:
        run_protonation_instance(ctx, tautomer, ctx["config"]["protonation_program_2"])
