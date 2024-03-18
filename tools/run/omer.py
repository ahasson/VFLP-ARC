import time
import logging
import shutil
from pathlib import Path
import selfies

from run.chemaxon import (
    run_chemaxon_stereoisomer_generation,
    run_chemaxon_tautomer_generation,
)
from run.utils import write_file_single
from run.obabel import (
    run_obabel_tautomerization,
    obabel_generate_pdb,
    obabel_check_energy,
    obabel_generate_targetformat,
)
from run.rdkit import run_rdkit_stereoisomer_generation
from run.conformer import run_conformation_generation
from run.protonation import run_protonation_generation
from run.utils import get_intermediate_dir_ligand
from run.tranche import generate_attributes, tranche_assignment


def get_intermediate_dir_tautomer(ctx, tautomer):

    output_dir: Path = get_intermediate_dir_ligand(ctx) / tautomer["key"]
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir


def process_stereoisomer(ctx, ligand, stereoisomer, completion_ligands):

    stereoisomer["status"] = "failed"

    # Tautomer generation
    step_timer_start = time.perf_counter()
    try:
        tautomer_generation(ctx, stereoisomer)
    except RuntimeError as error:
        logging.error(
            f"    * Warning: The tautomerization has failed. (error: {str(error)}, smi: {str(stereoisomer['smi'])}"
        )
        stereoisomer["status_sub"].append(
            ["tautomerization", {"state": "failed", "text": f"Failed {str(error)}"}]
        )

        # Can we go on?
        if ctx["config"]["tautomerization_obligatory"] == "true":
            logging.warning(
                "    * Warning: stereoisomer will be skipped since a successful tautomerization is required according to the controlfile."
            )
            stereoisomer["timers"].append(
                ["tautomerization", time.perf_counter() - step_timer_start]
            )
            return
        else:
            stereoisomer["tautomer_smiles"] = [stereoisomer["smi"]]

    stereoisomer["timers"].append(
        ["tautomerization", time.perf_counter() - step_timer_start]
    )

    # In some cases this will generate additional tautomers, which we need to process -- in other cases it will just be
    # the same ligand
    number_of_tautomers = len(stereoisomer["tautomer_smiles"])

    # Loop through every tautomer
    for index, tautomer_smile_full in enumerate(stereoisomer["tautomer_smiles"]):
        tautomer_timer_start = time.perf_counter()

        print(
            f"      ** Processing tautomer ({index+1} of {number_of_tautomers}) for {stereoisomer['key']}"
        )

        # If the SMILES string has a space, take the first part
        tautomer_smile = tautomer_smile_full.split()[0]
        tautomer_key = f"{stereoisomer['key']}_T{index}"
        logging.debug(
            f"processing {index}, tautomer_key:{tautomer_key} smile:{tautomer_smile}"
        )

        tautomer = {
            "key": tautomer_key,
            "smi": tautomer_smile,
            "smi_stereoisomer": stereoisomer["smi"],
            "smi_original": ligand["smi"],
            "remarks": ligand["remarks"].copy(),
            "status": "failed",
            "status_sub": [],
            "index": index,
            "seconds": 0,
            "timers": [],
        }

        #
        completion_ligands[tautomer_key] = tautomer

        try:
            process_tautomer(ctx, ligand, completion_ligands[tautomer_key])
        except RuntimeError as error:
            logging.error(f"Failed processing {tautomer_key} (error: {str(error)})")

        tautomer_timer_end = time.perf_counter()
        tautomer["seconds"] = tautomer_timer_end - tautomer_timer_start

    stereoisomer["status"] = "success"


#######################
# Step 3: Stereoisomer Generation


def stereoisomer_generation_instance(ctx, ligand, program):
    step_timer_start = time.perf_counter()
    valid_programs = ("rdkit", "cxcalc")

    if program == "none":
        raise RuntimeError(f"No stereoisomer generation program remaining")
    elif program not in valid_programs:
        raise RuntimeError(f"Stereoisomer generation program '{program}' is not valid")

    try:
        if program == "cxcalc":
            logging.info("Starting the stereoisomer generation with cxcalc")
            run_chemaxon_stereoisomer_generation(ctx, ligand)
            ligand["status_sub"].append(
                ["stereoisomer", {"state": "success", "text": ""}]
            )
        elif program == "rdkit":
            logging.info("Starting the stereoisomer generation with rdkit")
            # Need to include the option for only one (True option) or multiple stereoisomers (currently multiple)
            run_rdkit_stereoisomer_generation(ctx, ligand, False)
            ligand["status_sub"].append(
                ["stereoisomer", {"state": "success", "text": ""}]
            )
    except RuntimeError as error:
        ligand["timers"].append(
            [
                f"{program}_stereoisomer_generation",
                time.perf_counter() - step_timer_start,
            ]
        )
        return error

    ligand["timers"].append(
        [f"{program}_stereoisomer_generation", time.perf_counter() - step_timer_start]
    )


def stereoisomer_generation(ctx, ligand):

    if ctx["config"]["stereoisomer_generation"] == "true":
        try:
            stereoisomer_generation_instance(
                ctx, ligand, ctx["config"]["stereoisomer_generation_program_1"]
            )
        except RuntimeError as error:
            stereoisomer_generation_instance(
                ctx, ligand, ctx["config"]["stereoisomer_generation_program_2"]
            )
    else:
        ligand["stereoisomer_smiles"] = [ligand["smi_neutralized"]]


#######################
# Step 3: Tautomerization


def tautomerization_instance(ctx, stereoisomer, program):
    step_timer_start = time.perf_counter()
    valid_programs = ("cxcalc", "obabel")

    if program == "none":
        raise RuntimeError(f"No tautomerization program remaining")
    elif program not in valid_programs:
        raise RuntimeError(f"Tautomerization program '{program}' is not valid")

    try:
        if program == "cxcalc":
            logging.info("Starting the tautomerization with cxcalc")
            run_chemaxon_tautomer_generation(ctx, stereoisomer)
            stereoisomer["status_sub"].append(
                ["tautomerization", {"state": "success", "text": ""}]
            )
        elif program == "obabel":
            logging.info("Starting the tautomerization with obtautomer")
            run_obabel_tautomerization(ctx, stereoisomer)
            stereoisomer["status_sub"].append(
                ["tautomerization", {"state": "success", "text": ""}]
            )
    except RuntimeError as error:
        stereoisomer["timers"].append(
            [f"{program}_tautomerize", time.perf_counter() - step_timer_start]
        )
        return error

    stereoisomer["timers"].append(
        [f"{program}_tautomerize", time.perf_counter() - step_timer_start]
    )


def tautomer_generation(ctx, stereoisomer):

    if ctx["config"]["tautomerization"] == "true":
        try:
            tautomerization_instance(
                ctx, stereoisomer, ctx["config"]["tautomerization_program_1"]
            )
        except RuntimeError as error:
            tautomerization_instance(
                ctx, stereoisomer, ctx["config"]["tautomerization_program_2"]
            )
    else:
        stereoisomer["tautomer_smiles"] = [stereoisomer["smi"]]


#### For each Tautomer -- Steps 4-9


def process_tautomer(ctx, ligand, tautomer):

    tautomer["status"] = "failed"
    tautomer["smi_protomer"] = tautomer["smi"]
    tautomer["intermediate_dir"] = get_intermediate_dir_tautomer(ctx, tautomer)

    if ctx["config"]["protonation_state_generation"] == "true":
        try:
            run_protonation_generation(ctx, tautomer)
        except RuntimeError as error:
            tautomer["status_sub"].append(
                ["protonation", {"state": "failed", "text": f"{str(error)}"}]
            )

            logging.warning("* Warning: Both protonation attempts have failed.")

            if ctx["config"]["protonation_obligatory"] == "true":
                logging.error(
                    "* Warning: Ligand will be skipped since a successful protonation is required according to the controlfile."
                )
                raise
            else:
                logging.error(
                    "* Warning: Ligand will be further processed without protonation, which might result in unphysiological protonation states."
                )
                tautomer["remarks"][
                    "protonation"
                ] = "WARNING: Molecule was not protonated at physiological pH (protonation with both obabel and cxcalc has failed)"
        else:
            tautomer["status_sub"].append(
                ["protonation", {"state": "success", "text": ""}]
            )

    # Update remarks used in output files

    tautomer["remarks"]["basic"] = "Small molecule (ligand)"
    tautomer["remarks"]["compound"] = f"Compound: {tautomer['key']}"

    tautomer["remarks"]["smiles_original"] = f"SMILES_orig: {tautomer['smi_original']}"
    tautomer["remarks"][
        "smiles_current"
    ] = f"SMILES_current: {tautomer['smi_protomer']}"

    # Generate data for the attributes as needed
    step_timer_start = time.perf_counter()
    tautomer["attr"] = generate_attributes(ctx, ligand, tautomer)
    tautomer["timers"].append(
        ["attr-generation", time.perf_counter() - step_timer_start]
    )

    # If there are specific attributes to place in remarks, do that now
    tautomer["remarks"]["additional_attr"] = []
    for attribute_type in ctx["config"]["attributes_to_generate"]:
        tautomer["remarks"]["additional_attr"].append(
            f"{attribute_type}: {tautomer['attr'][attribute_type]}"
        )

    # Assign the tranches if needed
    step_timer_start = time.perf_counter()
    if ctx["config"]["tranche_assignments"] == "true":
        try:
            tranche_assignment(ctx, ligand, tautomer)
        except RuntimeError as error:
            tautomer["status_sub"].append(
                ["tranche-assignment", {"state": "failed", "text": f"{str(error)}"}]
            )
            tautomer["timers"].append(
                ["tranche-assignment", time.perf_counter() - step_timer_start]
            )
            logging.error(f"tranche_assignment failed for {tautomer['key']}")
            logging.error(
                "* Error: The tranche assignments have failed, ligand will be skipped."
            )
            raise RuntimeError(
                "The tranche assignments have failed, ligand will be skipped"
            ) from error

    tautomer["status_sub"].append(
        ["tranche-assignment", {"state": "success", "text": ""}]
    )
    tautomer["timers"].append(
        ["tranche-assignment", time.perf_counter() - step_timer_start]
    )

    # If SELFIES are requested, generate them here since we will place them in the
    # pdb file as well

    if "selfies" in ctx["config"]["target_formats"]:
        try:
            tautomer["status_sub"].append(["selfies", {"state": "success", "text": ""}])
            tautomer["selfies"] = selfies.encoder(tautomer["smi_protomer"])
            tautomer["remarks"]["additional_attr"].append(
                f"selfies: {tautomer['selfies']}"
            )
        except selfies.exceptions.EncoderError as error:
            tautomer["status_sub"].append(
                ["selfies", {"state": "failed", "text": f"{str(error)}"}]
            )
            raise RuntimeError("Selfies generation failed")

    # 3D conformation generation

    # Where to place any PDB output (either from conformation or generation step)
    tautomer["pdb_file"] = str(tautomer["intermediate_dir"] / f"gen.pdb")

    conformation_success = "true"
    if ctx["config"]["conformation_generation"] == "true":
        step_timer_start = time.perf_counter()

        try:
            run_conformation_generation(ctx, tautomer, tautomer["pdb_file"])
        except RuntimeError as error:
            tautomer["status_sub"].append(
                ["conformation", {"state": "failed", "text": f"{str(error)}"}]
            )
            conformation_success = "false"
            logging.error(f"conformation_generation failed for {tautomer['key']}")
            if ctx["config"]["conformation_obligatory"] == "true":
                tautomer["timers"].append(
                    ["conformation", time.perf_counter() - step_timer_start]
                )
                raise RuntimeError(
                    f"Conformation failed, but is required error:{str(error)}"
                ) from error
        else:
            tautomer["status_sub"].append(
                ["conformation", {"state": conformation_success, "text": ""}]
            )

        tautomer["timers"].append(
            ["conformation", time.perf_counter() - step_timer_start]
        )

    ## PDB Generation
    # If conformation generation failed, and we reached this point,
    # then conformation_obligatory=false, so we do not need to check this

    if (
        ctx["config"]["conformation_generation"] == "false"
        or conformation_success == "false"
    ):
        try:
            obabel_generate_pdb(ctx, tautomer, tautomer["pdb_file"])
        except RuntimeError as error:
            tautomer["status_sub"].append(
                ["pdb-generation", {"state": "failed", "text": f"{str(error)}"}]
            )
            logging.warning(
                "    * Warning: Ligand will be skipped since a successful PDB generation is mandatory."
            )
            logging.error(f"obabel_generate_pdb failed for {tautomer['key']}")
            raise RuntimeError("PDB generation failed, but is required") from error

    # At this point a valid pdb file will be at tautomer['pdb_file']

    # Checking the potential energy
    if ctx["config"]["energy_check"] == "true":
        logging.warning("\n * Starting to check the potential energy of the ligand")

        if not obabel_check_energy(
            ctx, tautomer, tautomer["pdb_file"], ctx["config"]["energy_max"]
        ):
            logging.warning(
                "    * Warning: Ligand will be skipped since it did not pass the energy-check."
            )
            tautomer["status_sub"].append(
                ["energy-check", {"state": "failed", "text": ""}]
            )
            raise RuntimeError("Failed energy check")
        else:
            tautomer["status_sub"].append(
                ["energy-check", {"state": "success", "text": ""}]
            )

    ## Target formats

    logging.debug(f"target formats is {ctx['config']['target_formats']}")
    step_timer_start = time.perf_counter()
    for target_format in ctx["config"]["target_formats"]:

        logging.debug(f"target_formats is {target_format}")

        # Put in the main job temporary directory, not the temp dir
        # for this ligand

        output_file_parts = [
            ctx["collection_temp_dir"].name,
            "complete",
            target_format,
            ctx["metatranche"],
            ctx["tranche"],
            ctx["collection_name"],
        ]

        output_dir = Path("/".join(output_file_parts))
        output_dir.mkdir(parents=True, exist_ok=True)

        if ctx["config"]["tranche_assignments"] == "true":
            output_file = (
                output_dir
                / f"{tautomer['tranche_string']}_{tautomer['key']}.{target_format}"
            )
        else:
            output_file = output_dir / f"{tautomer['key']}.{target_format}"

        try:
            generate_target_format(
                ctx, tautomer, target_format, tautomer["pdb_file"], output_file
            )
        except RuntimeError as error:
            logging.error(f"failed generation for format {target_format}")
            tautomer["status_sub"].append(
                [
                    f"targetformat-generation({target_format})",
                    {"state": "failed", "text": str(error)},
                ]
            )
        else:
            logging.debug(f"succeeded generation for format {target_format}")
            tautomer["status_sub"].append(
                [
                    f"targetformat-generation({target_format})",
                    {"state": "success", "text": ""},
                ]
            )

    tautomer["timers"].append(["targetformats", time.perf_counter() - step_timer_start])

    # Mark the complete tautomer as successfully finished if we get here
    # Note that the target format generation could have failed at this point

    tautomer["status"] = "success"


def generate_target_format(ctx, tautomer, target_format, pdb_file, output_file):

    if target_format == "smi":
        # We can just use the SMI string that we already have
        write_file_single(output_file, tautomer["smi_protomer"])
    elif target_format == "pdb":
        # The input to this function is already a pdb file, so we can
        # just copy it
        shutil.copyfile(pdb_file, output_file)
    elif target_format == "selfies":
        # We can just use the SELFIES string that we already have
        write_file_single(output_file, tautomer["selfies"])
    else:
        obabel_generate_targetformat(
            ctx, tautomer, target_format, pdb_file, output_file
        )
