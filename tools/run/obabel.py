################
# OBabel Components
#

# General functions for obabel and obtautomer
import subprocess
import re
import os
import time
from datetime import datetime
import logging

from rdkit import Chem

from run.utils import write_file_single

# Sometimes the coordinates do not end up as 3D. In order to
# verify that they are not completely, flat take the coordinates and
# remove '0 . + -' from the lines and then see if there
# is anything left other than space. If we find at least one
# line with non-zero coordinates that is sufficient (VERIFY)


def nonzero_pdb_coordinates(output_file):

    with open(output_file, "r") as read_file:
        for line in read_file:
            line = line.strip()

            if re.search(r"^(ATOM|HETATM)", line):
                line_parts = line.split()

                components_to_check = "".join(line_parts[5:7])
                components_to_check = re.sub(r"[0.+-]", "", components_to_check)
                if not re.search(r"^\s*$", line):
                    return True

    return False


def generate_remarks(remark_list, remark_order="default", target_format="pdb"):

    if remark_order == "default":
        remark_ordering = [
            "basic",
            "compound",
            "smiles_original",
            "smiles_current",
            "desalting",
            "neutralization",
            "stereoisomer",
            "tautomerization",
            "protonation",
            "generation",
            "conformation",
            "targetformat",
            "trancheassignment",
            "trancheassignment_attr",
            "additional_attr",
            "tranche_str",
            "date",
            "collection_key",
        ]
    else:
        remark_ordering = remark_order

    if target_format == "mol2":
        remark_prefix = "# "
    elif target_format == "pdb" or target_format == "pdbqt":
        remark_prefix = "REMARK    "
    else:
        raise RuntimeError(
            f"Invalid target_format type ({target_format})for generating remarks"
        )

    remark_parts = []
    for remark in remark_ordering:
        if remark in remark_list:
            if isinstance(remark_list[remark], list):
                for sub_remark in remark_list[remark]:
                    if sub_remark != "":
                        remark_parts.append(remark_prefix + " * " + sub_remark)
            elif remark_list[remark] != "":
                remark_parts.append(remark_prefix + remark_list[remark])
    remark_string = "\n".join(remark_parts)

    return remark_string


def debug_save_output(ctx, file, stdout="", stderr="", tautomer=None):
    if ctx["store_all_intermediate_logs"] == "true":
        if tautomer is not None:
            save_output(tautomer["intermediate_dir"] / file, stdout, stderr)
        else:
            save_output(ctx["intermediate_dir"] / file, stdout, stderr)


def save_output(save_logfile, save_stdout, save_stderr):
    with open(save_logfile, "wt", encoding="utf8") as write_file:
        write_file.write("STDOUT ---------\n")
        write_file.write(save_stdout)
        write_file.write("\nSTDERR ---------\n")
        write_file.write(save_stderr)


def file_is_empty(filename):
    with open(filename, "rt", encoding="utf8") as read_file:
        for _ in read_file:
            return False
    return True


def run_obabel_general(obabelargs, timeout=30, save_logfile=""):

    obabelargs_x = []
    for arg in obabelargs:
        if arg != "":
            obabelargs_x.append(arg)

    cmd = ["obabel", *obabelargs]

    try:
        ret = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as err:
        raise RuntimeError(f"obabel timed out") from err

    if ret.returncode != 0:
        raise RuntimeError(f"Return code from obabel is {ret.returncode}")

    for line in ret.stdout.splitlines() + ret.stderr.splitlines():
        if re.search(r"failed|timelimit|error|no such file|not found", line):
            raise RuntimeError(
                f"An error flag was detected in the log files from obabel"
            )

    return {"stderr": ret.stderr, "stdout": ret.stdout}


def run_obtautomer_general(obtautomerargs, output_file, timeout=30, save_logfile=""):

    obtautomerargs_x = []
    for arg in obtautomerargs:
        if arg != "":
            obtautomerargs_x.append(arg)

    cmd = ["obtautomer", *obtautomerargs]

    try:
        ret = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        write_file_single(output_file, ret.stdout)
    except subprocess.TimeoutExpired as err:
        raise RuntimeError(f"obtautomer timed out") from err

    if ret.returncode != 0:
        raise RuntimeError(f"Return code from obtautomer is {ret.returncode}")

    for line in ret.stdout.splitlines() + ret.stderr.splitlines():
        if re.search(r"failed|timelimit|error|no such file|not found", line):
            raise RuntimeError(
                f"An error flag was detected in the log files from obtautomer"
            )

    return {"stderr": ret.stderr, "stdout": ret.stdout}


def run_obabel_general_get_value(obabelargs, output_file, timeout=30):
    ret = run_obabel_general(obabelargs, timeout=timeout)

    if not os.path.isfile(output_file):
        raise RuntimeError(f"No output file from obabel")

    with open(output_file, "r") as read_file:
        lines = read_file.readlines()

    if len(lines) == 0:
        raise RuntimeError(f"No output in file from obabel")

    return lines[0].replace("\t", "").replace("\n", "").strip()


def run_obtautomer_general_get_value(obtautomerargs, output_file, timeout=30):
    ret = run_obtautomer_general(obtautomerargs, output_file, timeout=timeout)

    if not os.path.isfile(output_file):
        raise RuntimeError(f"No output file from obtautomer")

    with open(output_file, "r") as read_file:
        lines = read_file.readlines()

    if len(lines) == 0:
        raise RuntimeError(f"No output in file from obtautomer")

    return lines


# Step 1: Desalt (not performed by OBabel)
# Step 2: Neutralization by OBabel


def run_obabel_neutralization(ctx, ligand):

    input_file = f"{ctx['temp_dir'].name}/obabel.neutra_input.smi"
    output_file = f"{ctx['temp_dir'].name}/obabel.neutra_output.smi"

    # Output the SMI to a temp input file
    write_file_single(input_file, ligand["smi_desalted"])

    cmd = ["-ismi", input_file, "--neutralize", "-osmi", "-O", output_file]

    ligand["smi_neutralized"] = run_obabel_general_get_value(
        cmd, output_file, timeout=ctx["config"]["obabel_neutralization_timeout"]
    )
    ligand["neutralization_success"] = 1
    ligand["remarks"]["neutralization"] = "The compound was neutralized by Open Babel."
    ligand["neutralization_type"] = "genuine"


# Step 3a: Stereoisomer Generation by RDKit


def perform_isomer_unique_correction(ctx, sterio_smiles_ls):
    """
    When RDkit enumerates sterio-isomers, some of them can be redundant.
    To resolve this, we use obabel to convert them into 3D, convert back to
    canonical smiles, and, use this as the final list of stereoisomer smiles
    for a particular molecule.

    Parameters
    ----------
    sterio_smiles_ls : list of strings
        A list of valid SMILES strings.

    Returns
    -------
    isomers_canon : list of strings
        A list of valid SMILES strings, containing the corrected smiles!.
    """

    smi_file = ctx["intermediate_dir"] / "sterio.smi"
    sdf_file = ctx["intermediate_dir"] / "sterio.sdf"

    isomers_canon = []

    for smi in sterio_smiles_ls:
        with open(smi_file, "w") as f:
            f.writelines(smi)

        try:
            subprocess.run(
                ["obabel", smi_file, "--gen3D", "-O", sdf_file],
                capture_output=True,
                text=True,
                timeout=float(ctx["config"]["obabel_stereoisomer_timeout"]),
            )
            subprocess.run(
                ["rm", smi_file],
                capture_output=True,
                text=True,
                timeout=float(ctx["config"]["obabel_stereoisomer_timeout"]),
            )
            subprocess.run(
                ["obabel", sdf_file, "-O", smi_file],
                capture_output=True,
                text=True,
                timeout=float(ctx["config"]["obabel_stereoisomer_timeout"]),
            )

            with open(smi_file, "r") as f:
                new_smi = f.readlines()
            new_smi = new_smi[0].strip()

            subprocess.run(
                ["rm", smi_file],
                capture_output=True,
                text=True,
                timeout=float(ctx["config"]["obabel_stereoisomer_timeout"]),
            )
            subprocess.run(
                ["rm", sdf_file],
                capture_output=True,
                text=True,
                timeout=float(ctx["config"]["obabel_stereoisomer_timeout"]),
            )

            mol = Chem.MolFromSmiles(new_smi)
            new_smi_canon = Chem.MolToSmiles(mol, canonical=True)

            isomers_canon.append(new_smi_canon)

        except Exception:
            continue

    return list(set(isomers_canon))


# Step 3b: Tautomerization by OBabel


def run_obabel_tautomerization(ctx, stereoisomer):

    input_file = f"{ctx['temp_dir'].name}/obabel.tauto_input.smi"
    output_file = f"{ctx['temp_dir'].name}/obabel.tauto_output.smi"

    # Output the SMI to a temp input file
    write_file_single(input_file, stereoisomer["smi"])

    cmd = [input_file]

    stereoisomer["tautomer_smiles"] = []
    stereoisomer["tautomer_smiles"] = run_obtautomer_general_get_value(
        cmd, output_file, timeout=ctx["config"]["obabel_tautomerization_timeout"]
    )
    del stereoisomer["tautomer_smiles"][-1]
    stereoisomer["tautomer_smiles"] = [
        i.strip() for i in stereoisomer["tautomer_smiles"]
    ]

    stereoisomer["remarks"][
        "tautomerization"
    ] = "The tautomeric state was generated by Open Babel."


# Step 4: Protonation

# Step 4: 3D Protonation


def run_obabel_protonation(ctx, tautomer):

    input_file = f"{ctx['temp_dir'].name}/obabel.proto.{tautomer['key']}_input.smi"
    output_file = f"{ctx['temp_dir'].name}/obabel.proto.{tautomer['key']}_output.smi"

    # Output the SMI to a temp input file
    write_file_single(input_file, tautomer["smi"])

    cmd = [
        "-p",
        ctx["config"]["protonation_pH_value"],
        "-ismi",
        input_file,
        "-osmi",
        "-O",
        output_file,
    ]

    tautomer["smi_protomer"] = run_obabel_general_get_value(
        cmd, output_file, timeout=ctx["config"]["obabel_protonation_timeout"]
    )


# Step 5: Assign Tranche


def run_obabel_attributes(ctx, tautomer, local_file, attributes):

    step_timer_start = time.perf_counter()

    captured_attrs = {}
    obabel_attrs = []
    for attr in attributes:
        if attributes[attr]["prog"] == "obabel":
            obabel_attrs.append(attributes[attr]["prog_name"])

    cmd = ["obprop", local_file]
    try:
        ret = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except subprocess.TimeoutExpired as err:
        raise RuntimeError(f"obprop timed out") from err

    output_lines = ret.stdout.splitlines()

    if len(output_lines) == 0:
        raise RuntimeError(f"No output file from oprop")
    else:
        for line in output_lines:
            line_values = line.split(maxsplit=1)

            if len(line_values) == 2:
                line_key, line_value = line_values

                if line_key in obabel_attrs:
                    captured_attrs[line_key] = line_value
                    logging.debug(f"Got {line_key}:{line_value} for obabel")

    for attr in attributes:
        if attributes[attr]["prog"] == "obabel":
            if attributes[attr]["prog_name"] in captured_attrs:
                attributes[attr]["val"] = captured_attrs[attributes[attr]["prog_name"]]

    tautomer["timers"].append(
        ["obabel_attributes", time.perf_counter() - step_timer_start]
    )


def run_obabel_hbX(local_file, append_val):

    return_val = "INVALID"
    cmd = ["-ismi", local_file, "-osmi", "--append", append_val]

    ret = run_obabel_general(cmd)
    output_lines = ret["stdout"].splitlines()

    if len(output_lines) == 0:
        raise RuntimeError(f"No output file from obabel")
    else:
        line_values = output_lines[0].split()
        if len(line_values) >= 2:
            # line_key, line_value = line_values
            # return the last value -- sometimes there are two strings
            return line_values[-1]

    raise RuntimeError(f"Unable to parse attribute from obabel")


def run_obabel_hba(local_file, tautomer):
    step_timer_start = time.perf_counter()
    attr_val = run_obabel_hbX(local_file, "HBA1")
    tautomer["timers"].append(
        ["obabel_attr_hba", time.perf_counter() - step_timer_start]
    )
    return attr_val


def run_obabel_hbd(local_file, tautomer):
    step_timer_start = time.perf_counter()
    attr_val = run_obabel_hbX(local_file, "HBD")
    tautomer["timers"].append(
        ["obabel_attr_hbd", time.perf_counter() - step_timer_start]
    )
    return attr_val


#
# Step 6: 3D conformation generation
def obabel_generate_pdb_general(
    ctx, tautomer, output_file, conformation, must_have_output=1, timeout=30
):

    input_file = f"{ctx['temp_dir'].name}/obabel.conf.{tautomer['key']}_input.smi"
    write_file_single(input_file, tautomer["smi_protomer"])

    output_file_tmp = f"{output_file}.tmp"
    if conformation:
        cmd = ["--gen3d", "-ismi", input_file, "-opdb", "-O", output_file_tmp]
        logging.debug(f"Running obabel conformation on smi:'{tautomer['smi_protomer']}")
    else:
        cmd = ["-ismi", input_file, "-opdb", "-O", output_file_tmp]

    ret = run_obabel_general(cmd, timeout=timeout)
    output_lines = ret["stdout"].splitlines()

    if must_have_output and len(output_lines) == 0:
        raise RuntimeError(f"No output")

    if not os.path.isfile(output_file_tmp):
        raise RuntimeError(f"No PDB file generated")
    if file_is_empty(output_file_tmp):
        raise RuntimeError(
            f"obabel outputfile is empty (smi: {tautomer['smi_protomer']}"
        )
    if not nonzero_pdb_coordinates(output_file_tmp):
        raise RuntimeError(
            f"The output PDB file exists but does not contain valid coordinates."
        )

    # We are successful
    if conformation:
        tautomer["remarks"][
            "conformation"
        ] = f"Generation of the 3D conformation was carried out by Open Babel version {ctx['versions']['obabel']}"
        tautomer["remarks"][
            "targetformat"
        ] = "Format generated as part of conformation."
    else:
        tautomer["remarks"][
            "generation"
        ] = f"Generation of the the PDB file (without conformation generation) was carried out by Open Babel version {ctx['versions']['obabel']}"
        tautomer["remarks"]["targetformat"] = "Format generated as part of generation."

    first_smile_component = tautomer["smi_protomer"].split()[0]
    # tautomer['remarks']['smiles'] = f"SMILES: {first_smile_component}"

    local_remarks = tautomer["remarks"].copy()
    local_remarks.pop("compound")
    remark_string = generate_remarks(local_remarks) + "\n"

    # Modify the output file as needed #
    with open(output_file, "w") as write_file:
        write_file.write(f"COMPND    Compound: {tautomer['key']}\n")
        write_file.write(remark_string)

        with open(output_file_tmp, "r") as read_file:

            for line in read_file:
                if re.search(r"COMPND|AUTHOR", line):
                    continue

                if re.search(r"^\s*$", line):
                    continue

                line = re.sub(r" UN[LK] ", " LIG ", line)
                write_file.write(line)


def obabel_conformation(ctx, tautomer, output_file):
    obabel_generate_pdb_general(
        ctx,
        tautomer,
        output_file,
        conformation=1,
        must_have_output=0,
        timeout=ctx["config"]["obabel_conformation_timeout"],
    )


# Step 7: PDB Generation


def obabel_generate_pdb(ctx, tautomer, output_file):
    step_timer_start = time.perf_counter()
    obabel_generate_pdb_general(
        ctx, tautomer, output_file, conformation=0, must_have_output=0
    )
    tautomer["timers"].append(
        ["obabel_generation", time.perf_counter() - step_timer_start]
    )


# Step 8: Energy Check
# Format for last line should be:
# 	obv2:	TOTAL ENERGY = 91.30741 kcal/mol
# 	obv3:	TOTAL ENERGY = 880.18131 kJ/mol
# VERIFY


def obabel_check_energy(ctx, tautomer, input_file, max_energy):

    step_timer_start = time.perf_counter()

    with open(input_file, "r") as read_file:
        lines = read_file.readlines()

    try:
        ret = subprocess.run(
            ["obenergy", input_file], capture_output=True, text=True, timeout=30
        )
    except subprocess.TimeoutExpired as err:
        tautomer["timers"].append(["obenergy", time.perf_counter() - step_timer_start])
        raise RuntimeError(f"obprop timed out") from err

    tautomer["timers"].append(["obenergy", time.perf_counter() - step_timer_start])

    debug_save_output(
        stdout=ret.stdout,
        stderr=ret.stderr,
        tautomer=tautomer,
        ctx=ctx,
        file="energy_check",
    )

    output_lines = ret.stdout.splitlines()

    kcal_to_kJ = 4.184

    if len(output_lines) > 0:
        # obabel v2
        match = re.search(
            r"^TOTAL\s+ENERGY\s+\=\s+(?P<energy>\d+\.?\d*)\s+(?P<energy_unit>(kcal|kJ))",
            output_lines[-1],
        )
        if match:
            energy_value = float(match.group("energy"))
            if match.group("energy_unit") == "kcal":
                energy_value *= kcal_to_kJ

            if energy_value <= float(max_energy):
                return 1

    tautomer["status_sub"].append(
        ["energy-check", {"state": "failed", "text": "|".join(output_lines)}]
    )

    return 0


# Step 9: Generate Target Formats


def obabel_generate_targetformat(
    ctx, tautomer, target_format, input_pdb_file, output_file
):

    step_timer_start = time.perf_counter()

    output_file_tmp = tautomer["intermediate_dir"] / f"tmp.{target_format}"

    cmd = ["-ipdb", input_pdb_file, "-O", output_file_tmp]

    ret = run_obabel_general(cmd)

    if not os.path.isfile(output_file_tmp):
        logging.debug("no output file generated")
        tautomer["timers"].append(
            [f"obabel_generate_{target_format}", time.perf_counter() - step_timer_start]
        )
        raise RuntimeError("No output file generated")

    if target_format == "pdb" or target_format == "pdbqt":
        if not nonzero_pdb_coordinates(output_file_tmp):
            logging.debug(
                "The output PDB(QT) file exists but does not contain valid coordinates"
            )
            tautomer["timers"].append(
                [
                    f"obabel_generate_{target_format}",
                    time.perf_counter() - step_timer_start,
                ]
            )
            raise RuntimeError(
                "The output PDB(QT) file exists but does not contain valid coordinates."
            )

    # Slurp in the temporary file

    lines = ()
    with open(output_file_tmp, "r") as read_file:
        lines = read_file.readlines()

    if len(lines) == 0:
        logging.debug("Output file is empty")
        tautomer["timers"].append(
            [f"obabel_generate_{target_format}", time.perf_counter() - step_timer_start]
        )
        raise RuntimeError("The output file is empty.")

    # Setup the remarks information
    now = datetime.now()

    remark_string = ""
    if target_format == "pdb":
        remarks = tautomer["remarks"].copy()
        remarks.pop("compound")
        remarks["targetformat"] = (
            f"Generation of the the target format file ({target_format}) was carried out by Open Babel version {ctx['versions']['obabel']}"
        )
        remarks["date"] = f'Created on {now.strftime("%Y-%m-%d %H:%M:%S")}'
        remark_string += f"COMPND    Compound: {tautomer['key']}\n"
        remark_string += (
            generate_remarks(
                remarks,
                remark_order=["targetformat", "date"],
                target_format=target_format,
            )
            + "\n"
        )
    elif target_format == "pdbqt":
        remarks = tautomer["remarks"].copy()
        remarks.pop("compound")
        remarks["targetformat"] = (
            f"Generation of the the target format file ({target_format}) was carried out by Open Babel version {ctx['versions']['obabel']}"
        )
        remarks["date"] = f'Created on {now.strftime("%Y-%m-%d %H:%M:%S")}'
        remark_string += f"REMARK    Compound: {tautomer['key']}\n"
        remark_string += generate_remarks(remarks, target_format=target_format) + "\n"
    elif target_format == "mol2":
        remarks = tautomer["remarks"].copy()
        remarks["targetformat"] = (
            f"Generation of the the target format file ({target_format}) was carried out by Open Babel version {ctx['versions']['obabel']}"
        )
        remarks["date"] = f'Created on {now.strftime("%Y-%m-%d %H:%M:%S")}'
        remark_string = generate_remarks(remarks, target_format=target_format) + "\n"

    # Output the final file

    if target_format == "smi":
        # get rid of the source file comment that Open Babel adds
        smi_string = lines[0].split()[0]
        with open(output_file, "w") as write_file:
            write_file.write(f"{smi_string}\n")

    else:
        with open(output_file, "w") as write_file:
            if remark_string != "":
                write_file.write(remark_string)

            for line in lines:
                if target_format == "pdb" or target_format == "pdbqt":
                    if re.search(
                        r"TITLE|SOURCE|KEYWDS|EXPDTA|REVDAT|HEADER|AUTHOR", line
                    ):
                        continue
                    line = re.sub(r" UN[LK] ", " LIG ", line)
                    if re.search(r"^\s*$", line):
                        continue
                    if re.search(r"REMARK\s+Name", line):
                        continue

                # OpenBabel often has local path information that we can remove
                line = re.sub(rf"{input_pdb_file}", tautomer["key"], line)

                write_file.write(line)

    logging.debug(f"Finished the target format of {target_format}")
    tautomer["timers"].append(
        [f"obabel_generate_{target_format}", time.perf_counter() - step_timer_start]
    )
