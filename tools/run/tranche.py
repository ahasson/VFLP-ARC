#######################
# Step 5: Assign Tranche
import logging
import re

from run.obabel import (
    run_obabel_hba,
    run_obabel_hbd,
    run_obabel_attributes,
    write_file_single,
)
from run.chemaxon import run_cxcalc_attributes
from run.rdkit import run_rdkit_attributes
from run.props.base import (
    formalcharge,
    positivecharge,
    negativecharge,
    halogencount,
    sulfurcount,
    NOcount,
    electronegativeatomcount,
    get_file_data,
)


def get_mol_attributes():
    return {
        "mw_obabel": {"prog": "obabel", "prog_name": "mol_weight", "val": "INVALID"},
        "logp_obabel": {"prog": "obabel", "prog_name": "logP", "val": "INVALID"},
        "tpsa_obabel": {"prog": "obabel", "prog_name": "PSA", "val": "INVALID"},
        "atomcount_obabel": {
            "prog": "obabel",
            "prog_name": "num_atoms",
            "val": "INVALID",
        },
        "bondcount_obabel": {
            "prog": "obabel",
            "prog_name": "num_bonds",
            "val": "INVALID",
        },
        "mr_obabel": {"prog": "obabel", "prog_name": "MR", "val": "INVALID"},
        "mw_jchem": {"prog": "cxcalc", "prog_name": "mass", "val": "INVALID"},
        "logp_jchem": {"prog": "cxcalc", "prog_name": "logp", "val": "INVALID"},
        "hbd_jchem": {"prog": "cxcalc", "prog_name": "donorcount", "val": "INVALID"},
        "hba_jchem": {"prog": "cxcalc", "prog_name": "acceptorcount", "val": "INVALID"},
        "rotb_jchem": {
            "prog": "cxcalc",
            "prog_name": "rotatablebondcount",
            "val": "INVALID",
        },
        "tpsa_jchem": {
            "prog": "cxcalc",
            "prog_name": "polarsurfacearea",
            "val": "INVALID",
        },
        "atomcount_jchem": {
            "prog": "cxcalc",
            "prog_name": "atomcount",
            "val": "INVALID",
        },
        "bondcount_jchem": {
            "prog": "cxcalc",
            "prog_name": "bondcount",
            "val": "INVALID",
        },
        "ringcount": {"prog": "cxcalc", "prog_name": "ringcount", "val": "INVALID"},
        "aromaticringcount": {
            "prog": "cxcalc",
            "prog_name": "aromaticringcount",
            "val": "INVALID",
        },
        "mr_jchem": {"prog": "cxcalc", "prog_name": "refractivity", "val": "INVALID"},
        "fsp3": {"prog": "cxcalc", "prog_name": "fsp3", "val": "INVALID"},
        "chiralcentercount": {
            "prog": "cxcalc",
            "prog_name": "chiralcentercount",
            "val": "INVALID",
        },
        "logd": {"prog": "cxcalc", "prog_name": "logd", "val": "INVALID"},
        "logs": {"prog": "cxcalc", "prog_name": "logs", "val": "INVALID"},
        "doublebondstereoisomercount_jchem": {
            "prog": "cxcalc",
            "prog_name": "doublebondstereoisomercount",
            "val": "INVALID",
        },
        "aromaticproportion_jchem": {
            "prog": "cxcalc",
            "prog_name": "aromaticproportion",
            "val": "INVALID",
        },
        "qed_rdkit": {"prog": "rdkit", "prog_name": "qed", "val": "INVALID"},
        "scaffold_rdkit": {"prog": "rdkit", "prog_name": "scaffold", "val": "INVALID"},
    }


def generate_attributes(ctx, ligand, tautomer):

    attribute_dict = {}
    attributes_to_generate = {}

    if ctx["config"]["tranche_assignments"] == "true":
        for tranche_type in ctx["config"]["tranche_types"]:
            if tranche_type not in attributes_to_generate:
                attributes_to_generate[tranche_type] = 1
    for attribute_type in ctx["config"]["attributes_to_generate"]:
        if attribute_type not in attributes_to_generate:
            attributes_to_generate[attribute_type] = 1

    smi_file = f"{ctx['temp_dir'].name}/smi_tautomers_{tautomer['key']}.smi"
    write_file_single(smi_file, tautomer["smi_protomer"])

    attributes = get_mol_attributes()

    obabel_run = 0
    cxcalc_run = 0
    rdkit_run = 0

    for tranche_type in attributes_to_generate:
        if tranche_type in attributes:
            if attributes[tranche_type]["prog"] == "obabel":
                obabel_run = 1
            elif attributes[tranche_type]["prog"] == "cxcalc":
                cxcalc_run = 1
                nailgun_port = ctx["config"]["nailgun_port"]
                nailgun_host = ctx["config"]["nailgun_host"]
            elif attributes[tranche_type]["prog"] == "rdkit":
                rdkit_run = 1
        else:
            # We need to generate this one at a time
            attribute_dict[tranche_type] = generate_single_attribute(
                ctx, tranche_type, ligand, tautomer, smi_file
            )

    if obabel_run == 1:
        run_obabel_attributes(ctx, tautomer, smi_file, attributes)

    if cxcalc_run == 1:
        if "use_cxcalc_helper" in ctx["config"]:
            use_single = int(ctx["config"]["use_cxcalc_helper"])
        else:
            use_single = 0
        run_cxcalc_attributes(
            tautomer,
            smi_file,
            nailgun_port,
            nailgun_host,
            attributes_to_generate.keys(),
            attributes,
            use_single=use_single,
        )

    if rdkit_run == 1:
        run_rdkit_attributes(
            ctx,
            tautomer,
            tautomer["smi_protomer"],
            attributes_to_generate.keys(),
            attributes,
        )

    # Put all of the attributes from obabel, cxcalc, and rdkit into
    # the return dict

    for tranche_type in attributes_to_generate:
        if tranche_type in attributes:
            attribute_dict[tranche_type] = attributes[tranche_type]["val"]

    return attribute_dict


def generate_single_attribute(_ctx, attribute, ligand, tautomer, smi_file):

    if attribute == "enamine_type":
        return get_file_data(ligand, "enamine")
    elif attribute == "hba_obabel":
        return run_obabel_hba(smi_file, tautomer)
    elif attribute == "hbd_obabel":
        return run_obabel_hbd(smi_file, tautomer)
    elif attribute == "formalcharge":
        return formalcharge(tautomer["smi_protomer"])
    elif attribute == "positivechargecount":
        return positivecharge(tautomer["smi_protomer"])
    elif attribute == "negativechargecount":
        return negativecharge(tautomer["smi_protomer"])
    elif attribute == "halogencount":
        return halogencount(tautomer["smi_protomer"])
    elif attribute == "sulfurcount":
        return sulfurcount(tautomer["smi_protomer"])
    elif attribute == "NOcount":
        return NOcount(tautomer["smi_protomer"])
    elif attribute == "electronegativeatomcount":
        return electronegativeatomcount(tautomer["smi_protomer"])
    elif attribute == "mw_file":
        return get_file_data(ligand, "mw")
    elif attribute == "logp_file":
        return get_file_data(ligand, "logp")
    elif attribute == "hba_file":
        return get_file_data(ligand, "hba")
    elif attribute == "hbd_file":
        return get_file_data(ligand, "hbd")
    elif attribute == "rotb_file":
        return get_file_data(ligand, "rotb")
    elif attribute == "tpsa_file":
        return get_file_data(ligand, "tpsa")
    elif attribute == "logd_file":
        return get_file_data(ligand, "logd")
    elif attribute == "logs_file":
        return get_file_data(ligand, "logs")
    elif attribute == "heavyatomcount_file":
        return get_file_data(ligand, "heavyatomcount")
    elif attribute == "ringcount_file":
        return get_file_data(ligand, "ringcount")
    elif attribute == "aromaticringcount_file":
        return get_file_data(ligand, "aromaticringcount")
    elif attribute == "mr_file":
        return get_file_data(ligand, "mr")
    elif attribute == "formalcharge_file":
        return get_file_data(ligand, "formalcharge")
    elif attribute == "positivechargecount_file":
        return get_file_data(ligand, "positivecharge")
    elif attribute == "negativechargecount_file":
        return get_file_data(ligand, "negativechargeount")
    elif attribute == "fsp3_file":
        return get_file_data(ligand, "fsp3")
    elif attribute == "chiralcentercount_file":
        return get_file_data(ligand, "chiralcentercount")
    elif attribute == "halogencount_file":
        return get_file_data(ligand, "halogencount")
    elif attribute == "sulfurcount_file":
        return get_file_data(ligand, "sulfurcount")
    elif attribute == "NOcount_file":
        return get_file_data(ligand, "NOcount")
    elif attribute == "electronegativeatomcount_file":
        return get_file_data(ligand, "electronegativeatomcount")
    else:
        logging.error(
            "The value '%s'' of the variable used as an attribute is not supported.",
            attribute,
        )
        raise RuntimeError(
            f"The value '{attribute}'' of the variable used as an attribute is not supported."
        )


def tranche_assignment(ctx, _ligand, tautomer):

    tranche_value = ""
    tranche_string = ""

    tautomer["remarks"]["trancheassignment_attr"] = []

    string_attributes = ["enamine_type"]

    for tranche_type in ctx["config"]["tranche_types"]:

        if tranche_type in tautomer["attr"]:
            tranche_value = tautomer["attr"][tranche_type]
        else:
            logging.error(
                "The value '%s' of the variable tranche_types is not supported.",
                tranche_type,
            )
            raise RuntimeError(
                f"The value '{tranche_type}'' of the variable tranche_types is not supported."
            )

        tranche_value = str(tranche_value)

        if tranche_type in string_attributes:
            # process as string
            letter = assign_character_mapping(
                ctx["config"]["tranche_mappings"][tranche_type], tranche_value
            )
            logging.debug(
                "Assigning %s based on '%s' for type %s. String now: %s",
                letter,
                tranche_value,
                tranche_type,
                tranche_string,
            )
            tautomer["remarks"]["trancheassignment_attr"].append(
                f"{tranche_type}: {tranche_value}"
            )
        else:
            # Make sure we have a valid numerical value
            match = re.search(r"^([0-9+\-eE\.]+)$", tranche_value)
            if match:
                letter = assign_character(
                    ctx["config"]["tranche_partitions"][tranche_type], tranche_value
                )

                logging.debug(
                    "Assigning %s based on '%s' for type %s. String now: %s",
                    letter,
                    tranche_value,
                    tranche_type,
                    tranche_string,
                )
                tautomer["remarks"]["trancheassignment_attr"].append(
                    f"{tranche_type}: {tranche_value}"
                )
            else:
                logging.error(
                    "Invalid result from tranche_type:%s, value was: '%s'",
                    tranche_type,
                    tranche_value,
                )
                raise RuntimeError(
                    f"Invalid result from tranche_type:{tranche_type}, value was: '{tranche_value}'"
                )

        tranche_string += letter

    tautomer["tranche_string"] = tranche_string
    tautomer["remarks"]["tranche_str"] = f"Tranche: {tranche_string}"


def assign_character_mapping(string_mapping, tranche_value):

    if tranche_value in string_mapping:
        return string_mapping[tranche_value]

    return "X"


def assign_character(partitions, tranche_value):

    tranche_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    letter = ""

    logging.debug(
        "Assigning letter for partitions %s, tranche_value: %s",
        partitions,
        tranche_value,
    )

    # If less than the first one...
    if float(tranche_value) <= float(partitions[0]):
        return tranche_letters[0]
    # If larger than the largest partition specified
    if float(tranche_value) > float(partitions[-1]):
        return tranche_letters[len(partitions)]

    # Otherwise go through each one until we are no longer smaller
    for index, partition_value in enumerate(partitions):
        if float(tranche_value) > float(partition_value):
            letter = tranche_letters[index + 1]
        else:
            break

    return letter
