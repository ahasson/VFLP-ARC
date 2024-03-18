#######################
# Step 1: Desalt


def get_smi_string(ligand):
    return ligand["smi"].split()[0]


def desalt(ctx, ligand):

    ligand["number_of_fragments"] = 1
    ligand["remarks"]["desalting"] = ""
    if ctx["config"]["desalting"] == "true":
        # Number of fragments in SMILES

        smi_string = get_smi_string(ligand)
        smi_string_parts = smi_string.split(".")

        if len(smi_string_parts) == 1:
            # Not a salt
            ligand["status_sub"].append(
                ["desalt", {"state": "success", "text": "untouched"}]
            )
            ligand["smi_desalted"] = smi_string
            ligand["remarks"][
                "desalting"
            ] = "The ligand was originally not a salt, therefore no desalting was carried out."

        elif len(smi_string_parts) > 1:
            # need to desalt
            sorted_smi_string_parts = sorted(smi_string_parts, key=len)
            ligand["smallest_fragment"] = sorted_smi_string_parts[0]
            ligand["largest_fragment"] = sorted_smi_string_parts[1]
            ligand["number_of_fragments"] = len(smi_string_parts)
            ligand["status_sub"].append(
                ["desalt", {"state": "success", "text": "genuine"}]
            )
            ligand["remarks"][
                "desalting"
            ] = "The ligand was desalted by extracting the largest organic fragment (out of {len(smi_string_parts)}) from the original structure."
            ligand["smi_desalted"] = ligand["largest_fragment"]

        else:
            # this failed...
            raise RuntimeError(f"Failed to parse '{smi_string_parts}'")

    else:
        smi_string = get_smi_string(ligand)
        ligand["smi_desalted"] = smi_string
