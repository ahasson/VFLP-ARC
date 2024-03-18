import re

def halogencount(smi):
    return smi.count("F") + smi.count("Cl") + smi.count("Br") + smi.count("I")


# Q and J are the only two letters not in the periodic table
# TODO: Ask @Christoph about the Xx -> should this be Qq to avoid Xenon etc?
REPLACEMENT_LETTER: str = "Q"
def electronegativeatomcount(smi):
    smi = smi.replace("Na", "")
    smi = smi.replace("Cl", "X")
    smi = smi.replace("Si", "")

    return len(re.findall("[NnOoSsPpFfXxBbIi]", smi))


def NOcount(smi):
    smi = smi.replace("Na", "")
    return smi.count("N") + smi.count("O") + smi.count("n") + smi.count("o")


def sulfurcount(smi):
    smi = smi.replace("Si", "")
    return smi.count("S")


def positivecharge(smi):
    smi = smi.replace("+2", "++")
    return smi.count("+")


def negativecharge(smi):
    smi = smi.replace("-2", "--")
    return smi.count("-")


def formalcharge(smi):
    return positivecharge(smi) - negativecharge(smi)


def get_file_data(ligand, attribute):
    if attribute in ligand["file_data"]:
        return ligand["file_data"][attribute]
    else:
        raise RuntimeError(
            f"Asked for attribute '{attribute}' that does not exist in file_data"
        )