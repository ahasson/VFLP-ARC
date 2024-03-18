from pathlib import Path


def write_file_single(filename, data_to_write):
    with open(filename, "wt", encoding="utf8") as write_file:
        write_file.write(data_to_write)
        write_file.write("\n")


def get_intermediate_dir_ligand(ctx) -> Path:

    base_temp_dir = ctx["temp_dir"].name

    if ctx["store_all_intermediate_logs"] == "true":
        base_temp_dir = ctx["collection_temp_dir"].name

    # Intermediate log storage
    output_file_parts = [
        base_temp_dir,
        "intermediate",
        ctx["metatranche"],
        ctx["tranche"],
        ctx["collection_name"],
        ctx["ligand_key"],
    ]

    output_dir: Path = Path("/".join(output_file_parts))
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir
