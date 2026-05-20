from __future__ import annotations

from pathlib import Path
from typing import List

from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists,get_stem,check_filename_length,load_protein_structure,write_json_from_dict_inline_leaf_lists

from ..algorithms.clean_algorithms import check_cleaned_structure
from ..algorithms.dock_algorithms import dock_multiple_substrates_from_structure,save_docking_results_and_generate_dock_report
from ..utils.common_utils import get_optimized_filename


def _parse_float_triplet(value: str, parameter_name: str, logger: Logger) -> List[float] | None:
    if not isinstance(value, str) or not value.strip():
        logger.print(f"[ERROR] {parameter_name} is empty.")
        return None

    part_list = [x.strip() for x in value.split(",")]
    if len(part_list) != 3 or any(not x for x in part_list):
        logger.print(f"[ERROR] {parameter_name} must contain exactly 3 values separated by ','.")
        return None

    try:
        return [float(x) for x in part_list]
    except Exception:
        logger.print(f"[ERROR] {parameter_name} must contain numeric values.")
        return None


def run_dock_service(
    input_path: str | Path,
    substrate_names: str,
    substrate_dir: str | Path,
    output_dir: str | Path,
    max_docking_attempt_num: int = 20,
    early_stop: bool = False,
    exhaustiveness: int = 16,
    cpu: int = 0,
    min_rad: float = 1.8,
    max_rad: float = 6.2,
    min_volume: int = 50,
    catalytic_residue: int | None = None,
    catalytic_site_coord: str | None = None,
    box_size: str | None = None,
) -> bool:
    logger = Logger(output_dir)
    logger.print(f"[INFO] Dock processing started: {input_path}")

    if max_docking_attempt_num <= 0 or max_docking_attempt_num > 100 or exhaustiveness <= 0 or exhaustiveness > 64:
        logger.print(
            f"[ERROR] Invalid docking parameters. Require: max_docking_attempt_num (1–100), exhaustiveness (1–64).")
        return False

    catalytic_site_coord_text = str(catalytic_site_coord).strip() if catalytic_site_coord is not None else ""
    box_size_text = str(box_size).strip() if box_size is not None else ""
    use_manual_box = catalytic_residue is not None or bool(catalytic_site_coord_text)

    if catalytic_residue is not None and catalytic_site_coord_text:
        logger.print("[ERROR] --catalytic_residue and --catalytic_site_coord cannot be used together.")
        return False

    if use_manual_box and not box_size_text:
        logger.print("[ERROR] --box_size is required when --catalytic_residue or --catalytic_site_coord is provided.")
        return False

    catalytic_site_coord_list: List[float] | None = None
    box_size_list: List[float] | None = None

    if catalytic_residue is not None and catalytic_residue <= 0:
        logger.print("[ERROR] --catalytic_residue must be a positive integer aa_id.")
        return False

    if catalytic_site_coord_text:
        catalytic_site_coord_list = _parse_float_triplet(catalytic_site_coord_text, "--catalytic_site_coord", logger)
        if catalytic_site_coord_list is None:
            return False

    if use_manual_box:
        box_size_list = _parse_float_triplet(box_size_text, "--box_size", logger)
        if box_size_list is None:
            return False
        if any(x <= 0.0 for x in box_size_list):
            logger.print("[ERROR] --box_size values must be positive.")
            return False
    elif min_rad < 1.2 or min_volume <= 20 or min_rad >= max_rad:
        logger.print(
            f"[ERROR] Invalid pocket detection parameters. Require: min_rad ≥ 1.2, max_rad > min_rad, min_volume > 20.")
        return False

    input_path = Path(input_path)
    substrate_dir = Path(substrate_dir)
    output_dir = Path(output_dir)

    if not file_exists(input_path):
        logger.print(f"[ERROR] Input not found: {input_path}")
        return False

    if not substrate_names or not str(substrate_names).strip():
        logger.print("[ERROR] substrate_names is empty.")
        return False

    if not substrate_dir.exists() or not substrate_dir.is_dir():
        logger.print(f"[ERROR] Invalid substrate_dir: {substrate_dir}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    name = get_stem(input_path)
    if not check_filename_length(name, logger):
        return False
    logger.print(f"[INFO] Protein name resolved: {name}")

    structure = load_protein_structure(input_path, name, logger)
    if structure is None:
        logger.print(f"[ERROR] Failed to load structure: {input_path}")
        return False
    logger.print("[INFO] Structure loaded")

    if not check_cleaned_structure(structure, logger):
        return False
    logger.print("[INFO] Structure checked")

    logger.print("[INFO] Docking workflow started")
    docking_result_list = dock_multiple_substrates_from_structure(
        struct=structure,
        substrate_names=substrate_names,
        substrate_dir=substrate_dir,
        logger=logger,
        max_docking_attempt_num=max_docking_attempt_num,
        early_stop=early_stop,
        exhaustiveness=exhaustiveness,
        cpu=cpu,
        min_rad=min_rad,
        max_rad=max_rad,
        min_volume=min_volume,
        catalytic_residue=catalytic_residue,
        catalytic_site_coord_list=catalytic_site_coord_list,
        manual_box_size_list=box_size_list,
    )
    if docking_result_list is None:
        return False

    logger.print(f"[INFO] Docking finished")

    logger.print("[INFO] Saving docking results and generating report")
    report = save_docking_results_and_generate_dock_report(
        docking_result_list=docking_result_list,
        struct=structure,
        protein_name=name,
        output_dir=output_dir,
        logger=logger,
    )
    if report is None:
        return False

    json_name = f"dock_report_{name}_{substrate_names}.json"
    json_name = get_optimized_filename(json_name)
    json_report_path = output_dir / json_name
    write_json_from_dict_inline_leaf_lists(report, json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Dock processing finished")
    return True