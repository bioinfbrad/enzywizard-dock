from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple
import re

from rdkit import Chem
from Bio.PDB.Structure import Structure
from ..utils.logging_utils import Logger
import copy
from ..utils.common_utils import get_optimized_filename

def get_sdf_atom_info_from_mol(mol: Chem.Mol, logger: Logger) -> List[Dict[str, Any]] | None:
    if mol is None:
        logger.print("[ERROR] Input Mol is None.")
        return None

    try:
        if mol.GetNumConformers() == 0:
            logger.print("[ERROR] Mol has no 3D conformer.")
            return None

        atom_info_list: List[Dict[str, Any]] = []

        for atom in mol.GetAtoms():
            atom_info_list.append(
                {
                    "atom_index": int(atom.GetIdx() + 1),
                    "atom_name": str(atom.GetSymbol()).upper(),
                }
            )

        return atom_info_list

    except Exception:
        logger.print("[ERROR] Failed to extract atom information from Mol.")
        return None


def get_pdbqt_atom_info_from_lines(lines: List[str]) -> List[Dict[str, Any]] | None:
    try:
        atom_info_list: List[Dict[str, Any]] = []

        for line in lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            if len(line) < 54:
                return None

            try:
                pdbqt_atom_index = int(line[6:11].strip())
                pdbqt_atom_name = line[12:16].strip().upper()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except Exception:
                return None

            atom_info_list.append(
                {
                    "pdbqt_atom_index": pdbqt_atom_index,
                    "pdbqt_atom_name": pdbqt_atom_name,
                    "x": x,
                    "y": y,
                    "z": z,
                }
            )

        return atom_info_list

    except Exception:
        return None


def get_pdbqt_index_mapping(pdbqt_path: str | Path,logger: Logger) -> List[Dict[str, Any]] | None:
    try:
        pdbqt_path = Path(pdbqt_path)

        if not pdbqt_path.exists() or pdbqt_path.stat().st_size <= 0:
            logger.print("[ERROR] Invalid input PDBQT file.")
            return None

        lines = pdbqt_path.read_text(encoding="utf-8", errors="replace").splitlines()

        atom_info_list = get_pdbqt_atom_info_from_lines(lines)
        if atom_info_list is None or len(atom_info_list) == 0:
            logger.print("[ERROR] Failed to read atom information from PDBQT file.")
            return None

        pdbqt_index_to_name: Dict[int, str] = {}
        for item in atom_info_list:
            pdbqt_index_to_name[int(item["pdbqt_atom_index"])] = str(item["pdbqt_atom_name"])

        mapping_numbers: List[int] = []
        for line in lines:
            upper_line = line.upper()
            if "REMARK" not in upper_line:
                continue
            if not upper_line.startswith("REMARK"):
                continue
            if "INDEX MAP" not in upper_line:
                continue

            nums = re.findall(r"\d+", line)
            mapping_numbers.extend(int(x) for x in nums)

        if len(mapping_numbers) < 2 or len(mapping_numbers) % 2 != 0:
            logger.print("[ERROR] Invalid INDEX MAP format.")
            return None


        mapping_info_list: List[Dict[str, Any]] = []

        for i in range(0, len(mapping_numbers), 2):
            original_atom_index = int(mapping_numbers[i])
            pdbqt_atom_index = int(mapping_numbers[i + 1])

            if pdbqt_atom_index not in pdbqt_index_to_name:
                logger.print("[ERROR] PDBQT atom index in mapping not found in atom records.")
                return None

            mapping_info_list.append(
                {
                    "original_atom_index": original_atom_index,
                    "original_atom_name": "",
                    "pdbqt_atom_index": pdbqt_atom_index,
                    "pdbqt_atom_name": pdbqt_index_to_name[pdbqt_atom_index],
                }
            )

        mapping_info_list.sort(key=lambda x: int(x["pdbqt_atom_index"]))

        return mapping_info_list

    except Exception:
        logger.print("[ERROR] Failed to parse index mapping from PDBQT file.")
        return None


def get_pose_ligand_block_list(pose_string: str,logger: Logger) -> List[List[str]] | None:
    if not isinstance(pose_string, str) or len(pose_string.strip()) == 0:
        logger.print("[ERROR] Invalid pose string.")
        return None

    try:
        lines = pose_string.splitlines()
        ligand_block_list: List[List[str]] = []
        current_block: List[str] = []
        in_ligand_block = False

        for line in lines:
            stripped = line.strip()
            upper_line = stripped.upper()

            if upper_line == "ROOT":
                if in_ligand_block:
                    logger.print("[ERROR] Nested ROOT found.")
                    return None
                current_block = [line]
                in_ligand_block = True
                continue

            if in_ligand_block:
                current_block.append(line)

                if upper_line.startswith("TORSDOF"):
                    ligand_block_list.append(current_block)
                    current_block = []
                    in_ligand_block = False

        if in_ligand_block:
            logger.print("[ERROR] Incomplete ligand block (missing TORSDOF).")
            return None

        if len(ligand_block_list) == 0:
            logger.print("[ERROR] No ligand blocks found in pose string.")
            return None

        return ligand_block_list

    except Exception:
        logger.print("[ERROR] Failed to split pose string into ligand blocks.")
        return None


def get_pose_for_substrate_atom_info(
    substrate_name: str,
    ligand_order_index: int,
    pose_string: str,
    original_atom_info_list: List[Dict[str, Any]],
    mapping_info_list: List[Dict[str, Any]],
    logger: Logger,
) -> Dict[str, Any] | None:
    if not substrate_name:
        logger.print("[ERROR] substrate_name is empty.")
        return None

    if ligand_order_index < 0:
        logger.print("[ERROR] ligand_order_index must be non-negative.")
        return None

    if not isinstance(original_atom_info_list, list) or len(original_atom_info_list)==0:
        logger.print("[ERROR] Invalid original_atom_info_list.")
        return None

    if not isinstance(mapping_info_list, list) or len(mapping_info_list) == 0:
        logger.print("[ERROR] Invalid mapping_info_list.")
        return None

    original_index_set = set(item["atom_index"] for item in original_atom_info_list)
    mapping_index_set = set(item["original_atom_index"] for item in mapping_info_list)

    if not mapping_index_set.issubset(original_index_set):
        logger.print("[ERROR] Mapping contains invalid original atom indices.")
        return None

    ligand_block_list = get_pose_ligand_block_list(pose_string, logger)
    if ligand_block_list is None:
        return None

    if ligand_order_index >= len(ligand_block_list):
        logger.print(f"[ERROR] ligand_order_index {ligand_order_index} is out of range.")
        return None

    try:
        original_index_to_name: Dict[int, str] = {}
        for item in original_atom_info_list:
            atom_index = int(item.get("atom_index", 0))
            atom_name = str(item.get("atom_name", "")).upper()

            if atom_index <= 0 or not atom_name:
                logger.print("[ERROR] Invalid original atom information.")
                return None

            original_index_to_name[atom_index] = atom_name

        enriched_mapping_info_list: List[Dict[str, Any]] = []
        for item in mapping_info_list:
            original_atom_index = int(item.get("original_atom_index", 0))
            pdbqt_atom_index = int(item.get("pdbqt_atom_index", 0))
            pdbqt_atom_name = str(item.get("pdbqt_atom_name", "")).upper()

            if (
                original_atom_index <= 0
                or pdbqt_atom_index <= 0
                or not pdbqt_atom_name
                or original_atom_index not in original_index_to_name
            ):
                logger.print("[ERROR] Invalid mapping information.")
                return None

            enriched_mapping_info_list.append(
                {
                    "original_atom_index": original_atom_index,
                    "original_atom_name": original_index_to_name[original_atom_index],
                    "pdbqt_atom_index": pdbqt_atom_index,
                    "pdbqt_atom_name": pdbqt_atom_name,
                }
            )

        enriched_mapping_info_list.sort(key=lambda x: int(x["pdbqt_atom_index"]))

        expected_index_set = set(int(item["pdbqt_atom_index"]) for item in enriched_mapping_info_list)

        expected_pdbqt_atom_name_list = [str(item["pdbqt_atom_name"]).upper() for item in enriched_mapping_info_list]

        block_lines = ligand_block_list[ligand_order_index]
        matched_atom_info_list = get_pdbqt_atom_info_from_lines(block_lines)
        if matched_atom_info_list is None:
            logger.print("[ERROR] Failed to parse atom information from pose ligand block.")
            return None

        if len(matched_atom_info_list) != len(mapping_info_list):
            logger.print("[ERROR] Atom count mismatch between pose and mapping.")
            return None

        matched_atom_info_list.sort(key=lambda x: int(x["pdbqt_atom_index"]))

        pose_index_set = set(int(item["pdbqt_atom_index"]) for item in matched_atom_info_list)
        if pose_index_set != expected_index_set:
            logger.print(f"[ERROR] PDBQT atom index mismatch for substrate: {substrate_name}")
            return None

        pose_pdbqt_atom_name_list = [
            str(item["pdbqt_atom_name"]).upper()
            for item in matched_atom_info_list
        ]
        if pose_pdbqt_atom_name_list != expected_pdbqt_atom_name_list:
            logger.print(f"[ERROR] PDBQT atom name mismatch for substrate: {substrate_name}")
            return None

        pdbqt_index_to_pose_atom: Dict[int, Dict[str, Any]] = {}
        for item in matched_atom_info_list:
            pdbqt_index_to_pose_atom[int(item["pdbqt_atom_index"])] = item

        docked_atom_info_list: List[Dict[str, Any]] = []
        for item in enriched_mapping_info_list:
            original_atom_index = int(item["original_atom_index"])
            pdbqt_atom_index = int(item["pdbqt_atom_index"])

            if pdbqt_atom_index not in pdbqt_index_to_pose_atom:
                logger.print("[ERROR] Matched ligand block is inconsistent with mapping.")
                return None

            pose_atom = pdbqt_index_to_pose_atom[pdbqt_atom_index]

            docked_atom_info_list.append(
                {
                    "original_atom_index": original_atom_index,
                    "original_atom_name": str(item["original_atom_name"]).upper(),
                    "pdbqt_atom_index": pdbqt_atom_index,
                    "pdbqt_atom_name": str(item["pdbqt_atom_name"]).upper(),
                    "x": float(pose_atom["x"]),
                    "y": float(pose_atom["y"]),
                    "z": float(pose_atom["z"]),
                }
            )

        docked_atom_info_list.sort(key=lambda x: int(x["original_atom_index"]))

        return {
            "substrate_name": substrate_name,
            "atom_info_list": docked_atom_info_list,
        }

    except Exception:
        logger.print(f"[ERROR] Failed to parse pose for substrate: {substrate_name}")
        return None

def split_vina_pose_string(pose_string: str,logger: Logger) -> List[str] | None:
    if not isinstance(pose_string, str):
        logger.print("[ERROR] pose_string must be a string.")
        return None

    text = pose_string.strip()
    if not text:
        logger.print("[ERROR] Vina returned empty pose string.")
        return []

    lines = text.splitlines()
    pose_string_list: List[str] = []

    current_block: List[str] = []
    in_model = False

    try:
        for line in lines:
            if line.startswith("MODEL"):
                if len(current_block) > 0:
                    logger.print("[ERROR] Found a new MODEL before closing previous ENDMDL.")
                    return None

                in_model = True
                current_block.append(line)
                continue

            if line.startswith("ENDMDL"):
                if not in_model:
                    logger.print("[ERROR] Found ENDMDL before MODEL.")
                    return None

                current_block.append(line)
                pose_string_list.append("\n".join(current_block).strip() + "\n")
                current_block = []
                in_model = False
                continue

            if in_model:
                current_block.append(line)

        if in_model or len(current_block) > 0:
            logger.print("[ERROR] Incomplete MODEL/ENDMDL block in pose_string.")
            return None

        if len(pose_string_list) == 0:
            logger.print("[ERROR] No valid MODEL/ENDMDL pose block found in pose_string.")
            return None

        return pose_string_list

    except Exception:
        logger.print("[ERROR] Failed to split Vina pose string.")
        return None

def get_substrate_sdf_path_group_dict(substrate_names: str,substrate_dir: str | Path,logger: Logger) -> Tuple[List[str], Dict[str, List[str]]] | None:
    if not isinstance(substrate_names, str):
        logger.print("[ERROR] substrate_names must be a string.")
        return None

    if not substrate_names.strip():
        logger.print("[ERROR] substrate_names is empty.")
        return None

    if not isinstance(substrate_dir, (str, Path)):
        logger.print("[ERROR] substrate_dir must be a str or Path.")
        return None

    substrate_dir = Path(substrate_dir)

    if not substrate_dir.exists():
        logger.print(f"[ERROR] substrate_dir does not exist: {substrate_dir}")
        return None

    if not substrate_dir.is_dir():
        logger.print(f"[ERROR] substrate_dir is not a directory: {substrate_dir}")
        return None

    try:
        substrate_name_list = [item.strip() for item in substrate_names.split(",")]
        substrate_name_list = [item for item in substrate_name_list if item]

        if len(substrate_name_list) == 0:
            logger.print("[ERROR] No valid substrate name was found in substrate_names.")
            return None

        if len(set(substrate_name_list)) != len(substrate_name_list):
            logger.print("[ERROR] Duplicate substrate names are not allowed.")
            return None

        sdf_path_list = sorted(substrate_dir.glob("*.sdf"))

        substrate_to_sdf_path_list_dict: Dict[str, List[str]] = {}

        for substrate_name in substrate_name_list:
            matched_pairs: List[Tuple[int, str]] = []

            substrate_file_stem = get_optimized_filename(substrate_name)
            if not substrate_file_stem:
                logger.print(f"[ERROR] Invalid substrate filename generated from substrate: {substrate_name}")
                return None

            pattern = re.compile(rf"^{re.escape(substrate_file_stem)}_(\d+)$")

            for sdf_path in sdf_path_list:
                stem = sdf_path.stem

                if stem == substrate_file_stem:
                    matched_pairs.append((0, str(sdf_path)))
                    continue
                m = pattern.match(stem)
                if m:
                    n = int(m.group(1))
                    matched_pairs.append((n, str(sdf_path)))

            if len(matched_pairs) == 0:
                logger.print(f"[ERROR] No SDF files were found for substrate: {substrate_name}")
                return None

            matched_pairs.sort(key=lambda x: x[0])

            matched_sdf_path_list = [path for _, path in matched_pairs]

            substrate_to_sdf_path_list_dict[substrate_name] = matched_sdf_path_list

        return substrate_name_list, substrate_to_sdf_path_list_dict

    except Exception:
        logger.print("[ERROR] Failed to group substrate SDF files.")
        return None

def compute_ligand_centroid(atom_info_list: List[Dict[str, Any]],logger: Logger) -> List[float] | None:

    if not isinstance(atom_info_list, list) or len(atom_info_list) == 0:
        logger.print("[ERROR] Invalid atom_info_list for centroid calculation.")
        return None

    try:
        x_sum, y_sum, z_sum = 0.0, 0.0, 0.0

        for item in atom_info_list:
            x_sum += float(item["x"])
            y_sum += float(item["y"])
            z_sum += float(item["z"])

        n = len(atom_info_list)

        return [
            x_sum / n,
            y_sum / n,
            z_sum / n,
        ]

    except Exception:
        logger.print("[ERROR] Failed to compute ligand centroid.")
        return None

def remove_hydrogens_from_structure(struct: Structure, logger) -> Structure | None:
    if not isinstance(struct, Structure):
        logger.print("[ERROR] struct must be a Bio.PDB Structure.")
        return None

    try:
        new_struct = copy.deepcopy(struct)

        for model in new_struct:
            for chain in model:
                for residue in chain:
                    atom_id_list_to_remove = []

                    for atom in residue:
                        element = getattr(atom, "element", None)
                        atom_name = str(atom.get_name()).strip()

                        if element == "H" or (element is None and atom_name.startswith("H")):
                            atom_id_list_to_remove.append(atom.id)

                    for atom_id in atom_id_list_to_remove:
                        residue.detach_child(atom_id)

        return new_struct

    except Exception:
        logger.print("[ERROR] Failed to remove hydrogens from structure")
        return None