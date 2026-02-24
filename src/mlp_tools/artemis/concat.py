"""Concatenate Artemis metadata next to POSCAR."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path


def _find_environment(cwd: Path) -> dict[str, str | bool | None]:
    dir_name = cwd.name
    parent_name = cwd.parent.name

    is_swap = False
    shift_vals_rel: str | None = None
    struc_dat_rel: str | None = None

    if parent_name == "DSHIFT":
        shift_vals_rel = "../shift_vals.txt"
        struc_dat_rel = "../../struc_dat.txt"
    elif parent_name == "DSWAP":
        is_swap = True
        shift_vals_rel = "../../../shift_vals.txt"
        struc_dat_rel = "../../../../struc_dat.txt"
    elif dir_name == "DSHIFT":
        shift_vals_rel = "shift_vals.txt"
        struc_dat_rel = "../struc_dat.txt"
    elif dir_name == "DSWAP":
        is_swap = True
        shift_vals_rel = "shift_vals.txt"
        struc_dat_rel = "../struc_dat.txt"

    return {
        "is_swap": is_swap,
        "shift_vals_relpath": shift_vals_rel,
        "struc_dat_relpath": struc_dat_rel,
    }


def _parse_shift_vals(filepath: Path, target_interface_num: str) -> dict | None:
    if not filepath.is_file():
        return None

    pattern = re.compile(r"(\S+)\s*\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^,]+)\s*\)")
    with filepath.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            match = pattern.match(line)
            if match and match.group(1).zfill(2) == target_interface_num.zfill(2):
                return {
                    "interface_num": match.group(1),
                    "shift_a": float(match.group(2).strip()),
                    "shift_b": float(match.group(3).strip()),
                    "shift_c": float(match.group(4).strip()),
                }
    return None


def _parse_struc_dat(filepath: Path) -> dict:
    data: dict = {}
    if not filepath.is_file():
        return data

    with filepath.open(encoding="utf-8") as handle:
        for line in handle:
            line_str = line.strip()
            if "vector mismatch (%)" in line_str:
                data["vector_mismatch_percent"] = float(line_str.split("=", 1)[1].strip())
            elif "angle mismatch" in line_str:
                data["angle_mismatch_deg"] = float(line_str.split("=", 1)[1].strip())
            elif "area mismatch (%)" in line_str:
                data["area_mismatch_percent"] = float(line_str.split("=", 1)[1].strip())
            elif "Lower crystal Miller plane:" in line_str:
                vals = line_str.split(":", 1)[1].split()
                data["lower_miller_plane"] = [int(v) for v in vals]
            elif "Upper crystal Miller plane:" in line_str:
                vals = line_str.split(":", 1)[1].split()
                data["upper_miller_plane"] = [int(v) for v in vals]

    return data


def main() -> None:
    cwd = Path(os.getcwd())
    env = _find_environment(cwd)

    shift_vals_rel = env["shift_vals_relpath"]
    struc_dat_rel = env["struc_dat_relpath"]
    is_swap = bool(env["is_swap"])

    current_dir = cwd.name
    if is_swap:
        shift_folder = cwd.parent.parent.name
        shift_interface_num = shift_folder.lstrip("D").zfill(2)
        swap_interface_num = current_dir.lstrip("D").zfill(2)
    else:
        shift_interface_num = current_dir.lstrip("D").zfill(2)
        swap_interface_num = None

    shift_data = None
    if shift_vals_rel is not None:
        shift_data = _parse_shift_vals(cwd / shift_vals_rel, shift_interface_num)
    struct_data = {}
    if struc_dat_rel is not None:
        struct_data = _parse_struc_dat(cwd / struc_dat_rel)

    result = {
        "is_swap": is_swap,
        "shift_vals_path": shift_vals_rel,
        "struc_dat_path": struc_dat_rel,
        "shift_data": shift_data,
        "struct_data": struct_data,
    }
    if swap_interface_num is not None:
        result["swap_number"] = swap_interface_num

    with (cwd / "poscar_data.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)


if __name__ == "__main__":
    main()
