"""VASP setup utility commands for MP1 workflows."""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

import numpy as np


DEFAULT_INCAR_TEMPLATE = """SYSTEM = MP1 interface calculation
ISTART = 0
ICHARG = 2
LWAVE = .FALSE.
LREAL = AUTO
ENCUT = 700
ENAUG = 900
ADDGRID = .TRUE.
LMAXMIX = 6
LASPH = .TRUE.
LMIXTAU = .TRUE.

NCORE = 16

NELM = 100
NBANDS = {NBANDS}
EDIFF = 1d-5
AMIX = 0.6
ISMEAR = 0
SIGMA = 0.01
NELECT = {NELECT}

ISPIN = 1
NUPDOWN = 1

LVTOT = .FALSE.
LVHAR = .TRUE.

ISIF = 6
IBRION = 2
NSW = 100
EDIFFG = -0.02
POTIM = 0.5
SYMPREC = 1e-6
ISYM = 0
"""


def _read_species(poscar_path: Path) -> list[str]:
    lines = poscar_path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 6:
        raise ValueError("POSCAR file too short, missing species line.")
    return lines[5].split()


def _read_atom_numbers(poscar_path: Path) -> list[int]:
    lines = poscar_path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 7:
        raise ValueError("POSCAR file too short, missing atom numbers line.")
    return [int(x) for x in lines[6].split()]


def _resolve_potcar_file(species_dir: Path) -> Path:
    for filename in ("POTCAR", "PSCTR"):
        candidate = species_dir / filename
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"No POTCAR/PSCTR file found in {species_dir}")


def _read_zvals(potcar_path: Path) -> list[float]:
    zvals: list[float] = []
    pattern = re.compile(r"ZVAL\s*=\s*([\d\.]+)")
    with potcar_path.open(encoding="utf-8") as handle:
        for line in handle:
            if "ZVAL" not in line:
                continue
            match = pattern.search(line)
            if match:
                zvals.append(float(match.group(1)))
    return zvals


def _read_lattice_vectors(poscar_path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = poscar_path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 5:
        raise ValueError("POSCAR file too short, missing lattice vectors.")
    scale = float(lines[1].strip())
    lattice = np.array(
        [
            [float(x) for x in lines[2].split()],
            [float(x) for x in lines[3].split()],
            [float(x) for x in lines[4].split()],
        ]
    )
    lattice *= scale
    reciprocal = 2.0 * np.pi * np.linalg.inv(lattice).T
    return lattice, reciprocal


def _generate_kpoints_grid(reciprocal_lattice: np.ndarray, kspacing: float) -> list[int]:
    lengths = np.linalg.norm(reciprocal_lattice, axis=1)
    grid = np.maximum(np.round(lengths / kspacing).astype(int), 1)
    return grid.tolist()


def _load_required_int(path: Path, label: str) -> int:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label} file: {path}")
    try:
        return int(path.read_text(encoding="utf-8").strip())
    except ValueError as exc:
        raise ValueError(f"Invalid integer value in {path}") from exc


def _write_text(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def _cmd_generate_potcar(args: argparse.Namespace) -> int:
    workdir = args.workdir.resolve()
    poscar_path = workdir / "POSCAR"
    output_path = workdir / "POTCAR"

    if not poscar_path.exists():
        raise FileNotFoundError(f"POSCAR not found: {poscar_path}")

    potcar_dir = args.potcar_dir
    if potcar_dir is None:
        env_value = os.getenv("VASP_POTCAR_DIR")
        if env_value:
            potcar_dir = Path(env_value)
    if potcar_dir is None:
        raise ValueError("POTCAR directory required: set --potcar-dir or VASP_POTCAR_DIR")

    potcar_root = potcar_dir.resolve()
    if not potcar_root.exists():
        raise FileNotFoundError(f"POTCAR directory not found: {potcar_root}")

    species_list = _read_species(poscar_path)
    potcar_contents: list[str] = []

    for species in species_list:
        species_dir = potcar_root / species
        if not species_dir.exists():
            raise FileNotFoundError(f"Species directory not found for '{species}': {species_dir}")
        source_file = _resolve_potcar_file(species_dir)
        potcar_contents.append(source_file.read_text(encoding="utf-8"))

    _write_text(output_path, "".join(potcar_contents))
    print(f"Generated POTCAR at {output_path}")
    return 0


def _cmd_calculate_nbands(args: argparse.Namespace) -> int:
    workdir = args.workdir.resolve()
    poscar_path = workdir / "POSCAR"
    potcar_path = workdir / "POTCAR"

    if not poscar_path.exists() or not potcar_path.exists():
        raise FileNotFoundError(f"POSCAR/POTCAR not found in {workdir}")

    atom_counts = _read_atom_numbers(poscar_path)
    zvals = _read_zvals(potcar_path)

    if len(atom_counts) != len(zvals):
        raise ValueError(f"Mismatch: {len(atom_counts)} atom types vs {len(zvals)} ZVAL entries")

    total_valence = int(round(sum(n * z for n, z in zip(atom_counts, zvals))))
    min_bands = total_valence // 2 + (1 if total_valence % 2 else 0)
    nbands = int(math.ceil(min_bands * 1.2))

    _write_text(workdir / "nbands.txt", f"{nbands}\n")
    _write_text(workdir / "nelect.txt", f"{total_valence}\n")

    print(f"Total valence electrons: {total_valence}")
    print(f"Minimum bands: {min_bands}")
    print(f"Recommended NBANDS (+20%): {nbands}")
    return 0


def _cmd_generate_incar(args: argparse.Namespace) -> int:
    workdir = args.workdir.resolve()
    poscar_path = workdir / "POSCAR"
    if not poscar_path.exists():
        raise FileNotFoundError(f"POSCAR not found: {poscar_path}")

    nbands = _load_required_int(workdir / "nbands.txt", "NBANDS")
    nelect = _load_required_int(workdir / "nelect.txt", "NELECT")

    if args.incar_template is not None:
        template_path = args.incar_template.resolve()
        if not template_path.exists():
            raise FileNotFoundError(f"INCAR template not found: {template_path}")
        template = template_path.read_text(encoding="utf-8")
    else:
        template = DEFAULT_INCAR_TEMPLATE

    incar_content = template.format(NBANDS=nbands, NELECT=nelect)
    output_path = workdir / "INCAR"
    _write_text(output_path, incar_content)
    print(f"Generated INCAR at {output_path}")
    return 0


def _cmd_generate_kpoints(args: argparse.Namespace) -> int:
    workdir = args.workdir.resolve()
    poscar_path = workdir / "POSCAR"
    if not poscar_path.exists():
        raise FileNotFoundError(f"POSCAR not found: {poscar_path}")

    direct, reciprocal = _read_lattice_vectors(poscar_path)
    grid = _generate_kpoints_grid(reciprocal, args.kspacing)
    output_path = workdir / "KPOINTS"

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("Automatic KPOINTS generated\n")
        handle.write("0\n")
        if args.kcenter == "gamma":
            handle.write("Gamma\n")
        else:
            handle.write("Monkhorst-Pack\n")
        handle.write(f"{grid[0]} {grid[1]} {grid[2]}\n")
        handle.write("0 0 0\n")

    print(f"Direct lattice lengths: {np.linalg.norm(direct, axis=1)}")
    print(f"Reciprocal lattice lengths: {np.linalg.norm(reciprocal, axis=1)}")
    print(f"KPOINTS grid: {grid}")
    print(f"Generated KPOINTS at {output_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="MP1 VASP Utilities Manager")
    parser.add_argument("--workdir", type=Path, default=Path.cwd(), help="Calculation working directory")
    parser.add_argument("--potcar-dir", type=Path, default=None, help="POTCAR root directory")
    parser.add_argument("--incar-template", type=Path, default=None, help="Path to INCAR template file")

    parser.add_argument("--kcenter", choices=["gamma", "monkhorst"], default="gamma")
    parser.add_argument("--kspacing", type=float, default=0.2)

    parser.add_argument("--generate-potcar", action="store_true", help="Generate POTCAR from POSCAR")
    parser.add_argument("--calculate-nbands", action="store_true", help="Calculate NBANDS/NELECT")
    parser.add_argument("--generate-incar", action="store_true", help="Generate INCAR from template")
    parser.add_argument("--generate-kpoints", action="store_true", help="Generate KPOINTS")
    parser.add_argument("--all", action="store_true", help="Run all setup steps in legacy order")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        any_selected = any(
            [
                args.generate_potcar,
                args.calculate_nbands,
                args.generate_incar,
                args.generate_kpoints,
                args.all,
            ]
        )
        if not any_selected:
            parser.print_help()
            return 0

        if args.generate_potcar or args.all:
            _cmd_generate_potcar(args)
        if args.calculate_nbands or args.all:
            _cmd_calculate_nbands(args)
        if args.generate_incar or args.all:
            _cmd_generate_incar(args)
        if args.generate_kpoints or args.all:
            _cmd_generate_kpoints(args)
        return 0
    except Exception as exc:
        print(f"Error: {exc}")
        return 1
