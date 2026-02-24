"""Artemis workflow orchestration with dry-run safety."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


RUN_MACE_TEMPLATE = """#!/usr/bin/env python3
import os
import sys
from ase.io import read
from ase.optimize import FIRE
from mace.calculators import mace_mp


def run_mace_calculation(poscar_path):
    energy = None
    try:
        macemp = mace_mp(model="large", dispersion=True, default_dtype="float64")
        structure = read(poscar_path, format="vasp")
        structure.set_pbc(True)
        structure.calc = macemp
        opt = FIRE(structure)
        opt.run(fmax=0.01)
        structure.write("CONTCAR_MLP_MACE", format="vasp", direct=True)
        energy = structure.get_potential_energy()
    except Exception as exc:
        print(f"Error reading {poscar_path}: {exc}")
        sys.exit(1)
    return energy


if __name__ == "__main__":
    poscar_file = "POSCAR"
    if not os.path.exists(poscar_file):
        print("POSCAR file not found in the current directory.")
        sys.exit(1)

    energy = run_mace_calculation(poscar_file)
    with open("OUTCAR_MLP_MACE", "w", encoding="utf-8") as handle:
        handle.write(f"MACE Energy: {energy:.6f} eV\\n")

    print(f"MACE calculation complete: {energy:.6f} eV")
"""


SBATCH_MLP_TEMPLATE = """#!/bin/bash
#SBATCH -p pq
#SBATCH -A Research_Project-T127870
#SBATCH --job-name={name}MLPMACE
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -D .
#SBATCH --export=ALL
#SBATCH --mail-user=ah1365@exeter.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=MLPMACE.out
#SBATCH --error=MLPMACE.err

NNODES=$SLURM_JOB_NUM_NODES
PPN=16
np=$((NNODES*PPN))
export OMP_NUM_THREADS=1

module purge
module load intel/2023a

export MPI=$((NNODES*PPN))
export WORK_DIR=`pwd`
export JOB_TIME_START=`date -u`

cd ${{WORK_DIR}}

python run-mace.py

export JOB_TIME_END=`date -u`
echo "Job ended on: ${{JOB_TIME_END}}"
"""


SBATCH_DFT_TEMPLATE = """#!/bin/bash
#SBATCH -p pq
#SBATCH -A Research_Project-T127870
#SBATCH --job-name={name}DFTVasp
#SBATCH --time=24:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH -D .
#SBATCH --export=ALL
#SBATCH --mail-user=ah1365@exeter.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=DFTVasp.out
#SBATCH --error=DFTVasp.err

NNODES=$SLURM_JOB_NUM_NODES
PPN=16
np=$(($SLURM_JOB_NUM_NODES*16))
export OMP_NUM_THREADS=1

module purge
module load intel

export MPI=$(($NNODES*$PPN))
export WORK_DIR=`pwd`
export JOB_HOSTNAME=`hostname`
export JOB_TIME_START=`date -u`

cd ${{WORK_DIR}}
mpirun -genv IMPI_PIN_DOMAIN omp -envall -ppn 16 -np $np vasp_std
"""


CONCAT_SHIM = """#!/usr/bin/env python3
from mlp_tools.artemis.concat import main


if __name__ == "__main__":
    main()
"""


@dataclass(frozen=True)
class ArtemisPaths:
    root: Path
    interfaces_dir: Path
    run_root_mlp: Path
    run_root_dft: Path


def _resolve(path_like: str | None) -> Path:
    if path_like is None:
        return Path.cwd().resolve()
    return Path(path_like).expanduser().resolve()


def _resolve_paths(args: argparse.Namespace) -> ArtemisPaths:
    root = _resolve(args.root)
    interfaces_dir = _resolve(args.interfaces_dir) if args.interfaces_dir else (root / "DINTERFACES")
    run_root_mlp = _resolve(args.run_root_mlp) if args.run_root_mlp else (root / "DRUN_MLP_MACE")
    run_root_dft = _resolve(args.run_root_dft) if args.run_root_dft else (root / "DRUN_DFT_VASP")
    return ArtemisPaths(root=root, interfaces_dir=interfaces_dir, run_root_mlp=run_root_mlp, run_root_dft=run_root_dft)


def _print_paths(paths: ArtemisPaths) -> None:
    print(f"Resolved root: {paths.root}")
    print(f"Resolved interfaces dir: {paths.interfaces_dir}")
    print(f"Resolved MLP run root: {paths.run_root_mlp}")
    print(f"Resolved DFT run root: {paths.run_root_dft}")


def _require_safe_interfaces(path: Path) -> bool:
    resolved = path.resolve()
    if resolved == Path("/"):
        print("Refusing to operate on filesystem root '/'.")
        return False
    if not resolved.exists():
        print(f"Interfaces directory not found: {resolved}")
        return False
    return True


def _iter_poscar_dirs(interfaces_dir: Path):
    for root, _dirs, files in os.walk(interfaces_dir):
        if "POSCAR" in files:
            yield Path(root)


def _extract_miller(struc_dat_path: Path) -> tuple[str | None, str | None]:
    lower = None
    upper = None
    with struc_dat_path.open(encoding="utf-8") as handle:
        for line in handle:
            if "Lower crystal Miller plane:" in line:
                lower = "".join(re.findall(r"-?\d+", line.strip()))
            if "Upper crystal Miller plane:" in line:
                upper = "".join(re.findall(r"-?\d+", line.strip()))
            if lower and upper:
                break
    return lower, upper


def _nearest_struc_dat(path: Path) -> Path | None:
    for ancestor in path.parents:
        candidate = ancestor / "struc_dat.txt"
        if candidate.exists():
            return candidate
    return None


def _write_file(path: Path, content: str, dry_run: bool) -> None:
    if dry_run:
        print(f"[DRY-RUN] write {path}")
        return
    path.write_text(content, encoding="utf-8")


def _delete_file(path: Path, dry_run: bool) -> None:
    if not path.exists():
        return
    if dry_run:
        print(f"[DRY-RUN] delete {path}")
        return
    path.unlink()


def _run_command(cmd: list[str], cwd: Path, dry_run: bool) -> int:
    if dry_run:
        print(f"[DRY-RUN] cwd={cwd} cmd={' '.join(cmd)}")
        return 0
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout.strip())
        if result.stderr:
            print(result.stderr.strip())
    return result.returncode


def _setup_mlp(paths: ArtemisPaths, dry_run: bool) -> int:
    count = 0
    for poscar_dir in _iter_poscar_dirs(paths.interfaces_dir):
        struc_dat = _nearest_struc_dat(poscar_dir)
        if struc_dat is None:
            print(f"Skip (missing struc_dat.txt): {poscar_dir}")
            continue
        lower, upper = _extract_miller(struc_dat)
        if not lower or not upper:
            print(f"Skip (missing Miller planes): {struc_dat}")
            continue
        name = f"{lower}-{upper}"
        _write_file(poscar_dir / "SBATCH_MLP_MACE", SBATCH_MLP_TEMPLATE.format(name=name), dry_run)
        _write_file(poscar_dir / "run-mace.py", RUN_MACE_TEMPLATE, dry_run)
        if not dry_run:
            (poscar_dir / "run-mace.py").chmod(0o755)
        count += 1
    print(f"MLP setup completed in {count} POSCAR directories")
    return 0


def _setup_dft(paths: ArtemisPaths, dry_run: bool) -> int:
    count = 0
    for poscar_dir in _iter_poscar_dirs(paths.interfaces_dir):
        struc_dat = _nearest_struc_dat(poscar_dir)
        if struc_dat is None:
            print(f"Skip (missing struc_dat.txt): {poscar_dir}")
            continue
        lower, upper = _extract_miller(struc_dat)
        if not lower or not upper:
            print(f"Skip (missing Miller planes): {struc_dat}")
            continue
        name = f"{lower}-{upper}"
        _write_file(poscar_dir / "SBATCH_DFT_VASP", SBATCH_DFT_TEMPLATE.format(name=name), dry_run)
        rc = _run_command(["python3", "-m", "mlp_tools.vasp", "--workdir", str(poscar_dir), "--all"], poscar_dir, dry_run)
        if rc != 0:
            print(f"Skip (mp1-vasp failed): {poscar_dir}")
            continue
        count += 1
    print(f"DFT setup completed in {count} POSCAR directories")
    return 0


def _setup_structure(paths: ArtemisPaths, dry_run: bool) -> int:
    count = 0
    for poscar_dir in _iter_poscar_dirs(paths.interfaces_dir):
        _write_file(poscar_dir / "concat_artemis_data.py", CONCAT_SHIM, dry_run)
        if not dry_run:
            (poscar_dir / "concat_artemis_data.py").chmod(0o755)
        count += 1
    print(f"Structure setup completed in {count} POSCAR directories")
    return 0


def _run_concat(paths: ArtemisPaths, dry_run: bool) -> int:
    count = 0
    for poscar_dir in _iter_poscar_dirs(paths.interfaces_dir):
        if not (poscar_dir / "concat_artemis_data.py").exists():
            continue
        rc = _run_command(["python3", "concat_artemis_data.py"], poscar_dir, dry_run)
        if rc == 0:
            count += 1
    print(f"Concat run completed in {count} POSCAR directories")
    return 0


def _submit(paths: ArtemisPaths, run_root: Path, sbatch_name: str, output_suffix: str, dry_run: bool) -> int:
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    destination = run_root / timestamp
    print(f"Resolved submit destination: {destination}")

    if dry_run:
        print(f"[DRY-RUN] copytree {paths.interfaces_dir} -> {destination}")
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(paths.interfaces_dir, destination, dirs_exist_ok=True)

    submitted = 0
    for poscar_dir in _iter_poscar_dirs(destination if not dry_run else paths.interfaces_dir):
        files = set(p.name for p in poscar_dir.iterdir())
        if any(name.endswith(output_suffix) for name in files):
            continue
        sbatch_file = poscar_dir / sbatch_name
        if not sbatch_file.exists():
            continue
        rc = _run_command(["sbatch", sbatch_name], poscar_dir, dry_run)
        if rc != 0 and not dry_run:
            return rc
        submitted += 1
    print(f"Submission scan completed: {submitted} jobs")
    return 0


def _cleanup(paths: ArtemisPaths, dry_run: bool, remove_mlp: bool, remove_dft: bool, remove_structure: bool) -> int:
    for poscar_dir in _iter_poscar_dirs(paths.interfaces_dir):
        if remove_mlp:
            _delete_file(poscar_dir / "SBATCH_MLP_MACE", dry_run)
            _delete_file(poscar_dir / "run-mace.py", dry_run)
        if remove_dft:
            _delete_file(poscar_dir / "SBATCH_DFT_VASP", dry_run)
        if remove_structure:
            _delete_file(poscar_dir / "concat_artemis_data.py", dry_run)
            _delete_file(poscar_dir / "poscar_data.json", dry_run)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="MP1 Artemis Utilities Manager")
    parser.add_argument("--root", default=None, help="Root directory for Artemis workflows")
    parser.add_argument("--interfaces-dir", default=None, help="Override interfaces directory")
    parser.add_argument("--run-root-mlp", default=None, help="Override MLP run root")
    parser.add_argument("--run-root-dft", default=None, help="Override DFT run root")
    parser.add_argument("--dry-run", action="store_true", help="Print actions without mutating files or submitting")

    parser.add_argument("--setup-mlp-mace-jobs", action="store_true")
    parser.add_argument("--setup-dft-vasp-jobs", action="store_true")
    parser.add_argument("--setup-jobs", action="store_true")
    parser.add_argument("--setup-structure", action="store_true")
    parser.add_argument("--setup-all", action="store_true")

    parser.add_argument("--run-concat", action="store_true")
    parser.add_argument("--submit-mlp-mace-jobs", action="store_true")
    parser.add_argument("--submit-dft-vasp-jobs", action="store_true")
    parser.add_argument("--submit-jobs", action="store_true")

    parser.add_argument("--cleanup-structure", action="store_true")
    parser.add_argument("--cleanup-mlp-mace-jobs", action="store_true")
    parser.add_argument("--cleanup-dft-vasp-jobs", action="store_true")
    parser.add_argument("--cleanup-jobs", action="store_true")
    parser.add_argument("--cleanup-all", action="store_true")

    parser.add_argument("--all", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    paths = _resolve_paths(args)
    _print_paths(paths)

    if not _require_safe_interfaces(paths.interfaces_dir):
        return 1

    dry_run = bool(args.dry_run)
    selected = any(vars(args).values())

    if args.setup_mlp_mace_jobs:
        _setup_mlp(paths, dry_run)
    if args.setup_dft_vasp_jobs:
        _setup_dft(paths, dry_run)
    if args.setup_jobs:
        _setup_mlp(paths, dry_run)
        _setup_dft(paths, dry_run)
    if args.setup_structure:
        _setup_structure(paths, dry_run)
    if args.setup_all:
        _setup_mlp(paths, dry_run)
        _setup_dft(paths, dry_run)
        _setup_structure(paths, dry_run)

    if args.run_concat:
        _run_concat(paths, dry_run)

    if args.submit_mlp_mace_jobs:
        _submit(paths, paths.run_root_mlp, "SBATCH_MLP_MACE", "MLPMACE.out", dry_run)
    if args.submit_dft_vasp_jobs:
        _submit(paths, paths.run_root_dft, "SBATCH_DFT_VASP", "DFTVasp.out", dry_run)
    if args.submit_jobs:
        _submit(paths, paths.run_root_mlp, "SBATCH_MLP_MACE", "MLPMACE.out", dry_run)
        _submit(paths, paths.run_root_dft, "SBATCH_DFT_VASP", "DFTVasp.out", dry_run)

    if args.cleanup_mlp_mace_jobs:
        _cleanup(paths, dry_run, remove_mlp=True, remove_dft=False, remove_structure=False)
    if args.cleanup_dft_vasp_jobs:
        _cleanup(paths, dry_run, remove_mlp=False, remove_dft=True, remove_structure=False)
    if args.cleanup_jobs:
        _cleanup(paths, dry_run, remove_mlp=True, remove_dft=True, remove_structure=False)
    if args.cleanup_structure:
        _cleanup(paths, dry_run, remove_mlp=False, remove_dft=False, remove_structure=True)
    if args.cleanup_all:
        _cleanup(paths, dry_run, remove_mlp=True, remove_dft=True, remove_structure=True)

    if args.all:
        _cleanup(paths, dry_run, remove_mlp=True, remove_dft=True, remove_structure=True)
        _setup_mlp(paths, dry_run)
        _setup_dft(paths, dry_run)
        _setup_structure(paths, dry_run)
        _run_concat(paths, dry_run)

    if not selected:
        parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
