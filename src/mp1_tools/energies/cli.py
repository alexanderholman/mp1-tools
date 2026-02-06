"""Compute bulk/formation/surface/interface energies from VASP outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from typing import Any, Dict, Optional


def _import_pymatgen():
    try:
        from pymatgen.io.vasp.inputs import Kpoints, Poscar
        from pymatgen.io.vasp.outputs import Outcar, Vasprun
    except ImportError as exc:
        raise RuntimeError(
            "Missing dependency for mp1-energies. Install with: pip install pymatgen"
        ) from exc
    return Poscar, Kpoints, Vasprun, Outcar


def find_file(path: str, candidates):
    if os.path.isfile(path):
        return path

    for name in candidates:
        candidate_path = os.path.join(path, name)
        if os.path.isfile(candidate_path):
            return candidate_path
    raise FileNotFoundError(f"Could not find any of {candidates} under {path}")


def read_final_energy(calc_path: str) -> float:
    _poscar, _kpoints, Vasprun, Outcar = _import_pymatgen()
    try:
        vr_file = find_file(calc_path, ["vasprun.xml", "vasprun.xml.gz"])
        vr = Vasprun(vr_file)
        return float(vr.final_energy)
    except Exception:
        outcar_file = find_file(calc_path, ["OUTCAR"])
        outcar = Outcar(outcar_file)
        return float(outcar.final_energy)


def read_structure(calc_path: str):
    Poscar, _kpoints, _vasprun, _outcar = _import_pymatgen()
    poscar_file = find_file(calc_path, ["CONTCAR", "POSCAR"])
    poscar = Poscar.from_file(poscar_file)
    return poscar.structure


def load_db(path: str) -> dict:
    if path is None or not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def save_db(db: dict, path: str):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(db, handle, indent=2, sort_keys=True)


def ensure_label_not_in_db(db: dict, label: str):
    if label in db:
        raise ValueError(
            f"Label '{label}' already exists in DB. Choose a different label or remove it manually."
        )


def hash_file(path: str) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def get_vasp_setup(calc_path: str) -> Dict[str, Any]:
    _poscar, Kpoints, _vasprun, _outcar = _import_pymatgen()

    if os.path.isfile(calc_path):
        root = os.path.dirname(os.path.abspath(calc_path))
    else:
        root = os.path.abspath(calc_path)

    setup: Dict[str, Any] = {"calc_root": root}
    files = {"incar": "INCAR", "kpoints": "KPOINTS", "potcar": "POTCAR"}
    component_hashes = []

    for key, fname in files.items():
        fpath = os.path.join(root, fname)
        if os.path.isfile(fpath):
            setup[f"{key}_path"] = fpath
            sha = hash_file(fpath)
            setup[f"{key}_sha1"] = sha
            component_hashes.append(sha)

            if key == "kpoints":
                try:
                    kp = Kpoints.from_file(fpath)
                    style = getattr(kp.style, "name", str(kp.style))
                    setup["kpoints_style"] = style
                    if kp.kpts:
                        setup["kpts_grid"] = kp.kpts[0]
                    setup["num_kpts"] = kp.num_kpts
                except Exception:
                    pass

    if component_hashes:
        setup["setup_id"] = hashlib.sha1("".join(component_hashes).encode("utf-8")).hexdigest()
    else:
        setup["setup_id"] = None

    return setup


def compute_scaled_bulk_energy(bulk_energy: float, bulk_natoms: int, target_natoms: int) -> float:
    if bulk_natoms <= 0:
        raise ValueError("bulk_natoms must be positive.")
    return bulk_energy * (target_natoms / bulk_natoms)


def cmd_bulk(args):
    db = load_db(args.output)
    ensure_label_not_in_db(db, args.label)

    calc_path = os.path.abspath(args.calc)
    energy = read_final_energy(calc_path)
    structure = read_structure(calc_path)
    natoms = len(structure)

    db[args.label] = {
        "type": "bulk",
        "calc_path": calc_path,
        "natoms": natoms,
        "energy_eV": energy,
        "energy_per_atom_eV": energy / natoms,
        "vasp_setup": get_vasp_setup(calc_path),
    }
    save_db(db, args.output)


def cmd_formation(args):
    db = load_db(args.output)
    ensure_label_not_in_db(db, args.label)

    calc_path = os.path.abspath(args.calc)
    energy = read_final_energy(calc_path)
    natoms = len(read_structure(calc_path))

    bulk_path = os.path.abspath(args.bulk)
    bulk_energy = read_final_energy(bulk_path)
    bulk_natoms = len(read_structure(bulk_path))

    scaled_bulk_energy = compute_scaled_bulk_energy(bulk_energy, bulk_natoms, natoms)
    formation_energy = energy - scaled_bulk_energy

    db[args.label] = {
        "type": "formation",
        "calc_path": calc_path,
        "bulk_path": bulk_path,
        "natoms": natoms,
        "bulk_natoms": bulk_natoms,
        "energy_eV": energy,
        "bulk_energy_eV": bulk_energy,
        "scaled_bulk_energy_eV": scaled_bulk_energy,
        "formation_energy_eV": formation_energy,
        "formation_energy_per_atom_eV": formation_energy / natoms,
        "vasp_setup": get_vasp_setup(calc_path),
        "bulk_vasp_setup": get_vasp_setup(bulk_path),
    }
    save_db(db, args.output)


def cmd_surface(args):
    db = load_db(args.output)
    ensure_label_not_in_db(db, args.label)

    slab_path = os.path.abspath(args.calc)
    energy_slab = read_final_energy(slab_path)
    natoms_slab = len(read_structure(slab_path))

    bulk_path = os.path.abspath(args.bulk)
    bulk_energy = read_final_energy(bulk_path)
    natoms_bulk = len(read_structure(bulk_path))

    scaled_bulk_energy = compute_scaled_bulk_energy(bulk_energy, natoms_bulk, natoms_slab)
    formation_energy = energy_slab - scaled_bulk_energy

    if args.nsurfaces <= 0:
        raise ValueError("Number of surfaces must be positive.")
    gamma = formation_energy / (args.nsurfaces * args.area)

    db[args.label] = {
        "type": "surface",
        "calc_path": slab_path,
        "bulk_path": bulk_path,
        "natoms": natoms_slab,
        "bulk_natoms": natoms_bulk,
        "energy_eV": energy_slab,
        "bulk_energy_eV": bulk_energy,
        "scaled_bulk_energy_eV": scaled_bulk_energy,
        "formation_energy_eV": formation_energy,
        "area_A2": args.area,
        "n_surfaces": args.nsurfaces,
        "surface_energy_eV_per_A2": gamma,
        "vasp_setup": get_vasp_setup(slab_path),
        "bulk_vasp_setup": get_vasp_setup(bulk_path),
    }
    save_db(db, args.output)


def cmd_interface(args):
    db = load_db(args.output)
    ensure_label_not_in_db(db, args.label)

    int_path = os.path.abspath(args.calc)
    energy_int = read_final_energy(int_path)
    natoms_int = len(read_structure(int_path))

    bulk_a_path = os.path.abspath(args.bulk_a)
    energy_bulk_a = read_final_energy(bulk_a_path)
    natoms_bulk_a = len(read_structure(bulk_a_path))

    bulk_b_path = os.path.abspath(args.bulk_b)
    energy_bulk_b = read_final_energy(bulk_b_path)
    natoms_bulk_b = len(read_structure(bulk_b_path))

    if args.natoms_a + args.natoms_b != natoms_int:
        raise ValueError(
            f"n_atoms_a + n_atoms_b ({args.natoms_a + args.natoms_b}) != natoms_interface ({natoms_int})"
        )

    scaled_bulk_a = compute_scaled_bulk_energy(energy_bulk_a, natoms_bulk_a, args.natoms_a)
    scaled_bulk_b = compute_scaled_bulk_energy(energy_bulk_b, natoms_bulk_b, args.natoms_b)
    ref_energy = scaled_bulk_a + scaled_bulk_b

    formation_energy = energy_int - ref_energy
    gamma_int = formation_energy / args.area

    work_adhesion: Optional[float] = None
    if args.surface_energy_a is not None and args.surface_energy_b is not None:
        work_adhesion = args.surface_energy_a + args.surface_energy_b - gamma_int

    db[args.label] = {
        "type": "interface",
        "calc_path": int_path,
        "bulk_a_path": bulk_a_path,
        "bulk_b_path": bulk_b_path,
        "natoms_interface": natoms_int,
        "natoms_a": args.natoms_a,
        "natoms_b": args.natoms_b,
        "energy_interface_eV": energy_int,
        "bulk_a_energy_eV": energy_bulk_a,
        "bulk_b_energy_eV": energy_bulk_b,
        "scaled_bulk_a_energy_eV": scaled_bulk_a,
        "scaled_bulk_b_energy_eV": scaled_bulk_b,
        "reference_energy_eV": ref_energy,
        "formation_energy_eV": formation_energy,
        "interface_energy_eV_per_A2": gamma_int,
        "area_A2": args.area,
        "surface_energy_a_eV_per_A2": args.surface_energy_a,
        "surface_energy_b_eV_per_A2": args.surface_energy_b,
        "work_of_adhesion_eV_per_A2": work_adhesion,
        "vasp_setup": get_vasp_setup(int_path),
        "bulk_a_vasp_setup": get_vasp_setup(bulk_a_path),
        "bulk_b_vasp_setup": get_vasp_setup(bulk_b_path),
    }
    save_db(db, args.output)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute formation, surface, and interface energies from VASP runs "
            "and store in a JSON DB."
        )
    )
    parser.add_argument(
        "-o",
        "--output",
        default="energies.json",
        help="Output JSON file to store results (default: energies.json)",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_bulk = subparsers.add_parser("bulk", help="Register a bulk calculation (energy per atom).")
    p_bulk.add_argument("label", help="Label/key for this bulk entry in the JSON DB.")
    p_bulk.add_argument("calc", help="Path to bulk VASP calculation directory (or file).")
    p_bulk.set_defaults(func=cmd_bulk)

    p_form = subparsers.add_parser(
        "formation",
        help="Formation energy of a structure with respect to a bulk reference.",
    )
    p_form.add_argument("label", help="Label/key for this formation entry.")
    p_form.add_argument("calc", help="Path to structure VASP calc directory (or file).")
    p_form.add_argument("bulk", help="Path to bulk reference VASP calc directory (or file).")
    p_form.set_defaults(func=cmd_formation)

    p_surf = subparsers.add_parser(
        "surface",
        help="Surface energy of a slab with respect to a bulk reference.",
    )
    p_surf.add_argument("label", help="Label/key for this surface entry.")
    p_surf.add_argument("calc", help="Path to slab VASP calc directory (or file).")
    p_surf.add_argument("bulk", help="Path to bulk reference VASP calc directory (or file).")
    p_surf.add_argument("--area", type=float, required=True, help="Surface area in A^2.")
    p_surf.add_argument(
        "--nsurfaces",
        type=int,
        default=2,
        help="Number of equivalent surfaces in the slab (default: 2).",
    )
    p_surf.set_defaults(func=cmd_surface)

    p_int = subparsers.add_parser(
        "interface",
        help="Interface energy between A and B and (optionally) work of adhesion.",
    )
    p_int.add_argument("label", help="Label/key for this interface entry.")
    p_int.add_argument("calc", help="Path to interface VASP calc directory (or file).")
    p_int.add_argument("bulk_a", help="Path to bulk A VASP calc directory (or file).")
    p_int.add_argument("bulk_b", help="Path to bulk B VASP calc directory (or file).")
    p_int.add_argument("--natoms-a", type=int, required=True, help="Number of A atoms.")
    p_int.add_argument("--natoms-b", type=int, required=True, help="Number of B atoms.")
    p_int.add_argument("--area", type=float, required=True, help="Interface area in A^2.")
    p_int.add_argument("--surface-energy-a", type=float, help="Surface energy of A in eV/A^2.")
    p_int.add_argument("--surface-energy-b", type=float, help="Surface energy of B in eV/A^2.")
    p_int.set_defaults(func=cmd_interface)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
