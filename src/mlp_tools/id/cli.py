"""Layer ID plotting utility."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


DEFAULT_COLORS = {
    "H": "#CCCCCC",
    "C": "#222222",
    "N": "#3050F8",
    "O": "#FF0D0D",
    "F": "#90E050",
    "Si": "#F0C8A0",
    "Ge": "#6689CC",
    "Sn": "#668080",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="ID Utilities Manager")

    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=Path.cwd() / "POSCAR",
        help="Path to the input POSCAR. Default = ./POSCAR",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output directory. Default = <poscar>_layers/",
    )
    parser.add_argument(
        "--round",
        type=int,
        default=6,
        help="Number of decimal places for fractional coordinate rounding.",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=0.0,
        help="Layer width in fractional coordinates. Default = 0.0 (exact).",
    )
    parser.add_argument(
        "--range-mode",
        choices=["auto", "zero-one", "set"],
        default="auto",
        help="Global axis range: auto, zero-one, or set via --xrange/--yrange.",
    )
    parser.add_argument(
        "--per-axis-range-mode",
        nargs=3,
        default=None,
        help="Per-axis overrides for a b c, e.g. auto zero-one auto",
    )
    parser.add_argument(
        "--xrange",
        type=float,
        nargs=2,
        default=None,
        help="X-axis range: xmin xmax (requires range-mode=set)",
    )
    parser.add_argument(
        "--yrange",
        type=float,
        nargs=2,
        default=None,
        help="Y-axis range: ymin ymax (requires range-mode=set)",
    )
    parser.add_argument(
        "--padding",
        type=float,
        default=0.05,
        help="Padding for auto range mode (fraction of axis). Default = 0.05",
    )
    parser.add_argument("--square", action="store_true", help="Force square aspect ratio")
    parser.add_argument("--aspect", type=float, default=None, help="Force a custom aspect ratio")
    parser.add_argument(
        "--color-map",
        type=str,
        nargs="*",
        default=None,
        help="Override element colors, e.g. Si=#FF0000 Ge=#00FF00",
    )
    parser.add_argument(
        "--radius-scale",
        type=float,
        default=60.0,
        help="Scale factor for atomic radius in plotting",
    )
    return parser


def _import_dependencies():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from ase.io import read
    except ImportError as exc:
        raise RuntimeError(
            "Missing dependency for mp1-id. Install with: pip install ase matplotlib numpy"
        ) from exc
    return read, np, plt


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.input.exists():
        print(f"Input file '{args.input}' not found.")
        return 1

    read, _np, plt = _import_dependencies()

    if args.output is None:
        outdir = args.input.parent / f"{args.input.stem}_layers"
    else:
        outdir = args.output

    for axis in ["a", "b", "c"]:
        (outdir / axis).mkdir(parents=True, exist_ok=True)

    element_colors = DEFAULT_COLORS.copy()
    if args.color_map:
        for pair in args.color_map:
            if "=" in pair:
                el, col = pair.split("=", maxsplit=1)
                element_colors[el.strip()] = col.strip()

    atoms = read(str(args.input))
    scaled = atoms.get_scaled_positions()
    cart = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()

    layer_keys = {"a": {}, "b": {}, "c": {}}

    def layer_key(value: float):
        if args.width > 0:
            return round(value / args.width)
        return round(value, args.round)

    for idx, frac in enumerate(scaled):
        fa, fb, fc = frac
        layer_keys["a"].setdefault(layer_key(fa), []).append(idx)
        layer_keys["b"].setdefault(layer_key(fb), []).append(idx)
        layer_keys["c"].setdefault(layer_key(fc), []).append(idx)

    axis_layer_index = {"a": {}, "b": {}, "c": {}}
    for axis in ["a", "b", "c"]:
        for n, lkey in enumerate(sorted(layer_keys[axis].keys()), start=1):
            axis_layer_index[axis][lkey] = n

    axis_layer_inlayer = {"a": {}, "b": {}, "c": {}}
    for axis in ["a", "b", "c"]:
        axis_layer_inlayer[axis] = {
            lkey: {atom_idx: i + 1 for i, atom_idx in enumerate(atom_list)}
            for lkey, atom_list in layer_keys[axis].items()
        }

    atom_ids = []
    for idx, sym in enumerate(symbols, start=1):
        atom_ids.append(f"{sym.lower()} {idx}")

    per_axis_modes = {}
    if args.per_axis_range_mode:
        pa = args.per_axis_range_mode
        per_axis_modes = {"a": pa[0], "b": pa[1], "c": pa[2]}

    def get_range_mode(axis: str) -> str:
        if axis in per_axis_modes:
            return per_axis_modes[axis]
        return args.range_mode

    axis_global_ranges = {}
    for axis in ["a", "b", "c"]:
        if axis == "a":
            xs_all = scaled[:, 1]
            ys_all = scaled[:, 2]
        elif axis == "b":
            xs_all = scaled[:, 0]
            ys_all = scaled[:, 2]
        else:
            xs_all = scaled[:, 0]
            ys_all = scaled[:, 1]

        xmin, xmax = xs_all.min(), xs_all.max()
        ymin, ymax = ys_all.min(), ys_all.max()

        dx = xmax - xmin
        dy = ymax - ymin

        xmin -= dx * args.padding
        xmax += dx * args.padding
        ymin -= dy * args.padding
        ymax += dy * args.padding

        axis_global_ranges[axis] = (xmin, xmax, ymin, ymax)

    def plot_layer(axis: str, lkey, atom_indices: list[int]):
        coords = scaled[atom_indices]

        if axis == "a":
            xs, ys = coords[:, 1], coords[:, 2]
            xlabel, ylabel = "Fractional b", "Fractional c"
        elif axis == "b":
            xs, ys = coords[:, 0], coords[:, 2]
            xlabel, ylabel = "Fractional a", "Fractional c"
        else:
            xs, ys = coords[:, 0], coords[:, 1]
            xlabel, ylabel = "Fractional a", "Fractional b"

        lnum = axis_layer_index[axis][lkey]

        plt.figure(figsize=(6, 6))
        colors = [element_colors.get(symbols[a], "#000000") for a in atom_indices]
        sizes = [args.radius_scale * (atoms[a].number**0.3) for a in atom_indices]
        plt.scatter(xs, ys, s=sizes, c=colors)

        for coord_idx, atom_idx in enumerate(atom_indices):
            plt.text(xs[coord_idx], ys[coord_idx], atom_ids[atom_idx], fontsize=5)

        mode = get_range_mode(axis)
        if mode == "auto":
            xmin, xmax, ymin, ymax = axis_global_ranges[axis]
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
        elif mode == "zero-one":
            plt.xlim(0, 1)
            plt.ylim(0, 1)
        elif mode == "set":
            if args.xrange is None or args.yrange is None:
                raise ValueError("--xrange and --yrange required for --range-mode set")
            plt.xlim(*args.xrange)
            plt.ylim(*args.yrange)

        if args.square:
            plt.gca().set_aspect("equal", adjustable="box")
        elif args.aspect is not None:
            plt.gca().set_aspect(args.aspect)

        plt.title(f"Layer {axis.upper()}_{lnum} (key={lkey})")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()

        outfile = outdir / axis / f"{axis}_layer_{lnum}.png"
        plt.savefig(outfile, dpi=300)
        plt.close()

    def write_json(axis: str, lkey, atom_indices: list[int]):
        lnum = axis_layer_index[axis][lkey]
        data = {
            "axis": axis,
            "layer_index": lnum,
            "layer_key": lkey,
            "num_atoms": len(atom_indices),
            "axes": {
                "projection": {
                    "x": "b" if axis == "a" else "a",
                    "y": "c" if axis in ["a", "b"] else "b",
                },
                "range_mode": get_range_mode(axis),
                "global_range": {
                    "xmin": axis_global_ranges[axis][0],
                    "xmax": axis_global_ranges[axis][1],
                    "ymin": axis_global_ranges[axis][2],
                    "ymax": axis_global_ranges[axis][3],
                },
            },
            "atoms": [],
        }
        for ai in atom_indices:
            data["atoms"].append(
                {
                    "index": ai,
                    "symbol": symbols[ai],
                    "id": atom_ids[ai],
                    "frac": scaled[ai].tolist(),
                    "cartesian": cart[ai].tolist(),
                }
            )

        out = outdir / axis / f"{axis}_layer_{lnum}.json"
        with out.open("w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2)

    print("Generating images + JSON exports...")
    for axis in ["a", "b", "c"]:
        for lkey, atom_idxs in sorted(layer_keys[axis].items()):
            plot_layer(axis, lkey, atom_idxs)
            write_json(axis, lkey, atom_idxs)

    print(f"Output saved in: {outdir}")
    print("\nFull Atom ID Table:")
    for i, atom_id in enumerate(atom_ids, start=1):
        print(f"{i:3d}  {atom_id}")

    print("\nLayer summary:")
    for axis in ["a", "b", "c"]:
        print(f"\nAxis {axis.upper()}:")
        for lkey, atom_list in sorted(layer_keys[axis].items()):
            print(f"  Layer {axis_layer_index[axis][lkey]:2d} (key={lkey}) -> {atom_list}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
