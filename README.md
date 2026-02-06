# mp1-tools

`mp1-tools` is a Python package with command-line tools for MP1 workflows:

- `mp1-id`: layer ID plotting and JSON export utility for structure files.
- `mp1-energies`: bulk, formation, surface, and interface energy utilities for VASP runs.
- `mp1-vasp`: setup helper for generating `POTCAR`, `nbands.txt`, `nelect.txt`, `INCAR`, and `KPOINTS`.

## Quick start

Install in editable mode:

```bash
python -m pip install -e .
```

Show command help:

```bash
mp1-id --help
mp1-energies --help
mp1-vasp --help
```

`mp1-vasp` requires a configured POTCAR directory via `--potcar-dir` or `VASP_POTCAR_DIR`.

Current status is tracked in `WIP.md` and `tasks/wip.md`.

See [INSTALL.md](INSTALL.md) for details.
