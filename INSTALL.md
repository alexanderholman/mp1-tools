# Installation

## Prerequisites

- Python 3.9+
- `pip`

## From scratch

Install from the repository root:

```bash
python -m pip install .
```

For development work, use editable mode:

```bash
python -m pip install -e .
```

## Usage

```bash
mp1-id --help
mp1-energies --help
mp1-vasp --help
mp1-artemis --help
```

For `mp1-vasp`, configure POTCAR location with one of:

```bash
export VASP_POTCAR_DIR="/path/to/POTCARS"
```

or

```bash
mp1-vasp --generate-potcar --workdir /path/to/calc --potcar-dir /path/to/POTCARS
```

For `mp1-artemis`, run with `--dry-run` first to verify resolved paths and actions:

```bash
mp1-artemis --root /path/to/workspace --setup-all --dry-run
```

## Updating

Pull latest changes and reinstall:

```bash
python -m pip install -e .
```

## Removal and how-to uninstall

```bash
python -m pip uninstall mp1-tools
```
