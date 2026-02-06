from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def _env_with_src() -> dict[str, str]:
    env = dict(os.environ)
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(SRC) if not existing else f"{SRC}:{existing}"
    return env


def test_import_modules() -> None:
    env = _env_with_src()
    result = subprocess.run(
        [sys.executable, "-c", "import mlp_tools.id, mlp_tools.energies, mlp_tools.vasp"],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_mp1_id_help() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "mlp_tools.id", "--help"],
        cwd=ROOT,
        env=_env_with_src(),
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "ID Utilities Manager" in result.stdout


def test_mp1_energies_help() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "mlp_tools.energies", "--help"],
        cwd=ROOT,
        env=_env_with_src(),
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "Compute formation, surface, and interface energies" in result.stdout


def test_mp1_vasp_help() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "mlp_tools.vasp", "--help"],
        cwd=ROOT,
        env=_env_with_src(),
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "MP1 VASP Utilities Manager" in result.stdout
