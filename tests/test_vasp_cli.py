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


def test_generate_potcar_requires_configured_potcar_dir(tmp_path: Path) -> None:
    poscar = tmp_path / "POSCAR"
    poscar.write_text(
        "Si\n"
        "1.0\n"
        "5 0 0\n"
        "0 5 0\n"
        "0 0 5\n"
        "Si\n"
        "2\n"
        "Direct\n"
        "0 0 0\n"
        "0.5 0.5 0.5\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [sys.executable, "-m", "mp1_tools.vasp", "--generate-potcar", "--workdir", str(tmp_path)],
        cwd=ROOT,
        env=_env_with_src(),
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 1
    assert "POTCAR directory required" in result.stdout
