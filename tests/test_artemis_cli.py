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


def test_artemis_dry_run_requires_safe_interfaces(tmp_path: Path) -> None:
    root = tmp_path / "workspace"
    root.mkdir()
    interfaces = root / "DINTERFACES"
    interfaces.mkdir()
    poscar_dir = interfaces / "test" / "case"
    poscar_dir.mkdir(parents=True)
    (poscar_dir / "POSCAR").write_text("placeholder", encoding="utf-8")
    (poscar_dir.parent / "struc_dat.txt").write_text(
        "Lower crystal Miller plane: 1 1 1\nUpper crystal Miller plane: 1 1 2\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [sys.executable, "-m", "mlp_tools.artemis", "--root", str(root), "--setup-mlp-mace-jobs", "--dry-run"],
        cwd=ROOT,
        env=_env_with_src(),
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "[DRY-RUN] write" in result.stdout
    assert not (poscar_dir / "SBATCH_MLP_MACE").exists()


def test_artemis_setup_mlp_renders_job_name_without_placeholder_error(tmp_path: Path) -> None:
    root = tmp_path / "workspace"
    root.mkdir()
    interfaces = root / "DINTERFACES"
    interfaces.mkdir()
    poscar_dir = interfaces / "test" / "case"
    poscar_dir.mkdir(parents=True)
    (poscar_dir / "POSCAR").write_text("placeholder", encoding="utf-8")
    (poscar_dir.parent / "struc_dat.txt").write_text(
        "Lower crystal Miller plane: 1 1 1\nUpper crystal Miller plane: 1 1 2\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [sys.executable, "-m", "mlp_tools.artemis", "--root", str(root), "--setup-mlp-mace-jobs"],
        cwd=ROOT,
        env=_env_with_src(),
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    sbatch = (poscar_dir / "SBATCH_MLP_MACE").read_text(encoding="utf-8")
    assert "--job-name=111-112MLPMACE" in sbatch
