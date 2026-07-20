"""Test the offline D-residue survey reproduction script."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_reproduce_d_residue_errors_script():
    proc = subprocess.run(
        [sys.executable, str(ROOT / "benchmarks" / "reproduce_d_residue_errors.py")],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "PASS" in proc.stdout
    assert "12,573" in proc.stdout or "12573" in proc.stdout.replace(",", "")
    assert "29" in proc.stdout
