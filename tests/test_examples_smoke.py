"""Smoke tests: run each examples/*.py end-to-end and require a clean exit.

Examples double as regression tests for the public workflow. Plotting is
forced to a non-interactive backend (Agg) so plt.show() cannot block.

Split into three groups by current status:
  - working single-compartment examples (must pass now),
  - 1D diffusion examples (must pass; Bug 2 is fixed),
  - membrane_slab (3D slab + bulk, ~50 s; BLAS threads capped at 8).
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"

WORKING = ["degradation_ode.py", "birth_death.py", "Gillespie_AB_system.py"]
BROKEN_BUG2 = ["1D_diffusion.py", "1D_reac_diffusion.py", "1D_reac_diff_AB.py"]


def _run_example(name):
    env = dict(os.environ, MPLBACKEND="Agg", OPENBLAS_NUM_THREADS="8", OMP_NUM_THREADS="8")
    return subprocess.run(
        [sys.executable, str(EXAMPLES_DIR / name)],
        env=env, capture_output=True, text=True, timeout=300,
    )


@pytest.mark.parametrize("name", WORKING)
def test_working_examples_run(name):
    """Single-compartment ODE/Gillespie examples run to completion."""
    result = _run_example(name)
    assert result.returncode == 0, result.stderr[-2000:]


@pytest.mark.parametrize("name", BROKEN_BUG2)
def test_diffusion_examples_run(name):
    """1D diffusion examples run (ODE inter-compartment diffusion fixed in Step 2)."""
    result = _run_example(name)
    assert result.returncode == 0, result.stderr[-2000:]


def test_membrane_slab_runs():
    """membrane_slab.py (ODESystem, join3D with an anisotropic link, mol/L set_q) runs and conserves drug."""
    result = _run_example("membrane_slab.py")
    assert result.returncode == 0, result.stderr[-2000:]
    line = next(l for l in result.stdout.splitlines() if l.startswith("drug molecules"))
    n0, n1 = (float(x.split()[-1]) for x in line.split(":", 1)[1].split(","))
    assert n1 == pytest.approx(n0, rel=1e-6)