"""Negative and positive control of the quarter-split equilibration check (owner's review of N0,
2026-11-05).  Near the crossing (r = 4.1, L = 64) the run WITHOUT the psi-frame sigma cluster
(flags = 1, 0, 1) is known not to equilibrate in 1000 + 10000 sweeps (pilot v1: Ising tau of
hundreds to 1500 sweeps, drifting |m|); the check must flag it.  The same run with all three cluster
types must pass.  Runs are made on first use (~3 min each on 8 cores) and kept in runs/controls/.

    pytest -m slow tests/test_equilibration_control.py
"""
from pathlib import Path

import pytest

from sbaniso.analysis import equilibration_check, load_run, quarter_check
from sbaniso.run import run_config

REPO = Path(__file__).resolve().parents[1]
BASE = dict(model=dict(r=4.1), L=64, temperatures=dict(Tmin=6.3, Tmax=8.0, n=26), ntherm=1000,
            nmeas=10000, seed=5101, checkpoint_every=2000)
pytestmark = pytest.mark.slow


def _run(name, flags):
    return load_run(run_config(dict(BASE, name=name, clusters=flags,
                                    output=f"runs/controls/{name}_r4.1_L64.h5"),
                               verbose=False))


def test_equilibration_check_flags_run_without_sigma_psi_cluster():
    """The combined criterion must reject the negative control, and the block-length (tau) part
    must reject it on its own: that is the part added because the quarter test alone did not
    (z = 4.22 at a threshold of 4.22, docs/N0_report.md appendix)."""
    r = equilibration_check(_run("negctl_nopsi", [1, 0, 1]))
    assert not r["pass"], r
    assert "tau" in r["reasons"], r


def test_equilibration_check_passes_same_run_with_all_clusters():
    run = _run("posctl_all", [1, 1, 1])
    r = equilibration_check(run)
    assert r["pass"], r
    z, d = quarter_check(run)
    assert z.max() < d["threshold"]
