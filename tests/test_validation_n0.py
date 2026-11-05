"""N0 physics validation gates (docs/TASK.md, N0).  They read data/n0/validation_summary.json,
which scripts/n0_validation.py writes from the runs of configs/n0/ (python -m sbaniso.run
configs/n0/*.yaml; about 20 min on 8 cores).  Skipped when the summary does not exist.

    pytest -m validation
"""
import json
from pathlib import Path

import pytest

SUMMARY = Path(__file__).resolve().parents[1] / "data" / "n0" / "validation_summary.json"
pytestmark = [pytest.mark.slow, pytest.mark.validation,
              pytest.mark.skipif(not SUMMARY.exists(), reason="run scripts/n0_validation.py first")]


@pytest.fixture(scope="module")
def S():
    return json.loads(SUMMARY.read_text())


@pytest.mark.parametrize("fam", ["xy", "r0_kt", "j2zero_kt"])
def test_kt_weber_minnhagen(S, fam):
    w = S[fam]["weber_minnhagen"]
    assert abs(w["T"] - S[fam]["T_KT_ref"]) <= 0.01, w


@pytest.mark.parametrize("fam", ["xy", "r0_kt", "j2zero_kt"])
def test_kt_finite_size_helicity_matches_hasenbusch(S, fam):
    for r in S[fam]["hasenbusch"]:
        if "z" in r:
            assert abs(r["z"]) < 4, r


@pytest.mark.parametrize("fam", ["xy", "r0_kt", "j2zero_kt"])
def test_kt_crossings_decrease_toward_T_kt(S, fam):
    Tl = [c["T"] for c in S[fam]["crossings"]]
    assert all(b < a for a, b in zip(Tl, Tl[1:])) and Tl[-1] > S[fam]["T_KT_ref"], Tl


@pytest.mark.parametrize("fam", ["ising", "r0_ising", "j2zero_ising"])
def test_ising_binder_crossing(S, fam):
    b = S[fam]["binder"][-1]
    assert abs(b["T"] - S[fam]["T_ref"]) <= 0.01, b
    assert abs(b["U"] - 0.611) <= 0.01, b


@pytest.mark.parametrize("key", ["3.8", "4.1", "4.4", "4.1_random_vs_ordered"])
def test_pilot_regression(S, key):
    t = S["regression"][key]["total"]
    # 70 z-scores of correlated observables: no |z| > 4 and chi^2/dof compatible with 1
    assert t["max_abs_z"] < 4 and t["chi2dof"] < 2.0, t


def test_equilibration_diagnostics(S):
    diags = [d for k, v in S.items() if k != "regression" for d in v["diagnostics"]]
    diags += [v["diagnostics"] for v in S["regression"].values()]
    for d in diags:   # combined criterion of analysis.equilibration_check
        assert d["roundtrips"] >= 20, d
        assert d["swap_acc_min"] >= 0.2, d
        assert d["quarter_z_max"] < d["quarter_threshold"], d   # Bonferroni, p = 0.01 per run
        assert max(d["tau_max"].values()) <= d["nmeas"] / 200, d     # 20 blocks of >= 10 tau
