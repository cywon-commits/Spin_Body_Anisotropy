"""Single-body check of the axis-selection coefficient Gamma (theory_summary P2), for Fig. 2.

Recomputes the table of docs/N3_report.md section 1 (fixed-orientation relaxation vs formula) and
samples E(n) over random body orientations at the production anisotropy.  Needs the polyspin
anisotropy extension (pip install -e .[fullmodel], tag v0.2.0-aniso).

    python scripts/n4_gamma_check.py   -> data/n4/gamma_check.json
"""
import json
from pathlib import Path

import numpy as np
from polyspin.groundstate import gamma_formula, minimize_cluster_at_orientation, oriented_tetrahedron
from polyspin.model import Params, channel_coeffs

REPO = Path(__file__).resolve().parents[1]
P0 = dict(K=2.0, zG=0.08, zH4=0.80)


def gamma_numeric(P):
    """E(C3 along z) - E(S4 along z) = (2/3) Gamma for E = c - Gamma sum n^4."""
    e3, _ = minimize_cluster_at_orientation(P, oriented_tetrahedron([1, 1, 1], 0.3))
    e4, _ = minimize_cluster_at_orientation(P, oriented_tetrahedron([0, 0, 1], 0.3))
    return 1.5 * (e3 - e4)


def main():
    rows = []
    for a2, a4 in ((0, -0.02), (0, -0.01), (0, -0.005), (0.05, 0), (-0.05, 0)):
        P = Params.from_z(**P0, a2=a2, a4=a4)
        g, lamE, lamT = gamma_formula(P)
        rows.append(dict(a2=a2, a4=a4, formula=g, numeric=gamma_numeric(P)))
        print(rows[-1])
    rng = np.random.default_rng(20261119)
    P = Params.from_z(**P0, a4=-0.02)
    orient = []
    for nh in rng.normal(size=(40, 3)):
        nh /= np.linalg.norm(nh)
        e = minimize_cluster_at_orientation(P, oriented_tetrahedron(nh, rng.uniform(0, 2 * np.pi)))[0]
        orient.append(dict(K=float(np.sum(nh ** 4)), E=e))
    for nh in ([0, 0, 1], [1, 1, 1], [1, 1, 0]):
        nh = np.array(nh, float) / np.linalg.norm(nh)
        e = minimize_cluster_at_orientation(P, oriented_tetrahedron(nh, 0.3))[0]
        orient.append(dict(K=float(np.sum(nh ** 4)), E=e))
    kappa, _ = channel_coeffs(Params.from_z(**P0))
    out = dict(params=P0, lam_E=lamE, lam_T2=lamT, kappa=[float(k) for k in kappa],
               selection_combination=float(27 * kappa[0] - 63 * kappa[2] + 35 * kappa[3]),
               gamma_formula_a4=gamma_formula(P)[0], table=rows, orientations=orient)
    (REPO / "data" / "n4").mkdir(parents=True, exist_ok=True)
    (REPO / "data" / "n4" / "gamma_check.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
