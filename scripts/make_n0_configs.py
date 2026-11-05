"""Write the YAML configs of the N0 validation runs (docs/TASK.md, N0) to configs/n0/.

    python scripts/make_n0_configs.py

Families (units: couplings c = [J1, J2, J4] as in sbaniso.c3model)
  xy_L*      pure XY (J1 = J4 = 0, J2 = 1): Ups = 2T/pi crossings, Weber-Minnhagen T_KT = 0.8929
  ising_L*   pure Ising (J2 = J4 = 0, J1 = 1): Binder crossing T = 2.269, U* = 0.611
  r0_*       r = 0 (J4 = 0; sigma and phi decouple), rigid J1, J2:  T_I = 2.269 J1, T_KT = 0.8929 J2
  j2zero_*   J2 = 0 (sigma and psi decouple), J4 = 4 * (128/81):     T_I = 2.269 J1, T_KT = 0.8929 J4
  pilot_*    regression against data/pilot_v3 (L = 32, same ladders), random and ordered starts
"""
from pathlib import Path

import yaml

J1 = 160 / 81; J2 = 128 / 81; J4 = 4 * 128 / 81
OUT = Path(__file__).resolve().parents[1] / "configs" / "n0"


def cfg(name, couplings, L, Tmin, Tmax, n, clusters, ntherm=2000, nmeas=40000, seed=0, **kw):
    d = dict(name=f"n0_{name}", model=dict(couplings=[float(x) for x in couplings]), L=L,
             temperatures=dict(Tmin=float(Tmin), Tmax=float(Tmax), n=n), ntherm=ntherm, nmeas=nmeas,
             clusters=clusters, seed=int(seed),
             start="random", checkpoint_every=2000, output=f"runs/n0/{name}.h5")
    d.update(kw)
    return name, d


def all_configs():
    out = []
    reps = {16: 8, 24: 10, 32: 12, 48: 16, 64: 20}
    for i, L in enumerate((16, 24, 32, 48, 64)):
        out.append(cfg(f"xy_L{L}", [0, 1, 0], L, 0.82, 1.02, reps[L], [0, 0, 1], seed=1000 + i))
    for i, L in enumerate((16, 32, 64)):
        out.append(cfg(f"ising_L{L}", [1, 0, 0], L, 2.18, 2.36, {16: 6, 32: 8, 64: 12}[L],
                       [1, 0, 0], seed=2000 + i))
    for i, L in enumerate((16, 32, 64)):
        # decoupled limits: total-energy fluctuations include both subsystems -> more replicas
        k = {16: 10, 32: 14, 64: 22}[L]
        out.append(cfg(f"r0_kt_L{L}", [J1, J2, 0], L, 0.82 * J2, 1.02 * J2, k, [1, 1, 1], seed=3000 + i))
        out.append(cfg(f"r0_ising_L{L}", [J1, J2, 0], L, 2.18 * J1, 2.36 * J1, k, [1, 1, 1], seed=3100 + i))
        out.append(cfg(f"j2zero_kt_L{L}", [J1, 0, J4], L, 0.82 * J4, 1.02 * J4, k, [1, 1, 1], seed=3200 + i))
        out.append(cfg(f"j2zero_ising_L{L}", [J1, 0, J4], L, 2.18 * J1, 2.36 * J1, k, [1, 1, 1], seed=3300 + i))
    for i, r in enumerate((3.8, 4.1, 4.4)):
        out.append(cfg(f"pilot_r{r}_L32", [J1, J2, r * J2], 32, 6.3, 8.0, 14, [1, 1, 1],
                       ntherm=1000, nmeas=20000, seed=4000 + i))
    out.append(cfg("pilot_r4.1_L32_ordered", [J1, J2, 4.1 * J2], 32, 6.3, 8.0, 14, [1, 1, 1],
                   ntherm=1000, nmeas=20000, seed=4100, start="ordered"))
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for name, d in all_configs():
        (OUT / f"{name}.yaml").write_text(yaml.safe_dump(d, sort_keys=False))
    print(f"{len(all_configs())} configs in {OUT}")
