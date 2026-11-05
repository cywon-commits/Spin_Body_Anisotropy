"""Write the N1 YAML configs from configs/T1_runs.csv to configs/n1/, plus the ordered-start
controls of docs/N1_analysis_plan.md (r = 4.2, 4.4; L = 64), and the run order file
configs/n1/ORDER (priority, then L ascending, then r).

    python scripts/make_n1_configs.py
"""
import csv
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "configs" / "n1"


def main():
    rows = [r for r in csv.DictReader(open(REPO / "configs" / "T1_runs.csv")) if r["task"] == "N1"]
    OUT.mkdir(parents=True, exist_ok=True)
    order = []
    for i, r in enumerate(rows):
        name = f"r{r['r']}_L{r['L']}"
        d = dict(name=f"n1_{name}", model=dict(r=float(r["r"])), L=int(r["L"]),
                 temperatures=dict(Tmin=float(r["Tmin"]), Tmax=float(r["Tmax"]), n=int(r["nrep"])),
                 ntherm=int(r["ntherm"]), nmeas=int(r["nmeas"]), clusters=[1, 1, 1],
                 seed=10000 + i, start="random", checkpoint_every=1000, output=f"runs/n1/{name}.h5")
        (OUT / f"{name}.yaml").write_text(yaml.safe_dump(d, sort_keys=False))
        order.append((int(r["priority"]), int(r["L"]), float(r["r"]), name))
    for j, rr in enumerate((4.2, 4.4)):
        src = yaml.safe_load((OUT / f"r{rr}_L64.yaml").read_text())
        name = f"r{rr}_L64_ordered"
        src.update(name=f"n1_{name}", seed=20000 + j, start="ordered", output=f"runs/n1/{name}.h5")
        (OUT / f"{name}.yaml").write_text(yaml.safe_dump(src, sort_keys=False))
        order.append((1, 64, rr + 0.01, name))
    order.sort()
    (OUT / "ORDER").write_text("".join(f"configs/n1/{n}.yaml\n" for *_, n in order))
    print(f"{len(order)} configs; first: {order[0][-1]}, last: {order[-1][-1]}")


if __name__ == "__main__":
    main()
