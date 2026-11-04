# Spin_Body_Anisotropy — Paper II

Spin bodies (tetrahedral clusters of n = 4 Heisenberg spins) under uniaxial single-spin anisotropy.
Follow-up of *Deformable polyhedral spin bodies* (Paper I, repo `cywon-commits/Polyhedron_Spin_Group`).

Questions: which body axis the anisotropy selects (effective cubic coefficient Γ), and the 2D phase
transitions of the selected C₃ phase, whose effective model is the general XY–Ising model

    H = - sum_<ab> [ J1 s s' + (J2 + J4 s s') cos(phi - phi') ],   r = J4/J2 = (4/9) zH4/zG,

with an Ising (octupole sign, TR odd) and a KT (hexadecapole, TR even) line that cross at a
decoupled tetracritical point near r ≈ 4.3.

```
pip install -e .[test]
pytest -q                                   # ~20 s; includes the exact 6x6 Ising check
python -m sbaniso.pt --r 4.1 --L 32 --nrep 14 --Tmin 6.3 --Tmax 8.0 --ntherm 1000 --nmeas 10000 --seed 1
python src/sbaniso/analysis_pilot.py        # in a directory of pt_r*_L*.npz files
```

## Layout
```
src/sbaniso/c3model.py        effective model: energies, Metropolis, three Wolff cluster types, measurements
src/sbaniso/pt.py             parallel tempering driver (npz output; to be moved to YAML/HDF5, docs/TASK.md N0)
src/sbaniso/analysis_pilot.py Binder crossings, Upsilon = 2T/pi, Weber-Minnhagen fit, bimodality
tests/                        unit and statistical tests (must pass before any production run)
configs/T1_runs.csv           production run list for task N1
data/pilot_v3/                pilot data (r = 3.8, 4.1, 4.4; L = 32, 48, 64) and summary
docs/TASK.md                  work order for Claude Code (Korean)
docs/theory_summary.md        predictions P1-P10 and the formulas they come from
docs/pilot_report.md          pilot results and the two algorithm lessons
docs/note/                    derivation note (TeX + PDF)
reference/                    exploratory scripts used in chat (two-body relaxation, T=0 window scan)
```

## Conventions
Units G = 1 for the effective model (rigid couplings J1 = 160/81, J2 = 128/81). Full model: J = 1,
couplings via zG, zH4 as in polyspin; anisotropy V = sum_i [a2 P2(z.S_i) + a4 P4(z.S_i)].
Every production run records seed, parameters and git hash with its data.

License: MIT. Keep the repository private until submission.
