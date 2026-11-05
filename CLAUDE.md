# Working rules for this repository

- Read `docs/TASK.md` first; it is the work order. Background in `docs/theory_summary.md`, lessons in
  `docs/pilot_report.md`.
- Reports to the owner are in Korean (as in Paper I's docs/M*_report.md); code, comments, figure
  labels and the manuscript in English.
- Order of work: symmetry and analytic checks, then small exact checks, then large runs. A large run
  starts only after `pytest` passes, including the slow statistical tests.
- Monte Carlo rules
  1. Never stop cluster building on a size criterion within a sweep. The number of clusters per sweep
     is tuned during the first half of equilibration and then frozen.
  2. Always include the psi-frame sigma cluster (`wolff_sigma_psi`) for the C3 effective model.
  3. Equilibration: quarter-split averages agree within errors, random and ordered starts agree,
     >= 20 full PT round trips (bottom -> top -> bottom), swap acceptance >= 0.2 per attempt for every
     neighbour pair. Report tau_int for every temperature used in a fit.
- Decision rules in `docs/TASK.md` are fixed before looking at production data. If a result needs a
  rule change, write the reason in the report before rerunning the analysis.
- Grade every claim in reports: [proved] analytic, [numerical] checked with errors, [conjecture].
- Record each completed run in `docs/RESULTS_LEDGER.md` (date, config, git hash, outcome).
- Do not modify Paper I's repository except on a dedicated branch for the anisotropy extension
  (task N3); all Paper I tests must still pass there.
