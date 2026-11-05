# Results ledger

| date | task | config | git hash | outcome |
| --- | --- | --- | --- | --- |
| 2026-11-04 | pilot v3 | r = 3.8, 4.1, 4.4; L = 32, 48, 64; PT, fixed cluster counts | (chat, pre-repo) | see docs/pilot_report.md |
| 2026-11-04 | N0 validation | configs/n0/ (24 runs: XY L = 16–64; Ising L = 16–64; r = 0 and J2 = 0 limits L = 16–64; pilot regression r = 3.8, 4.1, 4.4 L = 32 + ordered start) | 8f092f6-dirty, src_sha256 da562a92ef403822 | all N0 gates pass; WM T_KT = 0.8906(2) vs 0.8929 (−0.25 % fit bias); Binder T = 2.2697(3), U* = 0.6102(5); see docs/N0_report.md |
| 2026-11-05 | N0 review controls | runs/controls: r = 4.1, L = 64, 1000 + 10000; flags (1,0,1) and (1,1,1) | 8f092f6-dirty (pre-commit of bcfeded) | negative control rejected by equilibration_check (tau 921, 2 round trips); positive passes; see N0 report appendix |
| 2026-11-05 | N1 | configs/n1/ (46 runs from T1_runs.csv + ordered starts r = 4.2, 4.4, L = 64) | 8c317a3 | all 48 pass equilibration_check; all pre-registered order rules pass; r_x = 4.340 ± 0.001 (stat) ± 0.024 (syst), T_x = 6.995 ± 0.001 ± 0.007; decoupled tetracritical criteria pass at r = 4.2, 4.4; no first-order signs; see docs/N1_report.md |
