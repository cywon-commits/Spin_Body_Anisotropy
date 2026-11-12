# Results ledger

| date | task | config | git hash | outcome |
| --- | --- | --- | --- | --- |
| 2026-11-04 | pilot v3 | r = 3.8, 4.1, 4.4; L = 32, 48, 64; PT, fixed cluster counts | (chat, pre-repo) | see docs/pilot_report.md |
| 2026-11-04 | N0 validation | configs/n0/ (24 runs: XY L = 16–64; Ising L = 16–64; r = 0 and J2 = 0 limits L = 16–64; pilot regression r = 3.8, 4.1, 4.4 L = 32 + ordered start) | 8f092f6-dirty, src_sha256 da562a92ef403822 | all N0 gates pass; WM T_KT = 0.8906(2) vs 0.8929 (−0.25 % fit bias); Binder T = 2.2697(3), U* = 0.6102(5); see docs/N0_report.md |
| 2026-11-05 | N0 review controls | runs/controls: r = 4.1, L = 64, 1000 + 10000; flags (1,0,1) and (1,1,1) | 8f092f6-dirty (pre-commit of bcfeded) | negative control rejected by equilibration_check (tau 921, 2 round trips); positive passes; see N0 report appendix |
| 2026-11-05 | N1 | configs/n1/ (46 runs from T1_runs.csv + ordered starts r = 4.2, 4.4, L = 64) | 8c317a3 | all 48 pass equilibration_check; all pre-registered order rules pass; r_x = 4.340 ± 0.001 (stat) ± 0.024 (syst), T_x = 6.995 ± 0.001 ± 0.007; decoupled tetracritical criteria pass at r = 4.2, 4.4; no first-order signs; see docs/N1_report.md |
| 2026-11-06 | N2 pilots | configs/n2/pilot/ (r = 5.0, 6.0; L = 32; 1000 + 5000) | 6d39759 | used only for C(T) and temperature ranges |
| 2026-11-06 | N2 | configs/n2/ (r = 5.0, 6.0; L = 32–96; designed ladders; + ordered start r = 6.0, L = 64) | 115dcdc | all 9 pass equilibration_check; r = 6.0 window confirmed (T_I 7.322, T_KT 8.940); r = 5.0 window exists (0.541 ± 0.005) but one pre-registered sub-criterion fails (U4 ratio 0.526 at q = 1/4; finite-size, supplementary xi_sigma ≈ 26); eta_psi 0.16–0.22, eta_phi 1.4–2.0; see docs/N2_report.md |
| 2026-11-06 | N3 polyspin extension | Polyhedron_Spin_Group branch anisotropy-extension | edda213 (tag v0.2.0-aniso) | 170 Paper I tests pass (157 old incl. slow + 13 new), mutation checks caught; Gamma formula 0.6 % at a4 = -0.02 |
| 2026-11-06 | N3 effective model | configs/n3/eff (relaxed couplings of zH4 = 0.72, 0.80, 0.88; L = 16-64) | 17aaba5 | all 12 pass equilibration_check; T_I = 0.1365, 0.1395, 0.1417; T_KT = 0.1302, 0.1386, 0.1463 J; crossing near r_eff 4.11 |
| 2026-11-06 | N3 full model | configs/n3/full (3 points x L = 16, 24, 32 + ordered start) | 17aaba5, polyspin v0.2.0-aniso | none pass equilibration_check (tau 300-1800, 2-10 round trips); single joint transition at T_c ~ 0.080-0.088 J, ~40 % below the effective model (robust: random/ordered starts bracket it within 0.0015); TASK rule: disagreement; see docs/N3_report.md |
| 2026-11-06 | N3 pilots | configs/n3/pilot | 17aaba5 era | first pilots discarded (stale numba cache, no anisotropy); valid pilots used only for windows |
