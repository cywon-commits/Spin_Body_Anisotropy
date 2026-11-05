# Results ledger

| date | task | config | git hash | outcome |
| --- | --- | --- | --- | --- |
| 2026-11-04 | pilot v3 | r = 3.8, 4.1, 4.4; L = 32, 48, 64; PT, fixed cluster counts | (chat, pre-repo) | see docs/pilot_report.md |
| 2026-11-04 | N0 validation | configs/n0/ (24 runs: XY L = 16–64; Ising L = 16–64; r = 0 and J2 = 0 limits L = 16–64; pilot regression r = 3.8, 4.1, 4.4 L = 32 + ordered start) | 8f092f6-dirty, src_sha256 da562a92ef403822 | all N0 gates pass; WM T_KT = 0.8906(2) vs 0.8929 (−0.25 % fit bias); Binder T = 2.2697(3), U* = 0.6102(5); see docs/N0_report.md |
