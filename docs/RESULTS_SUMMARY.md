# 결과 요약 (N0–N3f): 원고 입력용 (2026-11-19)

등급: [증명] 해석, [수치확인] 오차와 함께 확인, [추측]. 단위: 유효 모형은 G = 1, 전체 모형은 J = 1.
이론 번호(P1–P11)는 `docs/theory_summary.md`, 그림은 `figures/fig1–fig5` (아래 3절).

## 1. 핵심 수치

| # | 단계 | 결과 | 값 | 등급 | 이론 | 근거 | 그림 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 이론 | a₂의 1차 기여 없음 (정사면체 = 구면 2-디자인) | Σᵢ(ẑ·Sᵢ)² = 4/3 | [증명] | P1 | note 2절 | 1a |
| 2 | 이론·N3 | 축 선택 계수 Γ = (35/9)a₄ + 2a₂²(3/λ_E − 1/λ_T₂). Γ < 0이면 C₃ ∥ z | 단일 체 완화 대 공식: a₄ = −0.02에서 +0.6%, −0.01에서 +0.3%, −0.005에서 +0.15%. ±a₂ 평균 +0.06% | [증명] [수치확인] | P2 | N3_report 1절, `data/n4/gamma_check.json` | 2a, 2c |
| 3 | 이론 | a₂ 단독 선택 경계 3λ_T₂ − λ_E ∝ 27κ₁ − 63κ₃ + 35κ₄ = 27J + 9zG − 8zH₄ (κ₂ 상쇄). T 상(K ≤ 2)에서는 S₄ 선택이므로 C₃ 상에는 a₄ < 0이 필요 | (K = 2, zG = 0.08, zH₄ = 0.80): 21.32 > 0 | [증명] | P3 | note 4절 | 2b |
| 4 | 이론 | C₃ 유효 모형 (일반 XY–Ising). 강체 결합 J₁ = 160G/81, J₂ = 128G/81, J₄ = 512H₄/729 | 기계 정밀도로 수치 확인 | [증명] | 유효 모형 절 | note 5절 | 1c |
| 5 | N0 | 분석 사슬 검증: XY Weber–Minnhagen T_KT | 0.8906 (정확값 0.8929, 편향 −0.25%) | [수치확인] | — | N0_report | — |
| 6 | N0 | Ising 극한 T_c, U* | 2.2697, 0.6102 | [수치확인] | — | N0_report | — |
| 7 | N1 | Ising–KT 교차 (사전 등록 규칙 전부 통과) | r_× = 4.340 ± 0.001 (통계) ± 0.024 (계통), T_× = 6.995 ± 0.001 ± 0.007 G | [수치확인] | P6–P8 | N1_report, `data/n1/n1_summary.json` | 3a, 3b |
| 8 | N1 | ΔT = T_I − T_KT의 부호 변화 | r = 4.3: +0.033(2), r = 4.4: −0.021(2) | [수치확인] | P8 | N1_report | 3b |
| 9 | N1 | 분리 근사 예측 대 측정 | r_× 4.31 (계통 오차 안), T_× 7.07 (1% 높음) | [수치확인] (비교), 예측식은 근사 | P8 | N1_report | 3a |
| 10 | N1 | 분리된 사중임계점 (교차점 양쪽 r = 4.2, 4.4): Ising 지수와 KT가 따로 남음 | ν = 1.023(17), 1.004(15). U* = 0.61–0.62. C_max ∝ ln L (p = 0.43, 0.31), L² 기각. 히스토그램 단봉 | [수치확인] | P9 | N1_report | 4a–c |
| 11 | N1 | 결합 연산자 ε_σ(∂ψ)²는 무관 (차원 3 > 2) | — | [증명] (선형 RG) | P9 | theory_summary | 4 |
| 12 | N1 | 1차 합류 징후 없음 | 모든 r, L에서 히스토그램 단봉, C_max/L² 감소 | [수치확인] | P9 | N1_report | 4c |
| 13 | N2 | r = 6.0 뒤집힌 순서 구간 (기준 A–D 통과) | T_I = 7.322(2) < T < T_KT = 8.940(5) G | [수치확인] | — | N2_report | 3a |
| 14 | N2 | r = 5.0 구간 | T_I = 7.159(2), T_KT = 7.699(4) G. q = ¼에서 U₄ 비 0.526 (기준 < 0.5 미충족, ξ_σ ≈ 26 포화로 설명; 의뢰자 결정으로 L = 192 미실행) | 구간 존재 [수치확인], q = ¼ 기준 미충족 | — | N2_report 3–4절 | 3a |
| 15 | N2 | 구간 안 ψ는 준장거리, φ는 지수 감쇠 | η_ψ = 0.16–0.22, 유효 η_φ = 1.37–2.02 | [수치확인] | — | N2_report | — |
| 16 | N3 | polyspin 이방성 확장 | `v0.2.0-aniso`, Paper I 테스트 170개 통과 | [수치확인] | — | N3_report 1절 | — |
| 17 | N3 | 이완 결합의 r_eff (zH₄ = 0.72, 0.80, 0.88) | 3.677, 4.050, 4.418 | [수치확인] | 강체 근사 보정 | N3_report | — |
| 18 | N3 | 전체 모형 T_c가 같은 결합의 유효 모형보다 낮음 (a₄ = −0.02) | 약 40% 낮음 (z > 250). 순서 변화 미관측 | 40%는 [수치확인], 순서는 [추측] (평형 미달) | P10 | N3_report | 5a |
| 19 | N3s | 결합 축소 s로 ℓ_Γ를 줄이면 R = T_Υ^full/T_eff가 1로 다가감 | ℓ_Γ = 4.46, 2.44, 1.41, 0.77에서 R = 0.603, 0.692, 0.773, 0.866 | 경향 [추측] (사전 기준 미충족) | P10 | N3s_report | 5a |
| 20 | N3e | C₃ 축 기울기 직접 측정 (T_Υ^full) | ⟨θ²⟩ = 0.0973 (s = 1), 0.0261 (s = 0.03). 가우스 예측의 1.84, 1.28배 | [수치확인] (측정량) | P10 | N3e_report | — |
| 21 | N3e | 체별 DW 인자 1 − ½[l(l+1) − m²]⟨θ²⟩ | R = 1.094 (s = 1, 과보정), 0.988 (s = 0.03) | 사전 기준 미충족 | P10, P11 | N3e_report | 5a, 5c |
| 22 | 이론 | DW 계수 ½[l(l+1) − m²] (J₁: 6, J₂: 1.5, J₄: 5.5). 공통 기울기는 결합 불변. 정확한 결합 인자 ⟨d^l_mm(β_ab)⟩ | 단위 시험 (작은 β 전개, 독립 Wigner d, 공통 기울기 → 1) | [증명] | P11 | theory_summary P11, `tests/test_fullmodel.py` | 5c |
| 23 | N3f | 이웃 기울기 상관 | C_nn = 0.315 (s = 1), 0.104 (s = 0.03). ⟨β²⟩ ≈ 2⟨θ²⟩(1 − C_nn) | [수치확인] (측정량) | P11 | N3f_report | 5b |
| 24 | N3f | 상대 기울기의 정확한 DW 평균 | R = 0.801 (s = 1, 과소보정), 0.969 (s = 0.03) | 사전 기준 미충족 | P11 | N3f_report | 5a |
| 25 | N3 종합 | ℓ_Γ < 1 (κθ² ≲ 0.2)에서 DW 보정 유효 모형이 5% 안 | 세 보정(N3d, N3e, N3f)에서 0.96–0.99 | [추측] (전체 모형 평형 미달) | P10, P11 | N3f_report 3절 | 5a |
| 26 | N3 종합 | s = 1 (ℓ_Γ = 4.5)의 남는 차이는 결합 DW 감소 외 효과 (체 변형·결함, 2차 누적량) | 체별 1.094와 상대 0.801이 전체 모형을 감쌈 | [추측] | P10 | N3f_report 2절 | 5a |

**N3 원고 문구 (사전 고정, 실패 쪽):** "기울기 요동의 Debye–Waller 감소가 결손의 주원인이며, 선형 인자는 κθ² ≲ 0.2에서 5% 안에서 맞는다. 상대 기울기 평균으로도 s = 1의 차이가 남아, 결합 외 효과(체 변형·결함)를 시사한다."

## 2. 이론 항목과 수치 결과의 연결

| 이론 | 내용 | 등급 | 수치 결과 (위 표 번호) |
| --- | --- | --- | --- |
| P1 | 2-디자인: a₂의 1차 기여 0 | [증명] | 2 (a₂ 항이 2차임을 ±a₂ 대칭으로 확인) |
| P2 | Γ 공식, E_eff = −ΓΣn_α⁴ | [증명] | 2 (0.6% 이내, 차이는 다음 차수) |
| P3 | 선택 경계 27κ₁ − 63κ₃ + 35κ₄, κ₂ 상쇄 | [증명] | 3 |
| 유효 모형 | 강체 결합, ψ 틀, 등방 강성 | [증명] | 4, 7 |
| P6 | σ 질서 영역 T_KT ≈ 0.893(J₄ + J₂⟨σσ′⟩) | 구 코드 수치는 폐기, 식은 근사 | 7, 9 (해석선이 측정 T_KT보다 약 1.5% 높음, 그림 3a) |
| P7 | Ising 선 T_I ≈ 2.269(J₁ + J₂⟨cosΔψ⟩) | 근사 | 9 (교차점 근처에서만 맞음, 그림 3a) |
| P8 | 교차점 r_× = 4.31, T_× = 7.07 | 근사 (정확한 임계 상관 사용) | 7, 9 |
| P9 | 결합 연산자 무관 → 분리된 사중임계 교차 | [증명] (선형 RG) | 10, 12 |
| P10 | 유효 모형은 ⟨θ²⟩ → 0에서 정확, 선행 보정은 DW | [추측] | 18–21, 24–26 |
| P11 | DW 계수, 공통 기울기 불변, ⟨d^l_mm(β)⟩ | [증명] | 22–24 |
| 강체 근사 보정 | 두 체 완화로 r_eff ≈ 0.94 r_rigid | [수치확인] | 17 |
| 물리적 창 | T–Sq 경계 zH₄ = 0.692K + 0.923zG | [수치확인] | 3 (그림 2b) |

## 3. 그림과 캡션 (원고용, 영어)

모든 그림은 단일 단 3.4 in, Okabe–Ito 색맹 구분 팔레트, 공통 형식 `scripts/figstyle.py`. 각 스크립트는 저장소 데이터만 읽고 `figures/<name>.pdf`와 `.png`를 쓴다. 재현 명령은 README에 있다.

**Fig. 1** (`scripts/fig1_model.py`). A tetrahedral spin body under uniaxial anisotropy keeps one body axis along ẑ, and its remaining state is an Ising octupole sign σ times an XY hexadecapole phase ψ. (a) Spin body with the C₃ axis along ẑ (selected for Γ < 0) and one S₄ axis. (b) σ = sgn q₃₀; time reversal inverts the body, so σ is odd under 𝒯. (c) Top view: φ = 3α + π(1 − σ)/2 (q₃₃) and ψ = 3α (q₄₃); multipole weights and the couplings J₁, J₂, J₄ they produce.

**Fig. 2** (`scripts/fig2_axis_selection.py`, data `scripts/n4_gamma_check.py`). The anisotropy selects the body axis through one cubic coefficient Γ, which the second-order formula gives to within 0.6%. (a) Relaxed single-body energy at fixed orientation vs K = Σn_α⁴ (43 orientations, a₄ = −0.02): a line of slope −Γ. (b) Selection by a₂ alone, 27κ₁ − 63κ₃ + 35κ₄ = 27J + 9zG − 8zH₄ (κ₂ cancels); the boundary lies above the T–Sq boundary, so the C₃ phase needs a₄ < 0. (c) Deviation of the numerical Γ from the formula: linear in the amplitude (next order); the ±a₂ average agrees to 0.06%.

**Fig. 3** (`scripts/fig3_phase_diagram.py`). The Ising and KT lines of the C₃ effective model cross at r_× = 4.340(1)(24), T_× = 6.995(1)(7) G, beyond which the order of the two transitions is inverted. (a) T–r phase diagram: T_I (Binder crossings, L → ∞) and T_KT (Weber–Minnhagen) for r = 3.6–5.0 and r = 6.0 (open frames: inverted-order runs), the decoupling approximation (dashed), and the four phases. (b) T_I − T_KT vs r; open circles use T_KT × 1.0025 (Weber–Minnhagen bias bound); grey band: systematic uncertainty of r_×.

**Fig. 4** (`scripts/fig4_crossing_collapse.py`). On both sides of the crossing (r = 4.2, 4.4) the Ising and KT transitions keep their own scaling, as expected for a decoupled tetracritical point. (a) σ Binder cumulant vs (T − T_I)L (ν = 1): all sizes and both r on one curve. (b) ψ helicity modulus with the Weber–Minnhagen logarithmic factor vs (T/T_KT − 1)λ_L², λ_L = ln L + C/2 (C from the data at T_KT; the common point at 0 is by construction, the slopes converge with L). (c) C_max grows as ln L (fits p = 0.43, 0.31; L² rejected).

**Fig. 5** (`scripts/fig5_debye_waller.py`). The effective model with Debye–Waller-reduced couplings reproduces the full spin-body model within 5% for ℓ_Γ < 1, while at ℓ_Γ = 4.5 the per-body and relative-tilt prescriptions bracket it. (a) R = T_Υ^full/T_eff vs ℓ_Γ: bare couplings, per-body factor with the measured tilt, exact bond average ⟨d^l_mm(β_ab)⟩; grey band ±5%. (b) Distribution of the relative tilt β_ab at T_Υ^full with the Gaussian of variance 2⟨θ²⟩(1 − C_nn). (c) Coupling factors J^R/J per channel (l, |m|).

이전 초안 그림 (`figures/n1_*.png`, `n2_window.png`, `n3_full_vs_eff.png`, `n3s_ratio.png`, `n3d_ratio.png`)은 각 보고서가 참조하므로 그대로 둔다. 원고에는 fig1–fig5를 쓴다.

## 4. 남은 한계 (원고에 그대로 적을 것)
- 전체 모형은 국소 이동만 있어 평형 기준을 통과하지 못했다 (τ_int 수백–수천 스윕). 전체 모형 수치의 오차는 무작위·정렬 시작 폭의 절반이고, 관련 결론은 [추측] 등급이다.
- r = 5.0의 q = ¼ 기준 미충족 (N2), N3s·N3e·N3f 사전 기준 미충족은 그대로 보고한다.
- ΔT(r) 직선 맞춤의 χ²/자유도 = 526/8: r_×의 통계 오차는 모형 오차를 포함하지 않으며, 계통 오차 ±0.024가 이를 덮는다.
- 교차점이 Ising–KT 근접 영역에 있다는 것 자체가 선택 효과일 수 있다 (N3d_report의 단서: 이완 결합의 r^R − r_×^R = −0.32, −0.15, +0.01, +0.11).
