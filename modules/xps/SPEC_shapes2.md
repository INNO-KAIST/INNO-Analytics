# CasaXPS Cookbook 2026 선모양 (v1.9 · SPEC_shapes2)

정의 출처: *CasaXPS Cookbook 2026* (Casa Software Ltd), pp. 65–68 (식), p. 19 (LA(m) ≡ LA(1,701) ≡ LA(1,1,701)), 쓰임 예 p. 46 (ST(1)LA(1,143), TP(1,1,0)), p. 71 (SB(0)LA(1,43), TP(1,15,0)), p. 91 (DS(0.3,600,40)).
식만 옮겼고 CasaXPS 코드·자료는 쓰지 않았다. 엔진 `src/xps_core.js` §1b·§16, 기준 구현 `test/ref/ref_shapes2.py`, 시험 `test/test_shapes2.js`, `test/ui_xps_shapes2.py`.

## 0. 좌표·위치·폭
- Cookbook 식의 변수 x_c 는 FWHM 단위이고 **x_c = −x/F** (INNO x = E − E0, x > 0 = 높은 BE)로 옮긴다. 근거: LA 의 α(x_c ≤ 0), DS 의 긴 꼬리, E(j,k)(x_c < 0), ∫_{x_c}^∞ 가
  모두 x_c < 0 쪽이고, CasaXPS LA(1.53,2.43) 실측(SPEC_casa_v2 §2)이 α = 높은 BE 쪽이었다. 따라서 **꼬리·Shirley 계단은 모두 높은 BE 쪽**.
- 위치(BE) = 정의식의 x = 0 (CasaXPS 위치 칸과 같은 뜻으로 가정). 비대칭 모양의 최댓값은 x = 0 이 아니다 — 결과의 `beMax` 로 따로 보인다.
- F = 정의식의 단위 폭(DS 는 Lorentzian FWHM wL). 표의 FWHM 은 **실제 FWHM**(수치): 거친 표본 2401 점 → 최댓값·반높이 교차를 그 칸 안 65 점으로 다시 표본(선형 보간 오차 ~1e-9·폭). SB 는 FWHM 없음.
- **CasaXPS 정수 n → 폭**: LA 에서 실측한 `FWHM_G = n·F/144.87` (SPEC_casa_v2 §2.3, α = β = 1 에서만 측정)을 **다른 모양에도 가정**해 쓴다
  (TLA·GL(m,n)·SGL(m,n): Gaussian; DS: n·wL/144.87; TP: Lorentzian FWHM n·F/144.87; ST γ: n·F_바탕/144.87). `compToCasa` 는 모두 `exact:false` 와 "assumed (LA rule)" 로 표시.
  성분에 `casaN`(·`stN`)이 있으면 그 폭은 자유 매개변수에서 빠지고 F 와 함께 움직인다(화면 칸 잠김, 버튼 'n').

## 1. 모양 (넓이 1 = ℝ 적분, 단 §2 의 DS·SB 는 창 넓이)
| INNO | Cookbook | 정의 (x_c = −x/F, l(t) = 1/(1+4t²)) | 정규화 |
|---|---|---|---|
| `TLA` F, α, μ, wG | TLA(x:α,μ,w) p. 66 | T = l^α·[π/2 − atan(2x_c) + π/μ] ⊗ G(wG) | 닫힌 꼴 Z = 2F(π/2 + π/μ)·I₁, I₁ = (√π/4)Γ(α−½)/Γ(α). α > ½ 이어야 적분됨 → α ∈ [0.55, 20] |
| `GLN` F, m, wG [+꼬리] | GL(m,n) p. 67 | GL = exp(−(1−m)4ln2·x_c²)/(1+4m x_c²) ⊗ G(n) | ∫GL = glpNorm (§1.3) |
| `SGLN` F, m, wG [+꼬리] | SGL(m,n) p. 67 | SGL = m·l(x_c) + (1−m)exp(−4ln2 x_c²) ⊗ G(n) — **m 은 높이 비** (GLS 의 넓이 비 아님) | (1−m)√(π/4ln2) + mπ/2 |
| 꼬리 `tail:"T"` tk | GL(m,n)T(k) p. 67 | h = B + (1−B)·E(1,k), E(j,k) = j·e^{k x_c} (x_c < 0) | ∫h = ∫B + j/κ − j∫₀^∞B e^{−κu}du (Gauss–Kronrod, κ = k) |
| 꼬리 `tail:"PHI"` tj, tk | ·PHI(j,k) p. 67 | h = B + (1−B)·E(j, 13.8/k) | 같은 식, κ = 13.8/k |
| `TP` F, ta(α), tb(β), wL | TP(α,β,n) p. 68 | R = α + (β−α)(x_c+½), \|x_c\| ≤ ½ (높은 BE 끝 = α) ⊗ L(wL) | (α+β)F/2; ⊗L 은 atan·log 닫힌 꼴. α/β 비만 정해지므로 **ta 는 기본 고정** |
| `DS` α, wL, wG, [dsl] | DS(α,n), p. 66 | P = cos(πα/2 + (1−α)atan(2x_c))/(1+4x_c²)^{(1−α)/2} ⊗ G | §2 창 넓이 (기존 §1.5 와 같은 함수) |
| DS + `dsl` = l | DS(α,n,l) p. 66–67 | Q = P − k·S[P] (x_c ≥ −l), 0 (x_c < −l); S[P] = ∫_{x_c}^∞ P = `dsPrim`; k = P(−l)/S[P](−l) (Q(−l) = 0) | §2 창 넓이 |
| 수식어 `mod:"ST"` stMu, stG | ST(μ,γ)F p. 66 | S(μ)F = N·F·[Φ + A/μ], Φ(x) = ∫_{x_c}^∞ F = INNO 의 ∫_{−∞}^{x} F; ⊗ G(γ) | ∫F·Φ = A²/2 → **N = 1/(A²(½ + 1/μ)) 닫힌 꼴** |
| 수식어 `mod:"SB"` (sbArg 0) | SB(0)F p. 68 | Φ(x) (0 → 1, 높은 BE 쪽으로 올라가는 계단) | §2 창 넓이 |

- 수식어 바탕: GLS, GLP, LA, V, TLA, GLN, SGLN (넓이 유한한 것). DS·TP·AGT·GLT·MEAS 에는 무시. ST 의 γ 는 0 이면 기본 고정(= ST(μ), p. 46 의 쓰임).
- 하위 호환: 저장된 `AGT`(tu "fwhm") 성분은 그대로 계산된다. CasaXPS 표기 `GL(p)T(k)` 를 읽으면 이제 `GLN`(tail T, n 0) — AGT 와 같은 모양(창 넓이로 맞춘 뒤 차 1e-15), 정규화만 ℝ.
  `GL(p)`, `SGL(p)`, `LA(…)` 는 이전처럼 GLP·GLS·LA.

## 2. 창 넓이 (DS · DS-Shirley · SB)
α > 0 의 DS 는 높은 BE 쪽 ~|x|^{α−1} 로 ℝ 적분이 발산하고, Q 는 높은 BE 쪽은 l 에서 잘리지만 낮은 BE 쪽 −k·S[P] ~ −k|x|^{α−1} 이 발산(α ≥ 0), SB 는 계단이라 발산한다.
→ **넓이 A = 피팅 창 [Emin, Emax] 안의 적분** (SPEC §1.5 의 DS·§1.10 AGT 와 같은 규칙): f = p / ∫_{Emin}^{Emax} p. 결과의 Area 는 '창 안 넓이'이고, 정량 탭에 경고가 뜬다(창 폭이 같은 시료끼리만 비교).
CasaXPS 가 DS 넓이를 어떤 구간에서 재는지는 Cookbook 에 없어 **미확인**이다(CasaXPS 대조 대기). l 의 단위는 x_c (= wL), 자름 위치 x = +l·wL (높은 BE).
SB 는 Cookbook 의 '0' 인자 뜻이 적혀 있지 않아 값만 보존(`sbArg`)하고 계산에 쓰지 않는다.

## 3. 수치 (모두 `s2Eval`)
- 합성곱·누적이 없으면 해석 값(DS 창 넓이는 원시함수, DS-Shirley 는 적응 Gauss–Kronrod 7–15, 구간 경계 0·l).
- 격자: 마디 x = k·h (성분 중심이 마디), h = 0.04·2^{−k} ≤ min(s0/16 (꺾임 있으면 /32), σ/2) 로 양자화(폭이 조금 바뀌어도 마디 고정 → 수치 Jacobian 매끄러움).
  s0 = 가장 좁은 구조의 폭(F, wL, Voigt FWHM, 꼬리 길이 4F/κ). 여유 8(σ₁+σ₂). σ < 0.0015·s0 이면 합성곱 생략(상대 영향 ≲ 1e-5).
- Gaussian 합성곱 = 해석 무게의 사다리꼴 합(±8σ; 매끄러운 피적분에 지수 수렴). σ_eff ≥ s0/2 이면 출력은 r 칸마다(r·h ≤ σ_eff/8)만 계산.
- Φ = ∫_{−∞}^x F: 바탕의 Φ(0)(해석: ½, LA 의 β 쪽 넓이, TLA 는 [(π/μ)I₁ + ½J]/(Z/F), J(α) = ∫₀^{π/2} φ sin^{2α−2}φ dφ (GK, α 별 캐시), GL·SGL 은 ½∫B/Z)에서
  칸마다 Gauss–Legendre 4 점(해석 값) 누적; 바탕에 Gaussian 이 있으면 Φ ⊗ G (적분과 합성곱은 교환). GLS 는 erf·atan 닫힌 꼴.
- 데이터 점 값 = 6 점 Lagrange 보간; ST(γ 0)·SB(Gaussian 없음)은 바탕 값을 해석으로, Φ 만 보간(꺾인 바탕도 정확). 창 적분 = 6 점 칸 무게 (11, −93, 802, 802, −93, 11)/1440 + 끝 조각 GL4.
- 오차 예산(아래 §5 로 확인): 값 ≤ 2e-7(최댓값 대비), 넓이 ≤ 1e-8, FWHM ≤ 5e-8 — 요구 1e-5 · 1e-6 · 1e-5 보다 작다.
- 합성곱 핵 ±7σ(잘린 무게 < 3e-12), 대칭 접기, 출력은 필요한 범위(데이터·창 + 다음 단계 여유)만 계산. `compEval` 은 비싼(격자) 모양의 넓이 1 값을 성분마다 최근 4 개 캐시
  (키 = 모양 매개변수·위치·창·E 배열) → LM 의 넓이 열·넓이 초기화는 모양을 다시 계산하지 않음(값은 비트 단위로 같음).
- 속도(node, 2 코어, 2000 채널, 성분 1 개): TP 0.2 ms, TLA·GLN 0.5–0.6, DS 1.2–2.3, ST/SB 0.8–2.5 ms. 10 성분(새 모양 모두, 자유 = BE·넓이·주 폭) 피팅 1.7–2.1 s, Monte Carlo 1.2–1.5 s/회.
  주의: 모양 매개변수까지 모두 풀면(예: TP 의 wL) 퇴화 골짜기로 가 wL → 하한(가장자리가 채널보다 날카로움)에서 LM 이 수백 번 기어가며(500 반복, ~65 s) — 계산 비용 문제가 아니라 문제 설정 문제.
  CasaXPS 처럼 모양 문자열 매개변수는 고정하는 것을 권한다(`casaModelToINNO` 는 그렇게 들여옴).
- 화면: 피팅은 Web Worker 에서 돌고, LM 반복마다(≤ 150 ms 간격) 진행 줄에 '반복 k · χ²' 를 보인다.

## 4. CasaXPS 표기 (`casaToComp` / `compToCasa`)
문법: `[ST(μ[,γ]) | SB(a)]` + 바탕 + `[T(k) | PHI(j,k)]`. 바탕 `GL(m)`·`GL(m,n)`·`SGL(m)`·`SGL(m,n)`·`LA(m)`·`LA(α,n)`·`LA(α,β,n)`·`TLA(α,μ[,n])`·`DS(α,n[,l])`·`TP(α,β[,n])`.
예: `ST(1.7)LA(1,200)`, `ST(1)LA(1,143)`, `DS(0.1,50)`, `DS(0.3,600,40)`, `GL(30,50)`, `SGL(30,50)T(2)`, `GL(30)PHI(1,2)`, `TP(1,15,0)`, `SB(0)LA(1,43)` — 모두 왕복(같은 성분, 곡선 차 0).
- `TLA(α,μ,n)` 표기는 Cookbook 에 문자열로 나오지 않는다(식은 TLA(x:α,μ,w)) — LA(α,β,n) 에 맞춰 **가정**.
- 폭: CasaXPS fwhm 칸을 실제 FWHM 으로 보고 F 를 정한다(대칭 GL·SGL·LA 에서 실측으로 맞음). 비대칭(TLA·꼬리·DS·ST)은 미검증 → "refit the width" 표시.
  꼬리 있는 GL/SGL 은 F = fwhm 숫자(실측 GL(30)T(2): CasaXPS fwhm = 0.984 F, SPEC_casa_v2 §2.2). TP 는 F = fwhm(가정).
- 모르는·잘못된 표기: `casaToComp(s,{strict:true})` 는 분명한 문구로 throw(화면 'Casa' 칸이 이것을 씀), 기본 호출은 이전처럼 `null`.
- `casaModelToINNO`: 새 모양·수식어 성분은 CasaXPS 처럼 표기 안의 모양 매개변수(α, μ, m, k, j, ta, tb, μ_ST …)를 고정으로 들여온다(LA·GL·SGL 은 이전 동작).

## 5. 검증 (`node test/test_shapes2.js`, 기준 `python3 -I test/ref/ref_shapes2.py`)
기준 구현은 Cookbook 식을 x_c 좌표 그대로 쓰고(INNO 의 닫힌 꼴 상수는 쓰지 않음) 정규화는 quad 로 ℝ·창 적분, 합성곱은 quad(점마다) 또는 FFT(h = 1/1024 F) + 누적 Simpson(ST·SB·DS-Shirley),
꼬리 질량은 quad 로 정확히 더함. TLA 상수는 mpmath(30 자리, Beta 함수 + quad)로 따로.

| 항목 | 사례 | 최대 오차 | 요구 |
|---|---|---|---|
| 값 (최댓값 대비) | 30 (TLA 4, GLN 5, SGLN 3, ST 7, SB 3, DS 5, TP 3) | 1.5e-7 (DS-Shirley, l = 10) | ≤ 1e-5 |
| 넓이 (ℝ 사례: ∫_{−R}^{R} vs 기준; 창 사례: 창 적분 − 1) | 30 | 7.1e-9 | ≤ 1e-6 |
| 실제 FWHM (상대) | 26 | 4.6e-8 | ≤ 1e-5 |
| TLA Z·Φ(0) vs mpmath | α 0.55–2.5 | 1.8e-15 | — |
| ST 정규화 항등식 (기준의 수치 ∫S vs A²(½+1/μ)) | 7 | 5.1e-14 | — |
| n → 폭 규칙 | 9 | 0 (정의) | — |
| 합성 Poisson 회복 (200 회, 자유 매개변수 모두 \|p̂−p\| ≤ 3σ) | 9 모양 | §5.1 | ≥ 90 % |

### 5.1 회복 (200 회, 배경 c₀ 자유)
| 모양 (자유 매개변수) | 3σ 적중률 |
|---|---|
| TLA α 1.4 μ 2 wG 0.5 (be, A, F, μ, wG) | 99 · 100 · 100 · 100 · 99 % |
| GL(30,·)T(2) wG 0.4 (be, A, F, m, wG, k) | 모두 100 % |
| SGL(30,·)PHI(1,2) wG 0.3 (be, A, F, wG, j) | 모두 100 % |
| ST(1.5)LA(1,1,·) wG 0.5 (be, A, F, wG, μ) | 모두 100 % |
| ST(1, γ 0.3)GLS (be, A, F, μ, γ) | 모두 100 % |
| DS(0.15) wL 0.3 wG 0.5 (be, A, α, wL, wG) | 모두 100 % |
| DS(0.2,·,l 20) (be, A, α, wL, wG) | 100 … 99 % |
| LA + SB(0)LA 계단 (be, A, F, A_SB) | 모두 100 % |
| GLS + TP(1,3) ⊗ L 0.5 (be, A, F, A_TP, F_TP, β) | 모두 100 % |

(Shirley k 를 0 으로 고정하려면 `bg._auto:false` 가 필요하다 — 없으면 넓이 초기화가 고정된 k 도 자동값으로 바꾼다. 기존 동작, 시험에서만 영향.)

UI (`test/ui_xps_shapes2.py`, Chromium): 모양 목록·꼬리·수식어·n·Cut l·Casa 입력·오류 문구·칸 편집·피팅(독립 Python 으로 만든 TLA 자료에서 α 1.298/1.3, μ 2.92/3, F 0.896/0.9, 넓이 +0.05 %;
TP(1,1,0) 배경 성분(폭·위치 고정)이 상수 400 CPS 를 0.15 % 로 회복)·390 px 넘침 없음.

## 6. 미확인 (CasaXPS 대조가 필요한 것)
1. n → 폭 규칙을 LA(1,1,n) 밖(α ≠ 1, TLA, GL(m,n), SGL(m,n), DS, TP 의 Lorentzian, ST 의 γ)에 쓴 것 — 가정.
2. DS·DS(α,n,l)·SB 의 넓이를 CasaXPS 가 어느 구간에서 재는지 — 여기서는 피팅 창.
3. 비대칭 모양에서 CasaXPS 'fwhm' 칸의 정의(SPEC_casa_v2 §2.2 의 LA·GL(30)T(2) 와 같은 문제), TP 의 fwhm = 폭 가정.
4. 위치 칸 = 정의식의 x = 0 이라는 가정(TLA·ST 는 최댓값과 다름), SB(a) 의 인자 a 의 뜻, TLA 문자열 표기.
