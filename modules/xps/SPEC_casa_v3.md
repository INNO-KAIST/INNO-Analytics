# CasaXPS 규약 3차 판독(Round 3) — Cookbook 2026 선모양 정의 대조 (2026-10-29, 분석만)

자료: CasaXPS 2.3.27PR7-6 demo 판독 `test/ref/casa_readings_v3.md` 끝 표 (Round 3) + V2 행 일부.
합성 자료 `test/data/casaval/K_*.vms` (생성 `test/ref/make_casaval3.js`: INNO `casaToComp(str,{be:285,area:15000,fwhm})` → `compEval`, 평탄 1000).
기준 계산 `test/ref/ref_casa_v3.py` (+ `casa3_lib.py`, numpy/scipy, Cookbook pp. 65–67 식만; CasaXPS 코드 없음) → `ref_casa_v3.json`.
CasaXPS 피팅 흉내: 영역 280.03/297.13, Linear(Av.W 1), 채널 i0+1…i1−1, 가중 없는 최소제곱, RMS = √Σr² (SPEC_casa_v2 §1). 자유 = 위치·넓이·폭.
코드 변경 제안: `test/ref/casa_v3_patch.md` — **적용 상태는 §7** (2026-10-29, 실측으로 뒷받침되는 항목만 xps_core.js 에 넣음).

## 1. 결론 요약
| CasaXPS | INNO 정의와 | 넓이 (INNO/Casa−1) | 위치 | RMS (INNO / Casa) | 판정 |
|---|---|---|---|---|---|
| GL(30,100) | 같음 | +3.3e-5 | 같음 | 0.00 / 2.78 | **정함** (n → Gaussian FWHM n·F/144.87, F = GL 의 F) |
| SGL(30,100) | 같음 | +5.7e-4 | 같음 | 362.1 / 359.7 | **정함** (높이 혼합 ⊗ G, RMS 는 Lorentzian 꼬리 자르기) |
| TLA(1.2,3,0) | 같음 | +1.9e-4 | 같음 | 292.1 / 290.6 | **정함** (꼬리 높은 BE, +π/μ, 표기 TLA(α,μ,n) 받아들여짐) |
| TLA(1.2,3,100) | 같음 | +4.3e-4 | 1e-4 | 185.6 / 180.1 (3 %) | **정함** (n 규칙 같음; n 스캔 최소는 n≈104, LA(1,1,n) 스캔과 같은 수준의 차이) |
| DS(0.1,100) | 모양 같음 | 정규화 **다름** | **+wL·\|cot(π/(2−α))\|** | 1573.7 / 1570.2 | 모양·n·꼬리 쪽 정함, 위치 이동 정함(1 판독), 넓이 미정 |
| ST(1.7)LA(1,1,200) | 모양 같음 | +7.2e-4 | **+0.0984 eV 차** | 325.5 / 322.3 | 모양 정함(μ 1.4–3.0 · n 스캔에서 μ 1.7·n 200 이 최적), 위치·폭 숫자 미정 |
| LA(2,3,0) | 같음 | +1.4e-4 | 같음 | 5.6 / 22.9 | 모양 같음, 폭 숫자는 §3 |
| LA(1,1,250) | 같음 | +8.7e-4 | 같음 | 426.3 / 417.3 (2.2 %) | **같음** — RMS 는 Linear 배경의 Lorentzian 꼬리 자르기(V2 n-스캔과 같은 차이) |
넓이 차 +2e-4…+9e-4 (Lorentzian 성분이 클수록 큼)는 SPEC_casa_v2 §4 의 CasaXPS 쪽 수치 차이와 같은 부호·크기.

## 2. DoF 규칙 — 정함
**DoF = #{채널 i0+1…i1−1 : y_i ≠ B_i} − (자유 매개변수 수)** (B = 영역 배경). 배경과 정확히 같은 채널은 세지 않는다.
GL(30,100): 154 − 3 = **151**, GL(30) (Shape_GLP30): 159 − 3 = **156**, GL(30)T(2): 340 − 3 = **337** (Gaussian 꼬리가 6 자리 반올림으로 정확히 1000.000000 이 되는 채널·
배경 끝점 채널이 빠짐). V2 §4 의 'T 모양이 매개변수 하나를 더 센다' 가설은 틀림. 음수 잔차 채널은 센다(LA(1,1,250) 은 78 채널이 y < B 인데 338).
RMS 는 그 채널들의 잔차가 0 이라 영향 없음. 정확한 모양에서도 남는 CasaXPS RMS 바닥(GL(30) 8.3, GL(30,100) 2.8, GL(30)T(2) 6.8)은 CasaXPS 의 모양 계산 오차로 둔다.

## 3. CasaXPS 'fwhm' 숫자와 비대칭 지수 (asymmetric shapes)
정의(가설 Hw): x̂ = 봉우리 최댓값 위치, a < x̂ < b = 반높이 점(a 낮은 BE = 가파른 쪽), L = x̂ − a, H = b − x̂ 일 때
**fwhm_Casa = 2L**, **asym = 1 − L/H**. 대칭이면 실제 FWHM, asym 0 과 같다.
| 모양 | 2L / Casa | 1 − L/H / Casa | 함의된 x̂ 오프셋 δ (F 단위, W·asym 두 값에서 각각) |
|---|---|---|---|
| TLA(1.2,3,0) | 0.911758 / 0.911717 (−4.5e-5) | 0.1529 / 0.1524 | ≈ 0 |
| TLA(1.2,3,100) | 0.954677 / 0.955100 (+4.4e-4) | 0.0822 / 0.0817 | ≈ 0 |
| DS(0.1,100) | 0.907160 / 0.907537 (+4.2e-4) | 0.1180 / 0.1178 | ≈ 0 |
| LA(2,3,0) | 0.884163 / 0.892084 (+0.9 %) | 0.2075 / 0.1936 | 0.0023 → asym 0.1947 |
| LA(1.53,2.43,0) F1·F2 | 0.57463 F / 0.5800 F (+0.9 %) | 0.2408 / 0.2275 | 0.0027 · 0.0026 → 0.2280 · 0.2286 |
| GL(30)T(2) | 0.96211 F / 0.98376 F (+2.3 %) | 0.2445 / 0.2138 | 0.0108 (W) · 0.0109 (asym) |
| GL(30)T(1) · GL(0)T(2) (asym 만) | — | 0.4365 / 0.4151 · 0.2277 / 0.1973 | 0.0115 · 0.0109 |
| ST(1.7)LA(1,1,200) | 0.945902 / 0.951278 (+0.57 %) | 0.0184 / 0.0071 | 0.0054 → 0.0073 |
- **최댓값이 매끈한 모양(TLA, ⊗G 가 있는 모양, DS)에서 Hw 가 정확**(fwhm ≤ 4.4e-4, asym ≤ 5e-4; 대칭 모양의 CasaXPS fwhm 차이와 같은 수준) → 이 경우 **정함**.
- 최댓값이 곡률 불연속(LA α≠β, n = 0) 또는 한쪽이 평평한 꼭대기(GL/SGL T(k): 높은 BE 쪽 1 − c·t³)에서는 CasaXPS 의 x̂ 가 꼬리 쪽으로 δ 만큼 밀린 것처럼 보이며,
  **fwhm 과 asym 두 값이 같은 δ 로 함께 맞는다** (Hw 의 형식은 맞고 x̂ 수치 추정만 다름). δ 는 F 에 비례(LA(1.53,2.43) F1·F2 같음) — 채널 격자 탓이 아님.
  δ 를 하나의 규칙(채널 포물선 꼭짓점, 수준 λ 이상 무게중심·중점, 작은 Gaussian 평활)으로 설명하려 했으나 LA·GL T·ST 를 함께 맞추는 것은 없음 → **경험값**:
  LA(α≠β,0) δ ≈ 0.0025 F (fwhm ≈ 2L·1.009), GL/SGL T(k) δ ≈ 0.011 F (k 1–2, m 0–30), ST 미정.
- 'Casa fwhm = 실제 FWHM' 가정(현재 casaToComp 의 TLA·DS·꼬리 경로)은 틀림: TLA(1.2,3,0) 에서 F 가 8 % 작게 나옴 → patch §2.

## 4. DS(α,n) — 위치 이동 정함, 넓이 미정
- 모양: P(x) = cos(πα/2 + (1−α)atan(2x_c))/(1+4x_c²)^{(1−α)/2}, 꼬리 높은 BE, F = wL, Gaussian FWHM = n·wL/144.87 → RMS 0.2 % (n 80–120 스캔에서 100 최적).
- **CasaXPS 위치 = INNO 원점(x = 0) + wL·|cot(π/(2−α))|** (높은 BE 쪽). 예측 285.0524, 판독 285.0523. (이것은 DS 원함수 꼭대기 위치 γ·|cot(π/(2−α))| 의 정확히 2 배 = γ 대신 wL 을 쓴 꼴.)
  ⊗G 뒤의 실제 최댓값(+0.0414 eV)이나 반높이 중점(+0.072)은 아님. 판독 1 개 → α 다른 판독으로 확인 필요(§6 R1).
- 넓이: CasaXPS Area / (진폭·wL) = **3.13695** (진폭 = P(0) 의 계수, 즉 모형 = 진폭·P((E−E0)/wL)). π(1 − 1.5e-3) 에 가깝지만(Lorentzian 넓이 차이와 같은 부호·크기)
  α = 0 이면 π/2 여야 할 값의 2 배라 이상함; π·Γ(1−α)/2^α = 3.1324 등 후보 여럿 → **미정**. INNO 창-넓이(15000) 대비 1.4343 배.
  유한 구간 정규화라면 ±500 F 정도가 필요(비현실적). INNO 의 '창 안 넓이' 정의는 그대로 두고 변환 계수만 표시할 것.

## 5. ST(μ)F — 모양 정함, 위치·폭 미정
- S(μ)F = F·[∫_{x_c}^∞F + A/μ] (계단 높은 BE 쪽), F = LA(1,1,200) (⊗G 뒤에 S) 가 최적. 다른 순서(S 뒤 ⊗G: RMS 610), 계단 반대(2508), A·μ(757) 모두 기각.
- CasaXPS 위치 = INNO 원점 + 0.0984 eV (= 0.197 F_L = 0.103 W). 최댓값(+0.105), 반높이 중점(+0.110), 중앙값(+0.145), DS 식 유추(원함수 꼭대기 ×2 = 0.070 eV) 모두 아님 → **미정**.
- fwhm 0.951278: 실제 FWHM 0.954768 (−0.37 %), 2L 0.945902 (+0.57 %) → 미정. 현재 casaToComp 는 바탕 FWHM 으로 F 를 정해 4.5 % 어긋남(patch §4).

## 6. 더 받으면 정해지는 판독 (기본 영역, Linear, Av.W 1, Marquardt ×3; Area·fwhm·pos·RMS(DoF)·asym 6 자리)
생성: `make_casaval3.js` 의 `cases` 에 아래 이름:문자열을 넣고(`fw` 는 _F1_2 만 1.2) 실행 → `/tmp/k3/*.vms` 를 casaval 로. 각 파일은 같은 문자열로 맞춤.
- **R1 DS 넓이·위치**: `K_DS0_100_F1:'DS(0,100)'`, `K_DS02_100_F1:'DS(0.2,100)'`, `K_DS01_0_F1:'DS(0.1,0)'` + 기존 K_DS01_100_F1 을 **영역 End 290.03** 으로 줄여 다시 맞춤.
  예측(위치): DS(0) 이동 0, DS(0.2) ≈ +0.176·wL(≈ +0.10 eV), DS(0.1,0) ≈ +0.083·wL(≈ +0.07 eV). 넓이: α 에 따른 Area/(진폭·wL) 와 창 의존 여부가 정규화를 정함.
- **R2 ST 위치·폭**: `K_ST17_LA110_F1:'ST(1.7)LA(1,1,0)'`, `K_ST3_LA11200_F1:'ST(3)LA(1,1,200)'`, `K_ST17_GL30_F1_2:'ST(1.7)GL(30)'` — pos 이동이 μ·바탕·n 에 어떻게 변하는지.
- **R3 폭 숫자의 '쪽'**: `K_LA243_153_0_F1:'LA(2.43,1.53,0)'` (꼬리 낮은 BE). Hw(가파른 쪽 = 이제 높은 BE) 이면 fwhm ≈ 0.871, '낮은 BE 쪽' 이면 ≈ 1.147; asym 부호.
  `K_LA2_3_30_F1:'LA(2,3,30)'` (꼭대기 매끈) → Hw 예측 fwhm 0.9301 (정확해야 함), asym 0.1306.
- (선택) GL(30)T(1)·GL(0)T(2) 를 V2_GL30T2_F1.2 에 맞춘 행의 **fwhm** 숫자 — 예측 0.98519 F·1.02263 F (F 는 그 맞춤의 F, δ 0.011 규칙).

## 7. 적용 상태 (2026-10-29) — `src/xps_core.js` §16·`fitStats`, `src/xps_ui.js` casaEq1, 시험 `test/test_casa_v3.js`
측정으로 뒷받침되는 것만 넣었다. LM 멈춤 규칙(SPEC §3.3b, `LM_STALL`)은 건드리지 않음.
| 항목 | 상태 | 내용 |
|---|---|---|
| §2 DoF | **적용 (정함)** | `fitStats(...).casa.dof = #{i0+1…i1−1 : y ≠ B} − p` (B = 정적 영역 배경, `|y−B| > 1e-12·max(1,|y|)`; 능동 배경은 모든 채널). `casa.n` = 센 채널 수. `covariance()` 의 DoF 는 그대로(CasaXPS 오차 막대 DoF 미판독). 12 행 모두 판독과 같음 (151·156·337·338). |
| §3 fwhm·asym | **적용** | `trueFWHM`·`s2FWHM` 가 반높이 점 `lo`·`hi` 를 돌려줌. 새 `casaWidthAsym(c)` → `{fwhm, asym, exact, lo, hi, xmax, note}`: fwhm = 2(L+δ), asym = 1−(L+δ)/(H−δ). 꼭대기 매끈(TLA, DS, LA n>0, GL/SGL(m,n) 꼬리 없음, 대칭 모양) δ = 0 → `exact:true`. 꺾인/평평한 꼭대기 → 경험 δ, `exact:false` ('≈'): LA(α≠β, n=0) 0.0025·F, GL/SGL T(k)/PHI·AGT(tu fwhm) 0.011·F; ST 는 fwhm = 수식어 붙인 모양의 실제 FWHM, asym 은 δ 0.0054·F (판독 1). 일반 격자 경로(LA·GLP·AGT …)는 `refineHalfMax` 로 최댓값(황금분할)·반높이 점(이분)을 1e-12·W 까지 다듬음. |
| §3 보정 | 기록 | 기준 `ref_casa_v3.py` 의 `peak()` 는 LA 꺾인 꼭대기에서 스플라인 탓에 x̂ 를 +6.6e-5·F 옮겨 잡는다 → 정확한 반높이 점(LA(1.53,2.43): L 0.287289, H 0.378522; LA(2,3): 0.254913, 0.321797)으로 다시 구한 δ: LA(2,3,0) 0.00234 (fwhm)·0.00254 (asym), LA(1.53,2.43,0) 0.00276·0.00289. δ = 0.0025 는 그대로 두되 비대칭 LA 의 asym 허용을 2e-3 으로 (patch 의 'asym 0.2283' 은 0.2294 가 맞음). |
| §3 compToCasa | **적용** | 모든 결과에 `asym`, `posShift`(CasaXPS 위치 − INNO 원점), `fwhmExact` 추가. GL(m,n)/SGL(m,n) 꼬리 없음·TLA(α,μ,n): `exact` = n 정수 (n 규칙 '측정'으로 고침). 꼬리·비대칭 LA·AGT GL(p)T(k): fwhm·asym 은 경험값, `exact:false`. DS: fwhm·asym 정함, `posShift = wL·|cot(π/(2−α))|` ('measured (1 reading)'), 넓이 때문에 `exact:false`. ST: `posShift:null` (미정), `exact:false`. |
| §3 casaToComp | **적용** | 폭: CasaXPS fwhm → F (DS 는 wL) 를 단위 폭의 `casaWidthAsym` 으로 나눠 역산 ('Casa fwhm = 실제 FWHM' 가정 제거: TLA(1.2,3,0) F 8 % 작던 것 → 4.5e-5). 수식어를 먼저 붙이므로 ST 폭은 ST 모양으로 (4.5 % → 3.7e-3, 남은 차이는 미정 fwhm 숫자). GL(p)T(k) 의 'F = fwhm 숫자' 제거 (−1.6 % → 2e-4). DS 원점 = CasaXPS 위치 − wL·|cot(π/(2−α))| (`r.posShift`). `exact:true`: GL(m,n)·SGL(m,n)(꼬리 없음)·TLA. 노트 'assumed' → 'measured'. 새 선택 `o.innoUnits` (fw = 실제 FWHM, be = INNO 원점 — 화면에서 모양만 바꿀 때; 이전 화면 동작 유지). |
| §3 casaModelToINNO | **적용** | MFWHM 비·한계는 DS 에서 wL 로 (전에는 F 가 없어 1), FWHM link → wL link. 위치 한계 = CasaXPS 위치 − posShift; 위치 link 는 여전히 재현하지 않되(오프셋 표기 미확인), 짝의 posShift 가 다르면 알림에 'INNO 오프셋 = CasaXPS 오프셋 − 이동 차' 를 적음. |
| 화면 | **적용** | 성분 카드의 'CasaXPS:' 줄: 표기, fwhm, asym(≠ 0), DS 의 'pos INNO + Δ eV'; 정한 값은 그대로, 경험·미정 값은 '≈ ' 를 앞에 붙임 (ui_xps_shapes2 4c). |
| §4 DS 넓이 | **미정 (그대로)** | INNO 넓이 = 창 안 넓이. CasaXPS Area / INNO 창 넓이 = 1.4343 (DS(0.1,100)); 노트에 'not pinned'. R1 판독 대기. |
| §5 ST 위치·fwhm | **미정 (그대로)** | ST(1.7)LA(1,1,200): CasaXPS pos − INNO 원점 = +0.0984 eV, fwhm 0.951278 vs 실제 FWHM 0.954768 — 노트·`exact:false`·`posShift:null`. R2 판독 대기. |
| R3 (낮은 BE 꼬리) | 미정 | 'fwhm 은 낮은 BE 쪽' 그대로 (모든 판독이 높은 BE 꼬리). |

시험 `node test/test_casa_v3.js` (약 1 초, 69 항목): 12 행(K_* 8 + V2 4)을 INNO CasaXPS 규약 피팅으로 맞춰 (a) scipy 기준과 넓이·위치·폭·RMS ≤ 1e-4, DoF 규칙;
(b) CasaXPS 판독과 넓이 ≤ 1e-3 (DS 제외), 위치 ≤ 2e-4 eV (DS 는 posShift 뒤; ST 제외), RMS ≤ 3.5 % 또는 20 CPS (CasaXPS RMS 바닥), DoF 같음,
fwhm ≤ 5e-4 (정함) · 1e-3 (경험 δ) · 5e-3 (ST), asym ≤ 1.2e-3 (비대칭 LA 2e-3, ST 3e-4), `fwhmExact` 표시; asym 만 읽은 GL(30)T(1)·GL(0)T(2) ≤ 2e-3;
표기 왕복 (INNO → Casa → INNO 항등, CasaXPS 판독 → 맞춘 INNO 폭·원점), DS 원점 284.9994 (patch 기대 284.9996 ± 2e-4), casaModelToINNO DS 폭 비·위치 한계.
기존 `test_xps_casa.js` §5 의 비대칭 LA·GL(30)T(2) 항목은 '미확인 (fwhm null)' 대신 경험값 대조로, V2 DoF 는 CasaXPS 판독(337)과 비교로 바꿈; `test_shapes2.js` 의 GL(30)T(2) ≡ AGT 비교는 같은 F 로.

