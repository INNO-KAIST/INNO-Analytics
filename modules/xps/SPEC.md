# XPS 모듈 수식 명세 (단일 원본 · v1.7)

이 문서가 정의다. `src/xps_core.js`(JS 엔진)와 `test/ref/`(독립 Python 기준 구현)는 이 문서만 보고 따로 만든다.
코드와 문서가 다르면 문서를 고치고 코드를 맞춘다.

## 0. 표기
- 에너지 축은 결합에너지(BE, eV). 내부 배열은 **BE 오름차순**. 화면만 내림차순(XPS 관례)으로 그린다.
- `x = E − E0` (E0 = 성분 중심 BE). x>0 은 더 높은 BE 쪽.
- 모든 lineshape 는 별도 언급이 없으면 **실수 전체(ℝ)에서 넓이 1** 이다. 성분 세기 = `A · f(x)` (A = 넓이, 단위 CPS·eV).
- `ln2 = ln 2`. `erfcx(z) = exp(z²)·erfc(z)`.

## 1. Lineshape

### 1.1 기본 함수
- Gaussian: `G(x;F) = (2/F)·√(ln2/π)·exp(−4 ln2 · x²/F²)`
- Lorentzian: `L(x;F) = (2/(πF)) / (1 + 4x²/F²)`

### 1.2 `GLS` — Gaussian–Lorentzian **합** (넓이 기준 pseudo-Voigt; Avantage "Sum"·CasaXPS SGL 과 같은지는 §13·§13b)
매개변수: F (FWHM, eV), m ∈ [0,1] (Lorentzian 비율).
`f = (1−m)·G(x;F) + m·L(x;F)`. FWHM 은 정확히 F.

### 1.3 `GLP` — Gaussian–Lorentzian **곱** (Avantage "Product", CasaXPS GL(m))
매개변수: F, m ∈ [0,1].
`g(x) = exp(−4 ln2 (1−m) x²/F²) / (1 + 4m x²/F²)`, `f = g / (F·J(m))`
`J(m) = ∫ exp(−a t²)/(1+b t²) dt`, a = 4 ln2 (1−m), b = 4m.
- b>0: `J = (π/√b) · erfcx(√(a/b))`
- m=0: `J = √(π/a)`
주의: 0<m<1 이면 실제 FWHM ≠ F. 화면에는 실제 FWHM(수치로 구한 반치폭)을 따로 보인다.

### 1.4 `V` — Voigt (Gaussian ⊗ Lorentzian, 진짜 합성곱)
매개변수: wG (Gaussian FWHM), wL (Lorentzian FWHM).
σ = wG / (2√(2 ln2)), γ = wL/2, `z = (x + iγ)/(σ√2)`,
`f = Re[w(z)] / (σ√(2π))`, w = Faddeeva 함수.
극한: wG→0 이면 L(x;wL), wL→0 이면 G(x;wG).
표시용 전체 FWHM: Olivero–Longbothum `0.5346 wL + √(0.2166 wL² + wG²)`(실측 최대 상대오차 2.4e-4), 표에는 수치 반치폭.

### 1.5 `DS` — Doniach–Šunjić ⊗ Gaussian (금속의 비대칭)
매개변수: α ∈ [0, 0.45], wL (γ = wL/2), wG.
`ds(x) = cos(πα/2 + (1−α)·atan(−x/γ)) / (x² + γ²)^((1−α)/2)` — 꼬리는 **높은 BE 쪽**(x>0).
`p = ds ⊗ G(·;wG)` (wG=0 이면 합성곱 없음).
α>0 이면 ℝ 적분이 발산하므로 **넓이는 피팅 영역 [Emin, Emax] 안의 적분**으로 정의: `f = p / ∫_{Emin}^{Emax} p dE`.
(이 성분의 넓이 A 는 "영역 내 넓이"다. 결과표에 표시한다.)

### 1.6 `LA` — 비대칭 Lorentzian 거듭제곱 ⊗ Gaussian (CasaXPS LA 와 같은 형태)
매개변수: F, α (**높은 BE 쪽** 지수), β (낮은 BE 쪽 지수), α,β ∈ [0.55, 20], wG ≥ 0.
`b(x) = (1 + 4x²/F²)^(−α)` (x≥0), `(1 + 4x²/F²)^(−β)` (x<0)
(v1.1: 금속 문헌 템플릿 — 예: Ni metal LA(1.1, 2.2, ·), Biesinger et al., Appl. Surf. Sci. 257, 2717 (2011) — 이 높은 BE 꼬리를 주려면
α 가 높은 BE 쪽이어야 해서 v1 의 좌우를 바꿨다. CasaXPS 의 세 번째 인자는 정수 Gaussian 폭 지수이고 여기서는 wG(eV)다.)
ℝ 정규화 상수(해석식): `N = (F/4)·√π·[Γ(α−½)/Γ(α) + Γ(β−½)/Γ(β)]`
`f = (b/N) ⊗ G(·;wG)`. α=β=1, wG=0 이면 L(x;F).

### 1.7 `GLT` — GL 합 + 지수 꼬리 (Avantage 의 tail 과 같은 목적, 식은 여기서 정의)
매개변수: F, m, T ∈ [0,1] (꼬리 섞음 비율), τ > 0 (꼬리 길이, eV).
`e(x) = (1/τ)·exp(−x/τ)` (x ≥ 0), 0 (x<0) — 높은 BE 쪽 꼬리.
`f = (1−T)·GLS(x;F,m) + T·(GLS ⊗ e)(x)`. 각 항 넓이 1 → f 넓이 1.

### 1.8 합성곱 계산(DS, LA(wG>0), GLT)
균일 격자, 간격 h = 0.04·2^(−k) eV 중 `h ≤ min(관련 폭)/16` 인 가장 큰 값(하한 0.0005). **양자화**하는 까닭: 폭 매개변수가
조금 바뀌어도 마디가 그대로여야 목적함수가 매끄럽고 수치 Jacobian 이 맞다(연속 h 는 LM 정지·공분산 오류를 냈다 — 기준 구현 대조로 확인).
마디는 상대좌표 `x = k·h` (k 정수, 성분 중심이 항상 마디), 범위 [Emin − P, Emax + P], P = 8·(가장 큰 폭) + 5·σG (GLT 는 max(P, 12τ + 8F)).
Gaussian 핵은 ±5σ 절단 후 이산 합으로 정규화. 데이터 점으로는 **Catmull–Rom(C¹ 3차) 보간**. DS 영역 적분은 마디 Simpson + 끝 조각 3점 Simpson.
지수 꼬리 합성곱은 마디 사이 (v_i, v_{i−1}, v_{i−2}) 2차 보간의 정확 적분 점화식(v1.1).
LA(wG=0)·DS(wG=0)는 격자 없이 해석식으로 계산한다(DS 영역 정규화는 원시함수 `−r^α sin(αθ − πα/2)/α`, θ = atan(−x/γ), α→0 이면 π/2 + atan(x/γ)) — 격자는 Gaussian 합성곱이 필요할 때만.
기준 구현은 이와 무관하게 정확한 수치 적분(quad 또는 매우 세밀한 격자)으로 값을 낸다 → 허용 오차 1e-3(상대, 최댓값 기준)로 이산화 오차까지 검증.

### 1.8b `MEAS` — 측정 선모양 (v1.7, 정의는 §12)
매개변수: wG ≥ 0 (추가 Gaussian 넓힘 FWHM, eV). 모양 = 기준 스펙트럼에서 만든 넓이 1 template `R(x)` (§12.1), `f = R ⊗ G(·;wG)`.
성분 BE = template 의 기준 위치(eref, 배경 뺀 신호의 최댓값)가 놓이는 자리. doublet 은 쓰지 않는다(측정 모양이 이미 짝을 포함).

### 1.10 `AGT` — GL 곱 + Avantage 꼬리 함수 (v1.7)
정의: Thermo Avantage Data System manual 'Tail Function' (P. M. A. Sherwood, Appendix 3 in D. Briggs, M. P. Seah (eds.), *Practical Surface Analysis* (Wiley);
R. O. Ansell et al., *J. Electroanal. Chem.* 98, 79 (1979)). 매개변수 F, m(= L/G Mix/100), ET ∈ [0.001, 1] (/채널), TM, CT ∈ [0,1]; 고정 속성 ch = 채널 폭(eV, 성분을 만들 때 데이터 간격).
`GL(x) = exp(−4ln2(1−m)x²/F²)/(1+4m x²/F²)` (§1.3 의 정규화 전 곱 — Avantage 의 β = F/2), `T(x) = TM·CT + (1−TM)·exp(−ET·x/ch)`,
`h(x) = GL(x)` (x ≤ 0, 낮은 BE 쪽 꼬리 없음), `h(x) = GL + (1 − GL)·T` (x > 0). ET = 1 은 사실상 꼬리 없음, ET → 0 은 무한 꼬리.
상수 꼬리(TM·CT > 0)는 ℝ 넓이가 무한이므로 **넓이는 피팅 창 안 적분**(§1.5 DS 와 같은 규칙): `f = h / ∫_{Emin}^{Emax} h`, 적분은 [lo,0]·[0,hi] 각각 Simpson(간격 ≤ min(F/40, ch/4)).
검증: quad 대조 ≤ 2e-8; 합성(참 ET 0.06·TM 0.2·CT 0.08, 상수 배경 자유, 30 회) ET·넓이·F 95 % 구간 적중 30/30 · 30/30 · 30/30. 주의: 창 끝 직선 배경은 상수 꼬리를 배경으로 먹는다.

### 1.9 Doublet (스핀–궤도 쌍)
주 성분(낮은 BE, 큰 j) + 짝 성분: 중심 `E0 + Δ`, 넓이 비 `r = A_짝/A_주`, 폭 배율 `q`, 모양 매개변수(m, α, β, T, τ)는 같다.
q 는 **수명(Lorentzian) 폭에만** 곱한다: V·DS 는 wL, GLS·GLP·GLT·LA 는 F. wG(장비 Gaussian 폭)는 두 짝에 공통 — Coster–Kronig 폭 넓어짐은 core-hole 수명 효과다(v1.1, v1 은 wG 에도 곱했다).
매개변수 A 는 **쌍 전체 넓이**: `A_주 = A/(1+r)`, `A_짝 = A·r/(1+r)`.
기본 r = (2j₂+1)/(2j₁+1): p 1/2, d 2/3, f 3/4. 기본 q = 1 (Coster–Kronig 로 넓어지는 2p1/2 등은 q 를 풀어 피팅).

## 2. Background (영역 [Emin, Emax], 데이터 y)
끝점: `y_lo` = 가장 낮은 BE 쪽 n 점 평균, `y_hi` = 가장 높은 BE 쪽 n 점 평균(기본 n=5). 끝점 에너지 `E_lo`, `E_hi` = 그 n 점 BE 평균.
적분은 데이터 격자 위 사다리꼴 누적.

- `linear`: (E_lo, y_lo)–(E_hi, y_hi) 직선.
- `shirley` (Proctor–Sherwood 반복): `B(E) = y_lo + (y_hi − y_lo)·Q(E)/Q(Emax)`, `Q(E) = ∫_{Emin}^{E} (y − B) dE'`.
  B₀ = y_lo 에서 시작, `max|ΔB| < 1e-7·(|y_hi − y_lo| + 1e-12)` 또는 200 회까지.
  (y_hi − y_lo 가 0 이면 B = y_lo.) **200 회 안에 수렴하지 않으면 그대로 쓰되 `converged=false` 로 표시하고 화면에 경고**한다
  — 봉우리 밖 바탕이 기울어 있으면(창 양끝이 봉우리와 무관한 경사 위) 고정점이 없을 수 있다. 이때는 `active`(+ 기울기 s)를 권한다.
- `smart` (Shirley + "데이터보다 높지 않음" 제약): shirley 와 같은 반복이되 매 회 `B ← min(B, ys)`.
  ys = 국소 2차 다항식 평활(중심 창 2h+1, h = max(2, round(0.25 eV / 간격)), 가장자리는 가능한 점만 쓴 비대칭 창).
  적분식의 피적분은 `(y − B)` (자르지 않음 — 잡음에 치우치지 않게).
- `tougaard` (universal cross-section, 비반복):
  `K(T) = T/(C + T²)²` (U2, C=1643 eV²) 또는 `T/((C − T²)² + D·T²)` (U3, C·D 지정).
  `I(E) = ∫_{Emin}^{E} K(E − E')·(y(E') − y_lo) dE'`, `B(E) = y_lo + s·I(E)`,
  v1.7: 문턱 T0 ≥ 0 (기본 0) — K(T) = 0 (T ≤ T0) (CasaXPS U(x; B, C, D, T0), Fairley et al. 2023 SI Eq. 9; 절연체에서 손실 시작을 늦춤).
  s 는 `B(E_hi) = y_hi` 가 되게 정한다(I(E_hi) 는 E_hi 에서 선형 보간). **s < 0 이면 0 으로 두고 경고**(음의 손실 단면적은 비물리적 —
  y_hi < y_lo 인 좁은 창. Tougaard 는 봉우리 위로 30–50 eV 이상의 창을 전제한다).
- `active` (피팅 안에서 함께 정하는 Shirley): `B(E) = c₀ + k·∫_{Emin}^{E} P(E') dE' + s·(E − Emin)`,
  P = 모든 성분 합. c₀, k(≥0) 는 피팅 매개변수, s 는 선택(기본 0 고정).
  출발값: 넓이 초기화(NNLS)는 **직선 배경**(linear) 위에서 한다 — 정적 Shirley 가 수렴하지 않은 창에서 넓이 출발이 0 이 되고 k 가 발산(k·Q 와 넓이의 퇴화)하던 결함 때문.
  k 의 위 경계 `k_hi = max(1.5·(max y − min y)/Q(Emax), 2·k₀)` (Q 는 출발 성분 합의 누적 적분) — 퇴화 방향으로 달아나지 않게.
- `spline` (v1.7, CasaXPS 'Cubic Spline' 과 같은 목적): 매듭 = 끝점 둘 `(E_lo, y_lo)`, `(E_hi, y_hi)` (`ends:false` 면 없음) + 사용자 매듭 `{e, y?}`.
  y 를 주지 않은 매듭(자동)은 e 에 가장 가까운 데이터 점 i 의 ±kh 점 평균(기본 kh = 2, 창 끝에서 자름). 같은 e(1e-9 eV 안)의 끝점·사용자 매듭은 사용자 것.
  매듭 (x_k, y_k) (오름차순) 을 지나는 **자연 3차 spline** (양 끝 S″ = 0): 마디 2차 도함수 M 은 삼중대각
  `h_{k−1} M_{k−1} + 2(h_{k−1}+h_k) M_k + h_k M_{k+1} = 6[(y_{k+1}−y_k)/h_k − (y_k−y_{k−1})/h_{k−1}]`, M_0 = M_{K−1} = 0;
  `S(x) = A y_k + B y_{k+1} + [(A³−A) M_k + (B³−B) M_{k+1}] h_k²/6`, A = (x_{k+1}−x)/h_k, B = 1 − A.
  첫·끝 매듭 밖은 그 점의 기울기로 **직선 연장**(3차 연장은 창 밖에서 발산). 매듭 1 개 = 상수, 0 개 = 0.
  화면: '+ Knot' 뒤 그래프 클릭 = 매듭(자동 높이), 끌기 = 위치·높이(높이는 고정됨), 두 번 누르기 = 지움, 'Auto heights' = 끈 높이 되돌림.
  주의: 매듭은 봉우리 밖 배경 위에 둔다 — 봉우리 밑에서만 보이는 배경 모양(예: 계단)은 어떤 보간으로도 정할 수 없다(합성 시험에서 tanh 계단 배경은 −5…−10 % 넓이 오차).
  기준: Press et al., *Numerical Recipes* 3rd ed. §3.3; C. de Boor, *A Practical Guide to Splines* (1978).
  불확도: 매듭 높이(±kh 점 평균)의 잡음은 공분산 σ 에 들어가지 않는다(정적 배경 모두 같음). 합성 시험(휜 배경, 매듭 2 + 끝 2, 16 회)에서 넓이 치우침 없음(평균 +0.8 %),
  그러나 |오차|/σ통계 최대 12 → 보고는 §11 의 매듭 흔들기 '가능 범위'로(같은 시험에서 15/16 포함). 매듭 높이 잡음을 줄이려면 kh 를 늘린다.
- `none`: B = 0.

### 2b. CasaXPS 호환 영역 규약 (v1.7, 선택 `bg.conv = "casa"`, Linear·Shirley 에만; 기본은 위 INNO 규약 그대로)
CasaXPS 2.3.27 실측(연구실 PC, 2026-10-28 — test/ref/casa_readings.json, 자료 test/data/casaval)에서 역으로 정한 규약. 창 = 채널 i0 … i1 (CasaXPS 영역 Start/End).
- 끝점: `y_lo` = 채널 i0+1 을 중심으로 ±w 채널 평균, `y_hi` = 채널 i1−1 중심 ±w (w = CasaXPS 'Av. Width', 정수; 평균 창은 피팅 창 밖 채널도 포함 → 전체 스펙트럼에서 계산),
  직선은 E[i0+1], E[i1−1] 에 고정.
- Shirley: 수렴한 Proctor–Sherwood 이되 누적이 **배타적 직사각 합** `Q_i = Σ_{j=i0+1}^{i−1} (y_j − B_j)` (사다리꼴 아님), `B_i = y_lo + (y_hi − y_lo) Q_i / Q_{i1−1}`, i ≤ i0+1 은 y_lo, i ≥ i1−1 은 y_hi. 수렴 규칙은 §2 와 같음.
- 영역 숫자(CasaXPS Regions 표의 Area · Max Height · Min Height): 채널 i0+1 … i1−1 에서 (y − B) 의 사다리꼴 적분·최대·최소 (`casaRegionStats`).
- 결과: Max·Min Height 는 화면 자리수까지 같음(Linear w 0·1·3, Syn1·Syn2), 넓이 상대 ≤ 7e-5 (Linear; CasaXPS 넓이 적분 규칙은 남은 2–7e-5 차이만큼 미확인),
  Shirley(무잡음 Syn2) 넓이·높이 같음(1.7e-6); 실측 C1s Shirley 는 넓이 1.2e-3·높이 1–2 CPS 차 (원인 미확인);
  **기운 창(Syn1: 선형 경사 + Gaussian, Shirley 고정점 없음)은 CasaXPS 12545.3 vs 여기 10389.7** — CasaXPS 가 수렴하지 않는 반복을 어디서 멈추는지 모름(반복 1–2 회 사이 값). INNO 규약의 화면 경고(§2)는 그대로.
- Smart·Tougaard·Spline 은 CasaXPS 규약을 쓰지 않는다(Tougaard 는 아래 §13).

### 2c. Avantage 규약 Shirley (v1.7, `bg.conv = "avantage"`, Shirley 에만)
Thermo Avantage manual 'Shirley Background': 반복 수 상한 `it`(기본 10 — 설명서 권장; 1 = 단순 한 번 계산), 'Convergence limit' = 반복 사이 넓이 변화의 상한 `cl`(상대, 기본 0 = 모두 반복),
**높은 BE(낮은 KE) 끝이 낮은 BE 끝보다 높지 않으면 Linear 로 바꿈**(설명서). 끝점·적분은 §2 그대로(설명서에 끝점 평균 규칙은 없음 — Avantage 대조 대기). 반복 B_{k+1} = y_lo + (y_hi−y_lo)Q_k/Q_k(Emax), B_0 = y_lo.
**Smart 비교**(설명서 'Smart Background'): 둘 다 'Shirley + 배경 ≤ 데이터'. 차이 — Avantage 는 배경이 데이터보다 크면 '배경이 데이터를 따라가고 Shirley 시작점을 창 안으로 옮기며' 넓이를 최대로 하려 함;
여기 §2 smart 는 반복마다 점별로 `B ← min(B, 평활 데이터)` (시작점을 옮기지 않음, 평활한 데이터와 비교). 시작점 이동 규칙·평활 여부가 공개되지 않아 Avantage Smart 를 그대로 재현하지 않는다(대조 대기).

## 3. 피팅
### 3.1 목적함수
`χ² = Σ_i (y_i − μ_i)² / v_i`, μ = B + Σ 성분.
분산 모형(준-Poisson): `v_i = max(ŷ_i, v_min) / c`.
- c = 계수 환산(1 CPS 가 몇 count 인가) = dwell × scans (VGD 처럼 알 때). 모르면 c = 1 (척도는 φ 로 흡수).
- ŷ: 1차는 평활 데이터(2절의 ys), 수렴 뒤 ŷ = μ 로 바꿔 한 번 더 피팅(IRLS 2회, Pearson χ²).
- v_min = max(1/c, 1e-6·max ŷ).
- **일정 가중**(`weight = "const"`, v1.6): 배경을 이미 뺀 자료처럼 세기가 0·음수 근처라 위 분산 모형이 맞지 않을 때. `v_i = v` (상수 — 2차 차분 `(y_i − (y_{i−1}+y_{i+1})/2)²/1.5` 의 중앙값 ÷ 0.4549), IRLS 없음. 절대 척도는 φ 로 다시 정해지므로 σ 는 v 와 무관. AICc·BIC 의 χ²_d 도 같은 v.
  읽기에서 음수 세기가 있거나 10 % 넘는 점이 최댓값의 0.5 % 이하이면 '배경을 뺀 자료'로 보고 Background None + 일정 가중을 기본으로 둔다(경고 표시).
### 3.2 매개변수 제약
각 매개변수: `free`(상·하한 [lo,hi]), `fixed`, `link`(= 기준성분의 같은 매개변수 × k + b).
링크 기본: BE 는 `+b`(간격), 넓이·폭은 `×k`(비).
넓이 ≥ 0 은 항상. 한쪽만 유한한 한계는 √ 변환, 양쪽 유한은 sin 변환(MINUIT 방식)으로 무제약화.
### 3.3 최적화
Levenberg–Marquardt (Nielsen 감쇠 갱신), 수치 Jacobian(중앙차분, 넓이는 해석 열).
수렴: 상대 χ² 변화 < 1e-10 이 3회 연속, 또는 최대 500 회.
선택(v1.7): **Simplex**(Nelder–Mead, §14) — `prob.optimizer = "simplex"`. 같은 변환·같은 목적함수·같은 IRLS 2회, 공분산은 최솟값에서 §3.4 그대로.
초기 넓이는 고정한 비선형 매개변수에서 NNLS 로 구한다.
다중 시작: 비선형 free 매개변수를 [lo,hi] 안에서 Latin hypercube 로 뽑아 N 회(기본 24) → **1차 χ²**(분산 = 평활 데이터, 모든 시작점에 같은 가중)가 최소인 해 채택, 그 해와 Δχ² < 1 인 시작의 비율을 "해 유일성" 지표로 보고. 2차 분산 v=μ₁ 은 시작점마다 달라 2차 χ² 로 시작점을 고르면 안 된다(v1.1 — 두 번째 기준 구현 대조로 발견).
### 3.4 불확도
- 공분산: `C = (Jᵀ V⁻¹ J)⁻¹ · φ`, `φ = χ²/(N − p)` (p = free 매개변수 수). 한계에 닿은 매개변수는 고정으로 보고 제외(표시).
- 성분 값(BE·FWHM·넓이 …)이 어떤 자유 매개변수에도 의존하지 않으면(고정, 또는 고정된 것에 연결) σ = null — 화면 '± —', CSV 빈칸(0 이 아님).
- 잔차 자기상관: `ρ₁ = Σ r_i r_{i+1} / Σ r_i²`, r = (y−μ)/√v. ρ₁>0 이면 AR(1) 근사 보정 배수 `√((1+ρ₁)/(1−ρ₁))` 를 σ 에 곱한 값도 보인다.
- Monte Carlo: 최적 μ 에 AR(1) 잡음(분산 φ·v_i, 상관 ρ₁⁺ = max(ρ₁,0)) 을 더한 합성 데이터 M 개(기본 200)를 다시 피팅 → 표준편차·2.5/97.5 백분위.
### 3.5 적합도
χ²_ν = φ, R², 잔차 RMS(정규화), Durbin–Watson, ρ₁, 정규화 잔차 연속 부호(runs) z-값,
`AICc = N ln(χ²_d/N) + 2p + 2p(p+1)/(N−p−1)`, `BIC = N ln(χ²_d/N) + p ln N`. 여기서 χ²_d 는 **데이터만으로 정한 분산**(1차 v = 평활 데이터)으로 계산한 χ² —
2차 분산 v = μ 는 모형마다 달라 모형 비교 항이 섞이기 때문(v1.1). 같은 데이터·같은 영역·같은 c 사이에서만 비교.

## 4. 정량
`n_i = A_i / (S_i · T(KE_i, PE) · λ(KE_i))`, `at%_i = 100·n_i / Σ n`, `KE = hν − BE`.
- S: 감도계수 라이브러리(Scofield 1976 단면적, C1s = 1) 또는 사용자 값.
- T (transmission): `none`=1 / `power`: `(PE/KE)^x` / `thermo`: `ln T = Σ_{k=0}^{3} a_k (ln(KE/PE))^k` (VGD 의 계수 — 형태는 추정, 검증 필요 표시) / `table`(KE–T 쌍, log–log 보간, 표 밖은 끝값).
- λ (IMFP): `none` / `power`: KE^n (기본 n=0.6) / `tpp2m`:
  `λ = E / {Ep² [β ln(γE) − C/E + D/E²]}` (Å, E=KE eV), `Ep = 28.8 √(Nv ρ / M)`,
  `β = −0.10 + 0.944/√(Ep² + Eg²) + 0.069 ρ^0.1`, `γ = 0.191 ρ^(−1/2)`, `U = Nv ρ / M`, `C = 1.97 − 0.91U`, `D = 53.4 − 20.8U`.
- 불확도: 영역별 넓이는 서로 독립, at% 는 델타 방법(또는 MC 표본 결합).
- **T 의 Avantage 형태**(v1.7, 선택 `model: "avantage"`; Avantage manual 'Transmission Function'): `log10[A/(PE·XSF)] = a + b·x + c·x² + d·x³`, x = log10(KE/PE) → `T = PE·10^{Σ aₖ xᵏ}`.
  지금 기본(`thermo`: ln T = Σaₖ(ln(KE/PE))ᵏ, PE 곱 없음)과 다르다 — 1차 항만 있으면 두 형태의 기울기가 같지만 2·3차 항과 PE 곱(서로 다른 PE 를 섞을 때)이 다르다.
  화면 정량 탭 'Avantage form (log10)' 으로 고름(기본 꺼짐). VGD 계수가 어느 형태인지는 **같은 자료의 Avantage at% 대조로 확정**해야 한다(io 의 기본을 바꾸는 것은 그 뒤).
- **Monte Carlo**(v1.7, CasaXPS 식): 넓이를 독립 `A_i + σ_{A,i}·z` (z ~ N(0,1)) 로 M 번(화면 2000) 뽑아 at% 를 다시 계산 → 평균·sd·2.5/97.5 백분위.
  검출 안 됨(n.d.) 집합은 기준값으로 정해 고정. σ_A 가 크면(상대 > 10 %) at% 분포가 비대칭이라 델타 ±σ 대신 이 구간을 본다.
- **검출 안 됨**(v1.6): `A ≤ 0` 이거나 σ_A 를 알 때 `A < 2σ_A` 이면 n.d. — at% 를 내지 않고 Σ 에서 뺀다(음의 넓이가 다른 원소의 at% 를 부풀리지 않게).
- **행 기본 포함**(화면): 원소마다 Scofield(Al) 가 가장 큰 선 하나만(예: Co 2p 와 Co 3p 중 Co 2p). 사용자가 체크한 것은 그대로.
- **창 안 간섭 경고**(화면): 포함된 원소 + 같은 시료 survey 에서 채택한 원소의 주선(Scofield 최대 궤도의 subshell 라이브러리 BE)이 다른 원소 행의 피팅 창(보정 BE) 안이면 경고(예: Co 2p 창 안 Ba 3d).

## 5. 대전 보정
시료(같은 측정 묶음)마다 Δ 하나: `BE_보정 = BE_원래 + Δ`. 피팅은 원래 BE 로 저장하고 화면·표·라이브러리 대조는 보정값.
기준: C 1s C–C/C–H = 284.8 eV(기본, 바꿀 수 있음) 또는 직접 Δ 입력.
기준 성분(v1.6): Reference 에서 고른 성분, 아니면 자동 = 넓이가 C 성분 합의 10 % 이상인 것 중 가장 낮은 BE. 그래프에서 클릭해 둔 성분은 쓰지 않는다(쓰려면 'Use selected component' — 명시적).
고른 기준이 자동 기준과 다르면 알림에 경고. 기준 C 1s 가 피팅돼 있지 않으면 그 스펙트럼(화면 전환 없이)을 자동 제안 → 피팅한 뒤 맞춘다.
템플릿은 보정 BE 로 저장: 성분 BE 와 free BE 한계를 함께 +Δ, 불러올 때 함께 −Δ(`shiftModelBE`) — 왕복은 어떤 Δ 에서도 항등.
피팅·자동 제안은 계산 시작 때의 모형 판(ver)을 기억하고, 계산 중 사용자가 모형을 바꾸면 결과를 버린다.

## 6. 잡음 추정(데이터만으로)
주의(VGD 실측): VGData 는 정수가 아니다(검출기 채널→에너지 보간). 2차 차분 추정은 이웃 점 양의 상관 때문에 분산을 낮게 본다.
그래서 c(dwell·scans)는 분산의 **모양**(∝ 세기)에만 쓰고 절대 척도는 언제나 φ = χ²/(N−p) 로 다시 정하며, ρ₁ 보정·AR(1) Monte Carlo 로 상관을 반영한다.

2차 차분 `d_i = (y_i − (y_{i−1}+y_{i+1})/2)/√1.5`, 준-Poisson 비 `κ = median(d²/ys)/0.4549` (χ²₁ 중앙값 0.4549).
c 를 모를 때 κ 는 "1 CPS 당 분산" 추정(= 1/c_eff). 표시용 진단이며 피팅의 척도는 φ.

## 7. 서베이 원소 판별
봉우리(`findPeaks`): 국소 2차 평활(반창 hw ≈ 폭/간격, 2hw+1 점) 후 극대 · 음의 2차미분 > 0.
높이 h = **지형 prominence**(scipy.signal.peak_prominences 와 같은 정의): 양쪽으로 이 봉우리보다 높은 점(또는 ±max(3hw, 4 eV) 창 끝)까지의 구간 최솟값 중 큰 것을 기준선으로.
잡음 σ_h = 1.5·√2·√(κ·max(ys,1)/(2hw+1)) (평활 한 점의 σ ≈ 1.5σ/√N, 두 평활값의 차라 ×√2, κ = §6 의 분산/세기 비). SNR = h/σ_h ≥ 문턱(기본 4), `merge` eV 안은 큰 것만.
(v1.1 은 국소 기준선 대비 높이와 평활 전 σ 를 써서 큰 봉우리 꼬리의 잡음 굴곡을 봉우리로 잡고 SNR 을 약 2 배 과대평가했다.)
원소 e 의 광전자 선 {ℓ}: 위치 = Avantage 원소표(`idpos`) → 없으면 X-ray Data Booklet(단, 진공 기준 기체 원소 H He N O F Ne Cl Ar Kr Xe Rn 은 Booklet 미사용),
가중 w_ℓ = Scofield σ_ℓ / max σ. 같은 궤도 doublet 이 허용폭 tol 보다 가까우면 하나로 합침. Auger 선(Avantage KE, BE = hν − KE)은 위치만.
기본 제외: H, He, 방사성 원소(Tc Pm Po At Rn Fr Ra Ac Pa Np Pu Am Cm Bk Cf Es Fm).

**탐욕적 배정** (가장 큰 미설명 봉우리부터 — 수작업 판별 순서):
1. 미설명 봉우리 중 가장 큰 P. 각 원소의 센 선(w ≥ 0.3) ℓ* 를 P 에 맞추면 δ = P − pos(ℓ*), δ 를 [−4, +6] eV 로 자름.
2. `score = [Σ_{맞은 ℓ} w_ℓ (1 − ½ d_ℓ/tol) · (이미 설명된 봉우리면 0.3) − 0.6 Σ_{빠진 ℓ, w_ℓ>0.2} w_ℓ] · exp(−(δ/5)²) · w_{ℓ_P}`
   (한 봉우리는 한 선에만; w_{ℓ_P} = P 를 설명하는 선의 가중 — 큰 봉우리를 약한 선의 우연 일치로 설명하는 것을 억제).
3. 최고 점수 > 0.5 이면 채택하고 그 원소가 맞춘 봉우리(광전자·Auger)를 '설명됨'으로. 아니면 P 는 '미배정'.
4. 남은 원소는 모든 봉우리 기준 최고 점수로 정렬. 점수 > 0.5 이고 자기만의 봉우리가 있으면 '가능'(겹친 원소일 수 있음, 기본 미체크), 자기 봉우리가 모두 설명됐으면 '겹침'.

검증(test/survey_bench.js, 참값을 아는 합성 survey 150 개: 원소 3–8 개, 화학 이동 −1.5…+3.5 eV, 선 높이 ∝ 농도·σ·KE^−0.4, FWHM 1.6–2.4 eV,
비탄성 계단, Auger, Poisson 잡음, 문턱 4): 주선 SNR > 15 인 원소 채택 97.2 %, 없는 원소 오채택 0.05 개/스펙트럼(정밀도 98.7 %).
(v1.1 의 봉우리 찾기로는 97.7 %·0.59 개·87.6 %.)
비교(개발 중 예비 시험 — 40 개, 참값 기준 농도 > 0.08): 원소별 독립 점수(δ 공유만) 방식은 정밀도 72 %·오채택 1.8 개, 선 세기 일관성 벌점을 더한 방식은 61 %·2.7 개로 더 나빠(겹침 때문에 세기 추정이 흔들림) 버렸다.
'가능' 판정의 실제 적중은 3/6(150 개 중) — 그래서 기본 미체크.

## 8. 원소 라이브러리 (src/xps_lib.js — build/build_lib.py 가 생성, 직접 고치지 않음)
- 핵심 준위 BE: X-ray Data Booklet Table 1-1 (원소 상태, 금속은 Fermi 준위 기준). doublet 간격 = 두 subshell 차. Booklet 인쇄값의 순서 역전(Co 3p, Na 2p)은 `soWarn` 로 표시하고 간격을 쓰지 않음.
- 감도계수: Scofield 1976 (Al/Mg, C 1s = 1) — 디지털 표를 논문 판독값 263 개와 대조(1 %), 빠진 31 개는 논문값으로 채움, 1 개 교정. doublet 합은 항상 두 subshell 합으로 다시 계산(디지털 표에 합이 없거나 Ni 2p·Cl 2p 처럼 한 subshell 값이 들어간 것 포함 347 건).
- ALTHERMO1: Avantage 감도계수(Scofield 와 다름, 예: Ta 4f 9.09 vs 8.62). Avantage 의 transmission·energy compensation 과 짝이 맞아야 Avantage at% 를 재현.
- 빌드는 검사 하나라도 어긋나면 중단(조용한 대체 없음).

## 9. 자동 모형 (원소·궤도와 무관한 규칙)
원소·궤도는 파일 이름에서 읽는다(`parseRegionName`: 주기율표 기호 목록으로 검사, 대소문자 보정; 원소를 못 찾으면 `unknown` — 정보 탭에서 직접 지정).
광원: hν 가 1486.6 ± 3 eV 면 Al Kα, 1253.6 ± 3 eV 면 Mg Kα 의 Scofield 를 쓰고, 그 밖(단색 Ag 등)은 감도계수를 직접 입력하라고 경고한다.

**doublet 정보** `doubletInfo(el, orb)`: p·d·f 궤도의 두 subshell (p 3/2·1/2, d 5/2·3/2, f 7/2·5/2), 넓이 비 r = (2j₂+1)/(2j₁+1) (0.5·2/3·0.75).
Δ 는 Avantage 원소표 두 위치의 차가 있으면 그것, 없으면 Booklet(단 `soWarn` 이면 쓰지 않음). 두 출처가 모두 있으면 차(spread)를 같이 보관.
`dEConstraint`: Δ 는 **항상 범위 안에서 자유** — `[Δ − w, Δ + w]`, `w = max(0.3, 0.05Δ, spread + 0.1)` eV.
(Booklet 값은 원소 상태·인쇄 정밀도의 한계가 있어 고정하면 치우친다: 예 Si 2p Booklet 0.40 vs 실측 ≈ 0.61 eV.) Δ < 1.2 eV 이면 '분해되지 않음'으로 표시만 한다.

**출발 모형** `suggestModel`: 배경을 빼고(active·none 이면 Shirley) `findPeaks`(폭 0.5 eV, 문턱 4, merge 0.9 eV) 상위 8 개 중 가장 큰 봉우리 높이의 2 % 이상인 것을 BE 순으로 훑는다.
봉우리 q 의 짝 = q + Δ ± max(0.7, 0.08Δ) 안의, q 의 1.3 배보다 낮은 봉우리. 짝이 있거나 q + Δ 가 창 안(끝에서 0.5 eV 안쪽)이면 doublet.
연결 기준(anchor)은 **넓이가 가장 큰 성분**: 가장 큰 doublet 의 Δ 는 `dEConstraint`, 나머지 doublet 은 거기에 연결(같은 원소·궤도의 화학 상태는 Δ 가 같다는 가정).
r·q 고정, F 0.3–4 eV 자유, BE 는 창 안 [E_min, E_max] 에서 자유, m 은 가장 큰 성분에 연결. 넓이 출발 = 높이 × 모양의 높이→넓이 환산.
(v1.2 는 가장 낮은 BE 의 doublet 을 기준으로 삼아, 실측 Ba 4d 에서 0.7 % 높이의 작은 봉우리가 기준이 되고 피팅이 그것을 창 밖·넓이 0 으로 보내 퇴화했다.)

**다중항 준위** `multipletInfo(el, orb)`: 열린 껍질과 내각 구멍의 교환 상호작용이 스핀–궤도 갈라짐과 비슷하거나 큰 준위 —
Ti–Cu 의 3s·3p(열린 3d), Ce–Yb 의 4s·4p·4d·5s·5p(열린 4f). 원소·궤도 규칙(산화 상태를 모르므로 Ti⁴⁺ 3d⁰, Cu⁺/Cu⁰ 3d¹⁰, Yb²⁺ 4f¹⁴ 같은 닫힌 껍질 예외는 안내만).
이 준위에서 `suggestModel` 은 doublet 을 만들지 않고 넓은 성분(F 출발 2 eV, 0.3–8 eV)으로 모양(envelope)만 맞춘다 — 성분은 화학 상태가 아니고 넓이 합이 정량값.
**보이지 않는 주선**: 짝이 없는 봉우리 q 에 대해, 앞에 doublet 이 있고, q.be − Δ 자리가 기존 성분 중심·짝 위치와 0.5 eV 넘게 떨어져 있으며, 그 자리의 신호에서 기존 성분(짝 포함)의 예측을 뺀 나머지가 q 의 높이 이상이면 주선을 q.be − Δ 에 둔다(높이 = 그 나머지)
(주선이 짝보다 크다: 높이 비 1/r > 1. 예: Au 4f 두 상태 84.0·86.5 eV — 86.5 와 앞 doublet 의 4f5/2(87.67)가 한 덩어리라 90.17 만 따로 보임).

**잔차 보강** `refineByResidual` (판단은 항상 `active` Shirley + 기울기 s 자유 배경으로 — 정적 배경으로 판단하면 배경 모양 오차가 잔차 덩어리로 남아
참 1 성분 합성 스펙트럼 143 개 중 106 개에 가짜 성분이 더해졌다):
1. 맞춘 뒤 표준화 잔차를 ±0.6 eV 창에서 합하고 √(창 점수)로 나눈 z, 다시 AR(1) 팽창 √((1+ρ₁)/(1−ρ₁)) 로 나눈다(ρ₁ ∈ [0, 0.95]).
2. 한 걸음마다 후보 수를 모두 맞춰 보고 **BIC 를 가장 많이 줄이는 하나를, 10 이상 줄 때만** 받는다(없으면 멈춤):
   (가) 창 가장자리를 뺀 최대 z ≥ 4 인 BE 에 성분 추가(모양·m 은 가장 큰 성분, doublet 이면 가장 큰 doublet 에 Δ·r·q 연결, 최대 4 개);
   (나) doublet 넓이 비 r 풀기 [0.75, 1.33]×통계 값; (다) 폭 배율 q 풀기 [0.6, 1.6] — (나)(다)는 모든 화학 상태에 같은 값(가장 큰 doublet 에 연결), 분해되는 doublet(Δ ≥ 1.2 eV)만.
   통계 비 (2j+1) 는 얕은 준위에서 어긋나는 일이 흔하다(실측 Ba 4d: r = 0.69, ΔBIC 84). 비를 성분 추가와 경쟁시키는 까닭: 비를 풀지 않고 성분만 더하면 짝 자리 잔차를 가짜 성분이 메우고(합성 양성 대조: 성분 4 개),
   비를 늘 먼저 풀면 빠진 화학 상태를 비가 흡수할 수 있다 — 그래서 둘 중 데이터를 더 잘 설명하는 쪽을 고른다.
3. **죽은 성분 덜기**: 처음과 (더한 것이 있으면) 끝에, 넓이 몫 < 5 % 이거나 중심이 창 밖인 성분을 작은 것부터 빼 보고 BIC 가 2 넘게 나빠지지 않으면 뺀다.
   그 성분이 연결 기준이었으면 `removeComp` 가 남은 것 중 가장 큰 성분을 새 기준으로 삼고 뺀 성분의 제약(예: Δ 자유 범위)을 물려준다.
   (겹친 다중항 — 예: Au 4f5/2 와 다른 화학 상태의 4f7/2 — 의 위쪽 짝을 출발 모형이 따로 잡으면 피팅이 그 성분을 넓이 0 으로 보내 퇴화를 만든다.)
4. 성분 구성을 가지고 사용자 배경으로 다시 맞춘 값을 출발값으로 돌려준다. 화면은 이 결과를 성분 표에 넣을 뿐 피팅 결과로 취급하지 않는다(사용자가 '피팅').

검증(test/all_elements.js — Z 3–92 에서 제외 원소를 뺀 모든 원소의 Scofield 기준 센 궤도 2 개, BE 10–1300 eV, 143 경우;
참: GLS F 1.2·m 0.3, Δ_true = 라이브러리 ± 0.1 eV, active Shirley + 기울기 배경, Poisson 잡음 c = 0.3):
이름 해석 5 형식 × 143 모두 · Scofield Al 모두(Mg 빈칸 15) · Δ 없음 3(Ho·Er·Tm 4d) ·
단일 상태 자동 제안→피팅 BE ±0.06 eV·넓이 z < 4.5(active)/±8 %(smart)·Δ ±0.15 eV 143/143 · 잔차 보강이 가짜 성분을 더함 0/143 ·
화학 상태 2 개(δ = 2.5 eV, 넓이 2:1)를 자동 제안 + 잔차 보강만으로 성분 수·BE ±0.15 eV·넓이 ±10 % 143/143 ·
참이 통계 비인데 r·q 를 푼 경우 0/143 · 양성 대조(참 r = 1.2 × 통계 값, 두 상태, 분해되는 doublet 105 경우) r 을 풀어 ±5 % 회복·성분 수 2: 104/105(Au 4f: 4f5/2 가 다른 상태 4f7/2 와 1.2 eV 차로 겹쳐 r·q 가 함께 움직임).
정적(smart) 배경의 넓이 오차는 창 폭에 따라 커진다: ≤ 20 eV 최대 3 %, 20–40 eV 최대 8 %, > 40 eV 최대 14 % (화면 경고는 20 eV 초과).

## 10. 화학량론 비교 · Auger parameter
**화학식** `parseFormula`: 원소 기호(주기율표 목록으로 검사, 두 글자 먼저)·숫자(소수 가능)·괄호 `()`/`[]`(중첩 가능, 뒤 숫자를 안 전체에 곱함). 빈 값·모르는 기호·괄호 짝·0 개는 오류.
**비교** `stoichCompare(rows, ν, ref)`: rows 의 nᵢ = Aᵢ/(Sᵢ·Tᵢ·λᵢ)(정량 §4 와 같음). 화학식에 있고 측정된 원소 F 끼리만:
xᵢ = nᵢ/W, W = Σ_{j∈F} nⱼ; 기대 eᵢ = νᵢ/Σ_{j∈F} νⱼ; 측정/화학식 = xᵢ/eᵢ.
σ(xᵢ)² = Σⱼ ((δᵢⱼW − nᵢ)/W²)² σ(nⱼ)² (넓이끼리 독립). 기준 원소 비 (nᵢ/n_ref)/(νᵢ/ν_ref), 상대 σ = √((σᵢ/nᵢ)² + (σ_ref/n_ref)²).
화학식 밖 원소(C·O 오염 등)와 측정되지 않은 화학식 원소는 따로 보고하고 정규화에서 뺀다. 같은 원소 여러 행이면 첫 행만 쓰고 알림. 감도계수·T·λ 계통 오차는 포함하지 않음.
깊이: 화면은 TPP-2M 이면 3λ(nm), KE^n 이면 λᵢ/λ_ref, 선들의 KE 비가 1.8 넘으면 '보는 깊이가 다름' 경고.
검증(test §8): 해석 7 식·오류 7 식, 비교 값 정확(1e-12), σ(델타) vs Monte Carlo 20000 회 < 3 %(실측 0.9 %).

**Auger parameter** `augerParam`: α′ = BE_p + KE_A = BE_p + hν − BE_A,겉보기 (BE 눈금의 Auger 봉우리 위치). 대전 이동 δ 는 두 항에서 지워져 α′ 는 Δ 와 무관.
화면 좌표(Wagner plot): BE_corr = BE_p + Δ, KE_corr = hν − (BE_A + Δ). hν 는 Auger 스펙트럼의 것.
위치: 피팅된 성분의 BE(σ = 공분산 σ_BE, α′ 의 σ = √(σ_p² + σ_A²)) 또는 구간 최댓값 `peakMaxPos`(국소 2차 평활 반창 ≈ 0.25 eV, 가장 큰 점을 세 점 포물선으로; 최댓값이 구간 끝이면 null — σ 없음).
survey 기본 구간: 라이브러리(Avantage) Auger KE 중 가장 센 선의 hν − KE − Δ ± 6 eV. 문헌 α′ 표는 싣지 않는다(검증된 자료 없음) — 기준점은 사용자가 출처와 함께 넣음.
검증(test §8): 합성 Na 1s/Na KLL 에서 δ = 0·1.73·4.37·−0.6 eV 모두 α′ 참값(잡음 없음 < 5 meV; Poisson 잡음 100 회 평균이 4 SE 안), 구간 끝 최댓값 → null.

## 11. 모형 민감도 (분석 선택에 따른 흔들림)
목적: 성분 구성(성분 수·모양·제약)은 그대로 두고, 분석자가 합리적으로 고를 수 있는 선택을 바꿔 다시 맞춰 관심 값의 흔들림을 잰다 — '피팅의 주관'을 숫자로.
**요인**: 배경 종류(Smart · Shirley · 피팅 안 Shirley + 기울기 · 창 ≥ 30 eV 이면 Tougaard; 사용자가 고름), 끝점 평균 점 수(지금 정적 배경 종류에서만 3·10), 낮은/높은 BE 창 끝 ±dw(기본 0.5 eV).
v1.7: 지금 배경이 `spline` 이고 사용자 매듭이 있으면 **매듭 위치** 요인 `kn` — 매듭 하나씩 ±dk(기본 0.5 eV) 옮긴 변형(자동 높이는 새 자리에서 다시, 끌어 둔 높이는 그대로; 다른 요인은 기준 수준, 두 방식 모두 하나씩 → 2 × 매듭 수).
주효과 kn = {기준, 각 매듭 변형} 수준별 값의 범위.
**방식**: 빠르게 = 0 번(지금 설정) + 요인 하나씩(Smart 기준 9 개), 꼼꼼히 = 모든 조합(45 개). 창이 스펙트럼 끝을 넘으면 끝으로 자르고 같은 변형은 한 번만.
**다시 맞춤**: 각 변형은 지금 피팅된 성분 값에서 출발(넓이는 NNLS 로 다시), 성분 중심이 창 밖이면 건너뜀. 0 번은 직접 피팅과 같은 최솟값(차이 < 1e-6 σ).
**요약**: 성분·값(BE·넓이 %·넓이·FWHM)마다 기준·최소·최대·평균·sd(변형 같은 무게 — 창·배경이 달라 χ²·BIC 로 비교 불가).
주효과 = 다른 요인이 기준 수준인 변형들의 수준별 평균의 범위(배경 효과는 끝점 점 수를 따지지 않음).
**보고 구간(가능 범위)** = [최소 − 2σ통계, 최대 + 2σ통계] ∩ [0, 100]; σ통계 = max(공분산 σ × AR(1) 팽창 √((1+ρ₁)/(1−ρ₁)), Monte Carlo sd). 반범위 > 2σ통계 이면 '분석 선택이 지배'.
지금 배경과 같은 종류의 변형은 사용자의 배경 제약(예: 기울기 고정)을 그대로 쓴다(0 번 = 사용자 모형). 지금 창 밖에 있던 성분은 판정하지 않고, 지금 창 안의 성분 중심이 변형 창에서 잘리면 그 변형만 건너뜀.
보장 범위: 참 배경이 후보 안에 있을 때만(합성 시험의 참 배경 = 후보의 하나).
**검증**(test §9; 겹친 두 성분 2:1, 참 배경 = 피팅 안 Shirley + 기울기, Poisson 잡음): 기준 배경을 Smart 로(틀리게) 고르면 참 넓이 % 를 통계 ±2σ 는 5/40 만, 가능 범위는 37/40 포함.
기준 배경이 맞으면 통계 ±2σ 19/20, 가능 범위 20/20. (sd 로 합친 √(σ²+sd²) 는 20/40 — 모형 오차는 변형들의 sd 가 아니라 범위로 보고해야 한다.)

## 12. 측정 선모양 · 기준 스펙트럼 선형 결합 (LCF) (v1.7)
목적: CasaXPS 'Line-shapes defined from spectra'·Avantage 'NLLS fitting to reference spectra' 와 같은 일. 두 프로그램의 코드·자료는 쓰지 않았고 아래 식은 여기서 정의한다
(선형 결합 피팅의 일반형은 XAS 의 LCF — B. Ravel, M. Newville, *J. Synchrotron Rad.* 12, 537 (2005) — 와 같다).

### 12.1 template (기준 스펙트럼 → 넓이 1 모양)
입력: 기준 스펙트럼의 피팅 창 (E_i, y_i) (측정 BE, 오름차순), 그 시료의 대전 보정 Δ, 배경 B(그 스펙트럼의 정적 배경 §2; `active` 는 피팅 결과의 배경, 없으면 smart; `none` 은 0).
`z_i = y_i − B_i`, `x_i = E_i + Δ − eref`, `R(x) = interp(x; x_i, z_i) / I`, `I = ∫ interp(x; x_i, z_i) dx` (보간 함수의 정확한 적분, 12.2).
eref (보정 BE) = z 를 국소 2차 평활(§2 smart 와 같은 창)한 최댓값 점 + 세 점 포물선 보정. 마디 밖(x < x_0, x > x_{N−1})은 R = 0. I ≤ 0 이면 오류.

### 12.2 보간
- `pchip`(기본): 3차 Hermite, 마디 기울기 m_k = scipy `PchipInterpolator` 와 같은 정의(Fritsch–Butland):
  d_k = (z_{k+1}−z_k)/h_k; 안쪽 마디는 d_{k−1}·d_k ≤ 0 이거나 하나가 0 이면 m_k = 0, 아니면 `(w₁+w₂)/m_k = w₁/d_{k−1} + w₂/d_k`, w₁ = 2h_k + h_{k−1}, w₂ = h_k + 2h_{k−1};
  끝 마디 `m_0 = ((2h_0+h_1)d_0 − h_0 d_1)/(h_0+h_1)`, 부호가 d_0 와 다르면 0, d_0·d_1 < 0 이고 |m_0| > 3|d_0| 이면 3d_0 (끝도 같은 식).
  조각 적분 `∫ = h(z_k+z_{k+1})/2 + h²(m_k − m_{k+1})/12`. 기준: F. N. Fritsch, R. E. Carlson, *SIAM J. Numer. Anal.* 17, 238 (1980); F. N. Fritsch, J. Butland, *SIAM J. Sci. Stat. Comput.* 5, 300 (1984).
- `linear`: 선형, 적분은 사다리꼴.
- **보간 오차**(해석 모양을 간격 s 로 표본 → 보간 vs 참, 최대오차/최댓값; test/ref/ref_casa.py, JS 가 같은 값을 재현):

  | s (eV) | GLS F=1 m=0.3: PCHIP | 선형 | Voigt wG 0.8 wL 0.3: PCHIP | 선형 |
  |---|---|---|---|---|
  | 0.025 | 1.4e-4 | 4.8e-4 | 1.4e-4 | 4.8e-4 |
  | 0.05 | 5.6e-4 | 1.9e-3 | 5.7e-4 | 1.9e-3 |
  | 0.1 | 2.2e-3 | 7.3e-3 | 2.2e-3 | 7.5e-3 |
  | 0.2 | 8.0e-3 | 2.6e-2 | 8.4e-3 | 2.7e-2 |

  PCHIP(모양 보존, 넘침 없음)은 O(s²)·선형의 약 1/3. 측정 간격 0.05 eV·FWHM ≈ 1 eV 에서 잡음보다 훨씬 작다.

### 12.3 이동 · 넓힘
성분 세기 `A · R(E − be + 0) ⊗ G(·; wG)` — be 를 고정하면 이동 없음(LCF), 자유면 이동 δ = be + Δ_대상 − eref. wG = 0 이면 보간값을 바로 쓴다(be 에 대해 C¹).
wG > 0: 상대좌표 균일 격자 h = 0.04·2^(−k) ≤ min(wG/16, 기준 점 간격/4) (하한 0.0005; 양자화 이유는 §1.8), 범위 [x_0 − 5σ − h, x_{N−1} + 5σ + h],
격자에서 보간값 ⊗ ±5σ 절단·이산 정규화 Gaussian, 데이터 점으로는 Catmull–Rom(§1.8). 검증: quad 대조 최대오차/최댓값 ≤ 4e-5 (PCHIP·선형, wG 0.3·0.8), 넓힘 뒤 ∫ = 1 ± 1e-6.

### 12.4 LCF · 분율 · 공선성
모형 `μ = B + Σ_i a_i R_i(E − δ_i) ⊗ G(w_i)` — 기준마다 MEAS 성분 하나(화면 'LCF: add all' = be·wG 고정). a_i ≥ 0 은 넓이 하한(§3.2), 출발은 NNLS(§3.3), 그다음 LM(또는 Simplex)·IRLS 2회.
R_i 의 넓이가 1 이므로 a_i = 대상에서 그 기준이 차지하는 넓이(CPS·eV), **분율** `f_i = a_i / Σ a_j` (= 결과표 Area %; 공분산 σ 는 §3.4 델타).
**공선성**: X_j = R_j(지금 δ·w)/√v (데이터 점), (a) 척도 조건수 κ = √(λ_max/λ_min), λ = 열을 단위 길이로 맞춘 XᵀX 의 고윳값 (D. A. Belsley, E. Kuh, R. E. Welsch, *Regression Diagnostics* (Wiley, 1980));
(b) VIF_j = (R⁻¹)_jj, R = 평균을 뺀 열의 상관행렬. κ > 30 또는 VIF > 10 이면 경고(분율이 잘 정해지지 않음). numpy `linalg.cond`·`corrcoef` 대조 1e-8.

### 12.5 불확도 — 기준 잡음
공분산 σ(§3.4)는 대상 잡음만이다. Monte Carlo(§3.4)는 회마다 기준 원자료에 잡음을 더해 template 을 **다시 만든다**(배경 빼기·정규화 포함, eref 는 고정):
`y'_i = y_i + √(κ·max(ŷ_i, 10⁻⁶·max ŷ))·z_i` (ŷ = §2 평활, κ = §6 의 2차 차분 분산/세기 비; 일정 가중 자료는 §3.1 const 분산). `opt.refNoise = false` 면 끈다.
한계: 기준 잡음의 이웃 상관은 넣지 않는다(§6: 상관이 있으면 κ 가 작게 나와 구간이 좁아진다).

### 12.6 검증 (test/test_xps_casa.js ↔ test/ref/ref_casa.py)
PCHIP 값·기울기·적분 = scipy (≤ 1e-12) · template x·v·넓이 = 기준 (≤ 1e-12) · 선형 LCF 넓이 = scipy nnls (2회 IRLS, 5.5e-13) · 넓이 σ (0.000 %) ·
이동·넓힘 LCF = least_squares(TRF, 미세 격자 합성곱) |Δ|/σ ≤ 3e-4, σ 0.15 % ·
합성(참값, 기준 봉우리 ≈ 2.6k count·대상 30k count, 60 회): 분율 95 % MC 구간 적중 — 기준 잡음 포함 59/60, 기준 잡음 뺀 MC 47/60, 공분산 ±1.96σ 46/60 → **기준 잡음을 넣어야 구간이 맞다**;
측정 선모양(비대칭 LA⊗G, 두 상태 Δ 1.5 eV·넓이 비 0.5, 둘째 상태 추가 넓힘 0.6 eV, 30 회): ΔBE·비·w 95 % 적중 30/30 · 27/30 · 27/30 (w 평균 0.600 eV).

## 13. CasaXPS 의 선모양·배경 정의 조사 — 문헌 + CasaXPS 2.3.27 실측
**CasaXPS 문서의 식**(N. Fairley et al., *J. Vac. Sci. Technol. A* 41(1) (2023), doi:10.1116/6.0002196, Supplementary Information, Eq. 2–9; Casa Software tutorials 'Test Peak Model'):
l(x) = 1/(1+4x²), g(x; w) = exp(−4ln2 x²) (폭 w 로 척도), 일반화 Lorentzian lg(x; α, β) = l(x)^α (x ≤ 0), l(x)^β (x > 0), `LA(x; α, β, n) = N ∫ lg(τ; α, β) g(x − τ; f_G(n)) dτ`,
`LA(m) = LA(x; 1, 1, 1401 − (m/100)·1401)`, Shirley 성분 SB(0)F(x) = ∫_x^∞ F, Tougaard `T(E) = ∫_E^∞ F(E′ − E) J(E′) dE′`, `F(x) = U(x; B, C, D, T0) = Bx/((C − x²)² + Dx²)` (x > T0, 아니면 0),
B 는 봉우리보다 낮은 KE 의 세기에 맞춰 척도. 표기 'U 2 T[B, C]', 'U 3 T[C, D, T0]'(B 자동).
→ §1.6 LA 와 같은 꼴(이 모듈의 x > 0 = 높은 BE = CasaXPS 의 x ≤ 0 쪽(KE 축) — 실측으로 꼬리 쪽 확인). **f_G(n)(정수 n → Gaussian 폭)은 문서에 없다**; LA(1.53,2.43,0) 한 판독(n = 0)으로는 정할 수 없다(아래 추가 판독 목록).
§2 tougaard 는 같은 손실 핵(U2·U3)을 쓰고 v1.7 에서 T0 를 더했다; 차이는 척도 규칙(여기: y_lo 를 빼고 B(E_hi) = y_hi; CasaXPS: J 그대로, 낮은 KE 의 세기에 맞춤) — CasaXPS 판독 대기. 앞서의 'U 2 Tougaard 2866,1643,1,0' 판독은 D = 1 이 들어 있어 순수 U2 대조가 아니다(쓰지 않음). (v1.7, 2026-10-28)
CasaXPS 기능 목록(Casa Software, 'CasaXPS Feature List', rbdinstruments.com/phi/files/casa-xps/casaxps-feature-list.pdf)에 있는 이름:
배경 'Linear, Shirley, **Analytic Shirley**, **Kappa parameters Shirley**, Two/three parameter Tougaard, **Cubic Spline**', 선모양 '… Line-shapes defined from spectra', 최적화 '**Marquardt, Simplex**'.
식을 확인할 수 있었던 것과 없었던 것:
- **LA**: 공개된 정의 — l(x) = 1/(1+4x²) (x 는 FWHM 단위), g(x;w) = exp(−4 ln2 (x/w)²), Voigt LA(1,n) = N ∫ l(τ) g(x−τ; w(n)) dτ
  (B. Moeini, M. R. Linford, N. Fairley et al., "Definition of a new (Doniach–Sunjic–Shirley) peak shape …", *Surf. Interface Anal.* 54, 67 (2022), Eq. 2–4; OSTI 1976397).
  비대칭 LA(α,β,m): Lorentzian 을 양쪽에서 각각 α, β 제곱 — "a and b describe the spread of the tail to **higher and lower BE**" (M. J. Acres et al., "Core level photoemission line shape selection: atomic adsorbates on iron", Imperial College Spiral).
  → 여기 §1.6 `LA` 와 **같은 형태**(α = 높은 BE 쪽, β = 낮은 BE 쪽, Lorentzian 거듭제곱 ⊗ Gaussian)임을 확인. 차이: CasaXPS 의 세 번째 인자 m(정수)이 Gaussian 폭 w 로 바뀌는 식 w(n) 은
  공개 문서에 없다(위 논문도 'w = w(n)' 으로만 씀) → 여기서는 wG 를 eV 로 직접 둔다. 넓이 정규화는 ℝ 해석식(§1.6).
- **LF**(finite Lorentzian, LF(α,β,w,m)): 말로만 공개 — "w is a damping parameter to limit the range over which the asymmetry extends; m specifies the width of the Gaussian" (Acres et al.);
  "a parameter that suppresses the peak tail" (*J. Vac. Sci. Technol. A* 39, 013204 (2021)). **감쇠 함수의 식은 찾지 못함** — CasaXPS 매뉴얼/Cookbook PDF(casaxps.com)는 이 환경에서 읽을 수 없었다
  (WebFetch 가 PDF 를 이진으로만 받고, 직접 내려받기는 정책상 막힘). R. Lindsay 등 'Advanced line shapes in XPS II: The finite Lorentzian (LF) line shape' (*Vac. Technol. Coat.* 2020) 도 열 수 없음. **구현하지 않음** — 추측한 식으로 'CasaXPS LF' 라고 부르지 않는다.
- **Kappa parameters Shirley**: 이름만 확인, 정의 없음 → **구현하지 않음**.
- **Analytic Shirley**: 이름만 확인, 정의 없음 → **구현하지 않음** (반복 Shirley §2 는 그대로).
- **SGL(p)** (CasaXPS 의 GL 합): 실측 — 이 모듈 GLS(m 0.3, 넓이 혼합) 자료를 CasaXPS SGL(30) 으로 맞추면 RMS 882.7(GL(30) 은 8.3), 넓이 15279.4/15000, fwhm 1.177 → **혼합 정의가 다르다**.
  가설(미확인): SGL(p) = (1−p/100)·exp(−4ln2 x²/F²) + (p/100)/(1+4x²/F²) (높이 기준 혼합) — 그러면 넓이 혼합 m ↔ 높이 혼합 p/100 = m·L₀/(m·L₀ + (1−m)·G₀), L₀ = 2/π, G₀ = 2√(ln2/π) (m 0.3 → SGL(19.5)).
  확인 방법은 보고서의 CasaXPS 대조 목록. 확인되기 전에는 'CasaXPS SGL 과 같은 혼합' 이라고 하지 않는다(화면: 'CasaXPS SGL: different mixing').
- **GL(p)**: 실측 — 이 모듈 `GLP`(m = p/100) 와 같은 모양(넓이 14999.8/15000), CasaXPS 의 'fwhm' = **실제 FWHM**(1.15343 vs 여기 실제 FWHM 1.15363; F = 1.2). 화면 성분 카드에 'CasaXPS: GL(p), fwhm' 표시.
- **LA(α,β,m) 폭 숫자**: 실측 — 이 모듈 LA(F 1.0, α 1.53, β 2.43, wG 0) 자료를 CasaXPS LA(1.53,2.43,0) 로 맞추면 넓이 14997.0(여기 14998.3), 같은 꼬리 쪽 — 모양은 같다. 그런데 CasaXPS 의 폭 숫자는 0.580056
  (여기 F = 1.0, 실제 FWHM 0.66569). 단순 후보(F, 실제 FWHM, 한쪽 반폭 ×2, F√(2^{1/p}−1) 의 p = α, β, (α+β)/2, √(ln2/α) 꼴)는 모두 안 맞음 → **환산 미확인**. 추가 측정용 자료(대칭 α=β=1·2, F=2) 는 보고서에.
- **U 2 Tougaard "2866,1643,1,0"** 판독(C1s 넓이 37823.6): 위 Eq. 9 로 보면 D = 1 이 섞인 설정이라 순수 U2 대조가 아님 — 대조에서 뺌. 시험했던 4 가지 정의(34366–38898)는 기록만.
- **Cubic Spline** 배경: 식이 공개되지 않아 §2 `spline` 을 여기서 정의(자연 3차 spline). CasaXPS 와 같은 결과인지는 **CasaXPS 대조 대기**.
식이 들어 있는 CasaXPS 문서(매뉴얼 'Line shapes' 장 · Cookbook)의 해당 쪽을 받으면 LF·Kappa·Analytic Shirley 를 그 식대로 넣고 기준 구현을 더한다.

**CasaXPS Manual 2.3.15 rev 1.3 (Casa Software, 2009) p. 55·301** (식 그림):
`L(x: f, e, m) = 1/(1 + 4(m/100)((x−e)/f)²)`, `G(x: e, f, m) = exp(−4ln2(1 − m/100)((x−e)/f)²)`, `GL(m) = G·L` (→ §1.3 GLP 와 같음),
`LA(x: α, β, f, e) = L(x: f, e, 100)^α (x ≤ e), ^β (x > e)` (x = KE → α 는 높은 BE 쪽, f = Lorentzian FWHM = 여기 F; 2.3.15 에서 m 은 0–499 정수, 뒤 판은 1401),
`T(x: s, k, e, f) = s·exp(−k(x−e)/f)` (x ≤ e), 1 (그 밖), `GL(m)T(k) = GL + (1 − GL)·T(x: 1, k, e, f)`,
`SGL(m) = (1 − m/100)·G(x: e, f, 1) + (m/100)·L(x: e, f, m)` (인쇄된 그대로 — 안쪽 매개변수 1, m 은 오식일 수 있음).
- **LF**: "α 와 β 가 폭 매개변수 w 로 정해지는 매끄러운 함수로 일정한 값까지 커진다" 는 말뿐, 식 없음 → **구현하지 않음**.
- **GL(m)T(k)**: 인쇄된 지수 −k(x−e)/f 는 x ≤ e 에서 커지므로 꼬리는 `exp(−k|x−e|/f)` 로 읽었다(가정). §1.10 `AGT` 에 `tu = "fwhm"`(꼬리 길이 단위 = f), TM = 0, ET = k 로 같은 꼴 — CasaXPS 판독(V2_GL30T2) 대기.
- **SGL**: 인쇄형(V2_SGLprinted30)과 높이 기준 합(V2_SGLh30) 두 가설의 판독 자료를 만들었다. 프로젝트 쪽 대조: 인쇄형을 GLS 자료에 맞춰도 CasaXPS 2.3.27 판독(15279.4/1.17737)이 나오지 않음 → 미확인.
- **LA 폭 숫자**: 매뉴얼에서 f = Lorentzian FWHM 이므로 CasaXPS 의 LA 'fwhm' 칸 0.580056(여기 F = 1.0)은 f 가 아니라 다른 양으로 보인다
  (낮은 BE 쪽 반폭 ×2 = 0.5745(정확)·0.5779(0.05 eV 격자) 와 가깝지만 같지 않음) — V2_LA_* 판독 대기.

## 13b. Avantage 정의 대조 (Thermo Scientific Avantage Data System manual — 정의만 사용)
- **GL 곱(Product)**: Avantage `exp(−(1−M)ln2·(x−x₀)²/β²)/(1 + M(x−x₀)²/β²)`, β ≈ ½FWHM, M = L/G Mix/100 → §1.3 `GLP` 와 **같은 함수**(F = 2β). Avantage 피크 표의 FWHM 은 측정한 실제 FWHM(여기 결과표와 같음).
  (설명서 본문 식은 그림이라 위 꼴은 Sherwood(Briggs & Seah App. 3) 표준형과 설명서의 'β 는 거의 ½FWHM, M = 1 Lorentzian' 서술로 확인.)
- **Convolution(Voigt)**: Avantage 는 (fV, M) 로 보이고 fG·fL 은 M 에서 정한 비 fL/fG 와 Olivero–Longbothum 식 fV ≈ 0.5346fL + √(0.2166fL² + fG²) 로 푼다; 'fit param' FWHM = fG, M = 100 % 면 fL = fV.
  비 fL/fG(M) 의 식은 설명서에서 그림이라 읽지 못함 → 화면 표시는 **M = fL/(fL+fG) 가정**(`avantageVoigt`·`voigtFromAvantage`, 왕복 검증) — Avantage 대조 대기. 모양은 §1.4 V 와 같은 진짜 Voigt.
- **Sum**: 정의 식이 설명서 텍스트에 없음 → §1.2 GLS 와 같은지 **미확인**(이름만 'Avantage Sum' 이라 쓴 것은 v1.7 에서 '대조 대기'로 고침).
- **꼬리 함수**: §1.10 `AGT` 로 구현(GLT §1.7 는 다른 정의 — 넓이 1 지수 꼬리 합성곱).
- **Shirley·Smart**: §2c. **투과 함수**: §4 'Avantage 형태'. **ECF**: Scofield 라이브러리 KE^0.6(기본), Wagner KE^1.0 또는 TPP-2M 평균 — §4 의 `power` n = 0.6·`tpp2m` 과 같은 구조.
- **정량**: Normalised Area = A/(SF·TXFN·ECF), at% = 100·N_i/ΣN — §4 와 같음.

## 14. Simplex (Nelder–Mead) (v1.7)
변수: MINUIT 변환 u (§3.2) 를 매개변수마다 척도 s_j 로 나눈 z = (u − u₀)/s_j. s_j = |u(p₀ + dp_j) − u(p₀)|, dp = 0.1 eV (BE·Δ), 0.1|p| (그 밖), |p| 가 0 이면 0.05(hi − lo) 또는 10⁻³·max|y|.
첫 simplex: z₀ = 0, z₀ + e_j. 계수 반사 ρ = 1, 확장 χ = 2, 수축 γ = ½, 축소 σ = ½, 바깥/안쪽 수축 구분 — J. A. Nelder, R. Mead, *Comput. J.* 7, 308 (1965); J. C. Lagarias et al., *SIAM J. Optim.* 9, 112 (1998) (scipy 의 표준 계수와 같음).
정지: max_i‖z_i − z_0‖_∞ ≤ 10⁻⁷ 그리고 max_i|f_i − f_0| ≤ 10⁻¹²·max(1,|f_0|). 수렴하면 최적점에서 새 simplex(걸음 0.2)로 다시 시작, 개선이 10⁻¹¹ 상대 이하가 될 때까지(최대 6 회). 평가 한도 max(4000, 1500(p+1)).
검증: 같은 분산 v 의 GLS 2 성분 + 직선 배경에서 JS Simplex = least_squares 최솟값 |Δ|/σ 9.6e-7, χ² 상대 < 1e-9, = scipy Nelder–Mead(경계·재시작) 1.7e-7σ, = JS LM 3.4e-7σ;
기준 구현 ① 의 피팅 문제 8 개(GLS·GLP doublet·Voigt·DS·LA·GLT, 배경 5 종) 모두 Simplex = LM (|Δ|/σ ≤ 3e-3, χ² 상대 ≤ 1.2e-7; 평가 800–2700 회, 0.4–5 s).

## 15. 여러 스펙트럼에 같은 모형 (v1.7) — Avantage/CasaXPS 의 모형 전파와 같은 목적
`propagateModel(model, Δ_src, Δ_dst, [E_min, E_max])`: d = Δ_src − Δ_dst. 성분 BE 와 free BE 한계 +d (§5 shiftModelBE), spline 매듭 e +d, 창 [lo+d, hi+d] 를 대상 스펙트럼 끝에서 자름,
배경 종류·제약·가중·최적화는 그대로, active 의 c₀·k 는 다시 초기화. 넓이는 대상에서 NNLS 로 다시 출발(§3.3). 대상마다 Worker 하나(`batchJob` = makeProblem → fitOnce → finalize).
결과표(화면 'Batch results', CSV): 스펙트럼마다 χ²ν, 성분(이름으로 맞춤) 보정 BE·FWHM·넓이·넓이 % 와 σ, 그 원소의 at %(그 시료의 정량 §4 설정).
검증(합성 깊이 분포 6 개, 금속 분율 0.9→0.1, 시료마다 다른 대전 이동 −0.8…+1.2 eV, 참 배경 = 피팅 안 Shirley): 보정 BE 95 % 적중 11/12, 분율 5/6.

## 16. 화면 (v1.8) — Series · ARXPS · Imaging · Process
수식은 SPEC_multi.md · SPEC_proc.md. 여기는 화면이 엔진을 어떻게 부르는지만.
- **작업 공간**: 위쪽 전환 Spectra(기존) · Series · ARXPS · Imaging. 설정은 `S.ser`·`S.arx`·`S.imu`·`S.prc`(IndexedDB "cur"), 영상 stack 은 IndexedDB "img"(구조 복제), 계산 결과는 메모리(`RES`, `PRES`) — 다시 열면 다시 계산.
- **Worker**: `xps-core-src` 뒤에 `xps-multi-src` · `xps-arxps-src` · `xps-image-src` · `xps-proc-src` 글을 이어 붙인 하나의 Blob(`inject_all.py` 가 이 순서로 주입). 명령 `xj` = `xjobs()`(UI 파일 안, 글째로 Worker 에 들어감): PCA+MLPCA 잡음 시험, MCR-ALS, 기준 NNLS, MaxEnt, 영상 넓이·PCA·재구성·NNLS(화소 묶음을 여러 Worker 에), 디컨볼루션. 진행 막대(비율을 모르면 흐르는 막대 + 경과 시간) · Stop = Worker 종료. 주입되지 않은 환경에서는 주 스레드.
- **불러오기**: VAMAS MAP/SEM 블록(x, y 있음) → 목록이 아니라 Imaging stack. 여러 열 텍스트: 모든 열 이름에 서로 다른 숫자 → series 스펙트럼(`meta.mcol`), Imaging 화면에서는 영상 stack(열 = 화소, 행 우선; nx 를 바꿀 수 있음). 열 이름에 숫자가 없는 여러 열(Avantage 내보내기)은 전처럼 2 열 스펙트럼. 중복 판정에 VAMAS 블록 번호·열 번호를 더함(같은 첫 세기의 다른 단계를 버리지 않게).
- **Series**: 데이터 = 시료(같은 region 2 개 이상) 또는 목록에서 체크한 스펙트럼. 좌표 = VAMAS 실험 변수(자동) · take-off 각 · 열 이름 숫자 · 불러온 순서 · 직접 입력. 모든 계산은 시료마다의 대전 보정 BE 격자(`commonGrid`). region 마다 Integrate(`integrateSeries`) 또는 Fit model(단계마다의 피팅 성분 넓이; ‘Apply & fit all’ = Peak Fit 의 일괄 피팅). `compositionSeries` 는 정량 탭의 S 표·투과·IMFP 를 그대로. 정보 깊이 보정은 깊이 눈금(sputter rate)이 있을 때만, λe 기본값 = TPP-2M(정량 탭 물질)/10 × cos(take-off).
  다변량: 잡음 시험(MLPCA)은 늘 원자료 창으로(Poisson), TFA·MCR·기준 NNLS 는 배경을 뺀 자료(선택)로. n 기본값 = 잡음 한계 요인 수(준위마다 다른 배경은 요인 하나 — 뺐으면 1 줄임은 사용자가). 기준 NNLS 의 σ: 무가중이면 √φ, Poisson 이면 √max(φ,1) 를 곱함(엔진의 sFraction 도 같게).
- **기준 후보**(TFA · 기준 NNLS · ARXPS · 영상 NNLS 공용): 저장한 기준(template, `tplShape`) + 목록에서 선택한 스펙트럼의 피팅 성분(`compEval`, 그 시료의 Δ 로 보정 BE) + (영상) ROI 기준. 쓰는 격자에서 넓이 1 로 맞춤.
- **ARXPS**: 각도는 법선에서(‘from surface’ 면 90° − 값). 세기 = 기준 NNLS(배경을 빼면 무가중, σ×√φ) · 창 넓이 · 피팅 성분. 덮개층 `fitOverlayer`(σ×√max(φ,1)) + Hill 그림(L_o = L_s 일 때만 정확), 감쇠 길이 계산기(Seah–Dench · TPP-2M×(1−0.735ω)), 다층 `fitMultilayer`(종마다 L 하나를 모든 층에), RDP, MaxEnt — 종 세기를 I∞(다층 표)로 나누고 계수 σ 에 상대 σ 바닥(기본 1 %, 제곱합)을 더해 넣는다: 계수 σ(≈0.1 %)만으로는 정확한 합성 자료도 층 나눔 때문에 χ² = N 에 닿지 못함(arxps_truth 의 무잡음 세기, J = 20: σ 0.3 % → α 하한 1e-4 에서 χ² 74 > N = 8, σ 1 % → α 2.4 에서 도달).
- **Imaging**: 창 · 넓이(+σ) · at %(같은 격자의 stack 여럿, `atomicMaps`) · 상태 지도(`subtractBackgroundStack` 뒤 무가중 또는 Poisson IRLS `nnlsMap`, σ 지도 = σ×√χ²ν(무가중)) · PCA 점수(창 원자료, MLPCA 부분 표본) → 재구성 stack. ROI = 지도 위 끌기 → 평균 ± σ, ‘Add to spectra’(c = cfac × 화소 수) · ‘Use as reference’(창에서 Shirley 를 빼고 넓이 1). PNG(화소 그대로 확대 + 색 막대) · CSV.
- **Process 탭**(Peak Fit 창의 자료): 고르지 않은 격자는 같은 평균 간격으로 자동 다시 표본화. 디컨볼루션(Worker; Jansson · RL · Landweber, PSF = Gauss · Voigt · 측정 스펙트럼 · Fermi edge/기준 봉우리에서 저장) — 결과는 그래프에 ± σ 띠와 가장자리 R 점 음영으로, 해상도 이득 · 잡음 이득 · 넓이 비 표시, **보기·시작값 전용**(정량은 원자료 피팅). Fermi edge(E_F ± σ, 분해능) · 기준 봉우리 Voigt(w_G) → ‘Use as PSF’. 다시 표본화(미리 보기 · 새 스펙트럼으로), SG/Gauss 평활 · SG 미분(척도 맞춰 겹침).

