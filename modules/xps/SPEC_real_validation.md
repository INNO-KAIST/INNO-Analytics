# 실측 대조 — Avantage · CasaXPS 출력 ↔ INNO (v1.8, 2026-10-28)

자료: 사용자 실측(비공개, `test/data/real/kara/` — 공개 백업 제외). 시험 `test/test_real_kara.js` (자료 없으면 건너뜀, 44 항목).
A = Avantage Peak Table + 영역별 'Backgnd.'(Smart) 열 4 파일(26 영역, Kratos 계 장비: Hybrid lens · Slot, VAMAS 에 투과 함수 열 있음).
B = CasaXPS 2.3.25 Integrated.vms(CASA region/comp/Calib 주석) + 'Fitting' 내보내기 6 개(성분 곡선·Background·Envelope).
C = CasaXPS Cookbook 2026 (Casa Software Ltd) 정의만 인용(쪽 표시). 두 프로그램의 코드는 쓰지 않았다. 원자료 수치는 이 문서에 싣지 않는다.

## A. Avantage

### A0 읽기 (고침)
| 항목 | 전 | 후 |
|---|---|---|
| Kratos→Avantage .xls (위쪽 'Axis · Kinetic energy · Elements=' 줄) | 읽기 실패(머리칸을 그 줄로 잡음) | 8 스펙트럼 읽음 (`avantageSheets`: 머리칸 아래 6 줄 안에 숫자가 있어야 함) |
| 시트 이름 'Sheet5…' | 그대로 | Peak Table 의 Start/End BE 가 배경 열 구간(±½ 채널)과 같은 줄의 Name (26/26) |
| 추가 열(Backgnd. 등) 순서 | BE 를 뒤집을 때 안 뒤집혀 채널이 어긋남 | 같이 뒤집음 |
.xls 세기 = .vms 세기/(dwell·scans) (6 자리 반올림 안, ≤ 2.1e-6) → Avantage 는 Kratos 투과 열을 쓰지 않았고 Peak Table TXFN = 1(전부).
**따라서 이 자료로는 Avantage 투과 함수 형태(SPEC §4 'Avantage 형태')를 정할 수 없다.**

### A1 Smart 배경 — 정한 것
- 배경 구간 = Start/End BE 에 **가장 가까운 채널** (26/26).
- 끝값 = 끝 채널 **±4 채널(9 점) 평균**, 창 밖 채널 포함 (26/26, 6 자리 인쇄 오차 안).
- 평활 9 점 이동평균 S 위의 Shirley: 배경 증분 회귀에서 9 개 같은 가중(±4) — 원자료가 아니라 평활 자료를 적분.
- 높은 BE 끝 ≤ 낮은 BE 끝인 구간은 직선(설명서와 같음; C 1s 두 영역이 정확히 직선).
- 배경이 평활 자료를 넘으면 그 최대 초과 채널을 새 끝점(값 = S)으로 나눠 다시 계산 → 조각 Shirley/직선. 첫 분할점 = 전체 Shirley 의 argmax(B − S) 가 분할이 있는 17 영역 중 14 에서 맞음.

### A2 Smart 대조 (`bg:{type:"smart", conv:"avantage"}`, 새 선택; 기본 Smart 는 그대로)
| | 배경 max\|ΔB\|/봉우리 높이 (중앙 · 최대) | 넓이(P) \|Δ\| (중앙 · 최대) |
|---|---|---|
| INNO 기본 Smart (n 5, 점별 자름) | 3.2e-2 · 0.25 | 4.33 % · 54.6 % (Cu 2p3) |
| **INNO Avantage 규약 Smart** | 3.9e-3 · 5.3e-2 | **0.41 % · 4.12 %** |
남은 차이: Fe 2p(창 23.8 eV) −1.9 … −4.1 %, C 1s +1.7 … +2.9 %(파일 1·2), 나머지 ≤ 1.5 %.
**미확인:** 분할이 없는 영역(O 1s ×2, Sr 3d)에서도 수렴 Shirley 와 0.5–5.6 CPS(계단의 1–3e-3) 다르다.
피적분을 S − λB − (1 − λ)y_lo 로 두면 인쇄 오차 수준(0.03–0.18 CPS)이 되지만 λ = 0.926–0.932(Fe 2p 0.956)로 상수가 아니어서 쓰지 않았다
(반복 횟수 제한·시작 곡선·κ 고정 변형도 시험 — 해당 없음). Avantage 의 Shirley 고정점 식이 다르다고만 말할 수 있다.

### A3 넓이·정량
| 항목 | 결과 |
|---|---|
| Area(P) = Avantage 배경 위 사다리꼴 적분(끝 = 가장 가까운 채널) | 중앙 5.4e-5 · 최대 7.6e-4 |
| ECF = Area(P)/(SF·TXFN·Area(N)) vs KE^n | n = 0.750, 잔차 6.0e-3 ≫ 반올림 6.4e-4 → 거듭제곱 아님 |
| ECF vs TPP-2M λ(KE = hν − Peak BE, 대전 보정 전) | ρ 6.22, Nv·ρ/M = 0.339 (Eg 0) 에서 2.1e-4 (반올림 안). Eg 와 U 는 서로 바뀜(물질 하나로 못 정함) |
| INNO `quantify` (Area(P), SF, TPP-2M, tf none) at% | 맞춘 물질 0.0067 · INNO 기본 물질(5, 24, 150, 3) 0.011 · KE^0.6 0.72 at% |
| Atomic % = Area(N)/Σ, PP At.% = PP Hgt(N)/Σ | ≤ 0.007 at% |
→ INNO 의 TPP-2M ECF 는 Avantage 'Area (N) TPP-2M' 을 재현한다(고칠 것 없음). Avantage at% 와 맞추려면 정량 탭 IMFP 를 TPP-2M 으로(기본 KE^0.6 이면 최대 0.7 at% 차).

### A4 정의 미확인 (INNO 에 대응 없음 — 참고만)
- PP Hgt(N) = PP Height/(SF·ECF·k): 한 줄 영역 k = 5.00, 2p 영역 k = 2.22–3.34(원소마다 다름).
- PP Height: 5/26 만 영역 안 원자료 max − min 과 같음.
- Height · Peak BE · FWHM('Standard' 봉우리): 평활 최댓값과 Height 중앙 1.7 %(최대 11 %), Peak BE 중앙 0.09 eV(최대 0.49 eV).
  Sr 3d 쌍에서 Peak BE 가 두 봉우리 사이·FWHM 이 둘을 덮음 — 단일 Gaussian/GL 맞춤도 아님(시험).

## B. CasaXPS

### B0 VAMAS 주석 읽기 (새 기능: `XPS_IO.casaCommentParse` → `meta.casaModel`, `XPS_CORE.casaModelToINNO`)
- `Calib M = … A = … BE ADD` 여러 줄은 **합**이 이동량(D.I. O 1s 세 줄 합 0.5761 eV = 내보내기 축 차, 4/4).
- CASA region/comp 의 KE 는 **보정 눈금**: 측정 BE = hν − KE − shift (영역 끝 채널·성분 위치로 확인).
- region: `KEstart KEend RSF AvWidth StartOffset EndOffset c1–c4(Cross Section) (*tag*) mass …`; comp: `Area|MFWHM|Position v lo hi link k`,
  link = 성분 번호(−1 없음), k = 배율 (3/2 넓이 = 0.66666667 × 5/2, 폭 연결 k 1). 위치 연결은 자료에 없어(모두 −1) 표기 미확인 → 변환 안 함·알림.

### B1 Shirley 배경 (`conv:"casa"`)
- 끝점 = **영역 안 첫·끝 채널** 중심 평균, 주석 Av.Width **2 → ±1 채널**(3 점). (주석 칸과 화면 'Av Width' 의 관계는 미확인 — /2 또는 −1.)
- 누적 = 배타 직사각 합(SPEC §2b 와 같음). `casaModelToINNO(…, {step})` 가 창을 한 채널씩 넓혀 INNO 규약(창의 2 번째·끝에서 2 번째 채널)에 맞춤.
| 내보내기 | max\|ΔB\|/봉우리 높이 |
|---|---|
| Bare Ba 4d · Bare O 1s · D.I. Ba 4d · D.I. O 1s | 1.8e-4 · 2.6e-4 · 1.5e-4 · 1.4e-4 (4–9 CPS) |
| D.I. Ba 4d GL(30) #1 · LA(50) #1 | 1.5e-4 · 3.7e-5 |
남은 수 CPS 는 반복 수·시작 곡선·끝점 이동으로 설명되지 않음(시험).

### B2 선모양 (Casa 매개변수에서 INNO 곡선 vs 내보내기 성분 열)
| 모양 | 성분 | max\|Δ\|/최댓값 | 비고 |
|---|---|---|---|
| GL(30) | 14 | ≤ 3.8e-4 | CasaXPS fwhm = 실제 FWHM 으로 F 환산. fwhm 을 +1.4…2.2e-4 키우면 ≈ 2e-4 (SPEC_casa_v2 §2.2 의 −1.7e-4 와 같은 방향) |
| LA(50) | 2 | 1.4e-3 (n 700.5) · 1.3e-3 (n 701) | 최소 ≈ 1.1e-3 (n 705–710) — SPEC_casa_v2 §4 의 Lorentzian 수치 차이; 700.5/701 구별 안 됨 |
Area = ℝ 넓이(영역 적분: 내보내기 = INNO).

### B3 포락선 · 다시 맞춤 (Casa 규약, 무가중, 채널 = 영역 안)
| 내보내기 | RMS: Casa 저장 / INNO@Casa 값 / INNO 다시 맞춤 | 다시 맞춤 넓이 변화 |
|---|---|---|
| Bare Ba 4d | 3701.7 / 3715.5 / 3458.1 | −1.4 … +1.7 % |
| Bare O 1s | 2976.3 / 2979.2 / 2972.9 | −4.9 … +2.4 % |
| D.I. Ba 4d (GL(30)+LA(50)) | 6373.5 / 6375.8 / 3339.0 | ±50 % |
| D.I. O 1s | 5395.7 / 5395.9 / 5286.2 | −0.3 … +1.0 % |
| D.I. Ba 4d GL(30) #1 · LA(50) #1 | 4916.9 / – / 4256.5 · 5971.3 / – / 2850.1 | |
INNO 가 Casa 값에서 내는 RMS 는 Casa 포락선과 0.4 % 안(모양·배경·연결 재현). **저장된 CasaXPS 모형은 같은 목적함수의 최소점이 아니다**(6/6 에서 INNO 가 더 낮음;
Poisson 가중으로 봐도 같음) → 이 파일로는 맞춤 결과(넓이·폭·위치) 일치를 시험할 수 없다.

## C. CasaXPS Cookbook 2026 이 정한 것 (정의만)
| 쪽 | 내용 | INNO |
|---|---|---|
| 19 | LA(50) = LA(1,701) = LA(1,1,701) | `casaToComp` 는 1401 − 14.01m = 700.5 (Fairley SI) — 실측으로 구별 안 됨 |
| 65–66 | lg(x;α,β) = l(x)^α (x ≤ 0), l(x)^β (x > 0), LA = lg ⊗ g(w); LA(α,w) = LA(α,α,w); TLA·ST(μ,γ) | x 는 KE 방향 → α = 높은 BE 쪽 = SPEC §1.6 와 같음. TLA·ST 없음 |
| 66 | DS = P(x,α) ⊗ G(n); DS Shirley 형(Q, k) | DS ⊗ G 있음, Shirley 형 없음 |
| 67 | GL(x,m) 곱 = SPEC §1.3 GLP; SGL = **높이** 혼합 (m/100)·L + (1 − m/100)·G | SPEC_casa_v2 §2.1 판정과 같음 |
| 67 | GL(m,n)T(k) = (GL + (1 − GL)·E(1,k)) ⊗ G(n), E(j,k) = j·exp(kx), x < 0; PHI(j,k) = E(j, 13.8/k) | AGT(tu 'fwhm', TM 0) 와 같은 식(2.2e-16). ⊗G(n)·PHI 없음 |
| 68 | SB(0)F = ∫ₓ^∞ F; TP(α,β,n) = R ⊗ L(n) | 없음 |
| 69–70 | Shirley 반복식(연속형), Av Width = 끝 채널 평균, Start/End Offset = 끝값 % 이동 | 이산 합 규약·Av Width 칸 의미는 정하지 않음 (B1 실측) |
| 70 | U 4 Tougaard: U(x; B,C,D,T0) = Bx/((C − x²)² + Dx²), x > T0 | INNO u3 핵 + T0 와 같음; `casaBg` 가 'U 4 Tougaard' 도 읽게 함 |
**정하지 못한 것:** 비대칭 모양의 화면 'fwhm'(GL(30)T(2) F 1.2 → 1.18051)과 'asymmetry index'(0.2138 · 0.4151 · 0.1973 · 0.2275)는 Cookbook 에 정의가 없다 —
모양 식만으로는 이 숫자가 나오지 않는다(SPEC_casa_v2 §2.2·2.4 그대로 미확인).

## D. 더 받으면 정해지는 것
1. Avantage **TXFN ≠ 1** 인 Peak Table(Thermo 장비, 투과 보정 켬) — Avantage 투과 함수 형태(`tfValue model:"avantage"`)와 VGD 계수.
2. Avantage **Shirley**(Smart 아님) 배경 열 1–2 영역 — 고정점 식 차이(A2)와 끝점 9 점 평균이 Shirley 에도 같은지.
3. CasaXPS 에서 **수렴까지 맞춘** 뒤 내보낸 Fitting 파일(같은 영역) — 넓이·폭·위치 대조. 위치 연결(오프셋) 쓴 모형 하나와 Av Width 를 1·3 으로 바꾼 주석 한 줄.
