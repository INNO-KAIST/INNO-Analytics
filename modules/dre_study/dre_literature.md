# 란타나이드 UC/PA 나노입자용 DRE(differential rate equation) 시뮬레이션·피팅 플랫폼 기획을 위한 문헌 조사

작성일: 2026-09-29
범위: 공개 문헌만 사용. 미공개 조성·샘플 정보는 검색·서술에 포함하지 않음.

표기 규칙
- **[verified]**: 이번 조사에서 원문/출판사 페이지/공개 코드/공개 SI로 직접 확인한 값·서지.
- **[secondary]**: 원문이 아니라 리뷰·요약 페이지에서 확인한 값.
- **[unverified]**: 기억 또는 간접 정보에 기반, 원문 대조 전. 플랫폼 기본값으로 쓰기 전에 원문 확인 필요.
- **[derived]**: 문헌 값으로부터 본 보고서에서 직접 계산/유도한 것(문헌 주장 아님).

---

## 0. 요약 (Executive summary)

1. **PA DRE의 사실상 표준(reference implementation)은 Lee et al. (Nature 2021)의 Tm³⁺ 3-level mean-field 모델**이며, 저자 공개 Python 코드(GitHub `nawhgnahc/Photon_Avalanche_DRE_calculation`)에 모든 수치가 들어 있다. 이 코드의 값은 Szalkowski et al. (Opt. Mater. X 2021)에 표로 재현되어 있어 교차검증이 가능했다(σ_GSA = 6.0×10⁻²⁹ m², σ_ESA = 6.4×10⁻²⁵ m², s₃₁ = 10 240 s⁻¹ @ 8 % Tm 등). 단, (i) 코드 주석의 단위(cm²)와 실제 사용 단위(m²)가 다르고, (ii) Q₂₂/Q₂₃ 라벨이 Szalkowski 표와 코드 사이에서 뒤바뀐 것으로 보이며, (iii) 정상상태를 0.5 s 적분으로 근사하는데 임계 부근 rise time이 ~0.6 s까지 늘어나므로 **임계 근처 정상상태가 수렴되지 않을 수 있다**. 플랫폼은 이 세 가지를 명시적으로 처리해야 한다.
2. **정상상태 파워 곡선(S-curve)만으로는 구조적으로 식별 불가능한 방향이 최소 1개 존재**한다: 모든 rate 상수와 흡수 단면적을 같은 배수 c로 스케일하면 정상상태 population이 불변이다 [derived]. 여기에 채널별 검출 gain(미지 scale)과 방출 준위의 radiative rate가 곱으로 묶인다. 따라서 (a) 한 개 이상의 절대 rate(예: ³F₄ 또는 ³H₄ 수명, JO radiative rate)를 고정하거나 (b) time-resolved(rise/decay) 데이터를 동시 피팅해야 한다. 이는 systems biology의 relative-data 식별성 문제(Weber 2011, Loos 2018)와 동형이다.
3. **최적화 전략**: 대규모 kinetic 모델 벤치마크(Raue 2013; Villaverde 2019)는 "log-parameter + 해석적/감도 기반 gradient + multistart local(또는 scatter search 하이브리드)"가 순수 확률적 전역탐색(PSO 등)보다 우수함을 보였다. 여기에 **scale/offset 파라미터의 계층적(hierarchical) 해석적 소거**(Loos 2018; Schmiester 2020)를 결합하는 것이 채널별 gain 미지 문제에 가장 직접적으로 맞는다. 식별성 진단은 profile likelihood(Raue 2009)로.
4. **Bistability**: Skripka et al. (Nat. Photon. 2025)은 KPb₂Cl₅:Nd³⁺에서 **열 효과 없이 순수 광학 DRE만으로** hysteresis가 재현됨을 보였다(시뮬레이션에 열 항 미포함). 즉 PA DRE는 일반적으로 다중 정상해를 가질 수 있고, 정상상태 solver는 **up/down sweep + numerical continuation(pseudo-arclength)** 으로 해 가지(branch)를 추적해야 한다. Lee 2021은 Tm:NaYF₄에서 up/down scan에 hysteresis가 없음을 보고.
5. **Mean-field의 한계**: 평균 거리 기반 mean-field는 고농도·소형 입자·비균질 도핑에서 실패 사례가 보고됨(Villanueva-Delgado 2016: average RE model이 decay를 재현 못 함, 571-ion microscopic model 필요; Teitelboim 2019: 고농도에서 DRE 적분이 수렴하지 않아 KMC 사용). 10 nm β-NaYF₄ 코어의 RE 사이트는 ~7 200개, 8 % Tm이면 ~580 Tm 이온 수준이다 [derived].
6. **크기/표면 효과**는 PA 문헌에서 대부분 "³F₄(k₂)·³H₄(k₃) 추가 비방사 rate"라는 현상론적 항으로 들어간다(Dudek 2022; Kwock 2021). Kwock 2021은 W₂,NR = Γ₀·exp[−λ(d_shell)]·(R_core+d_shell)²/R_core² 형태(Γ₀ = 2 380 s⁻¹, λ = 0.9 nm⁻¹)를 제시하고 **임계값 ∝ W₂/σ_ESA** 관계를 사용했다.

---

## 1. 모델 계열(Model families) 정리

### 1.1 PA 3-level mean-field 모델 (Tm³⁺, 1064 nm) — Lee 2021 / Szalkowski 2021 계열

**준위**: n₁ = ³H₆ (ground), n₂ = ³F₄ (intermediate, reservoir), n₃ = ³H₄ (800 nm 방출). ³H₅, ³F₂,₃ 등은 빠른 MPR로 흡수·소거되는 중간 경유 준위로 처리.

**Szalkowski et al. 2021의 식 (Opt. Mater. X 12, 100102)** [verified, 원문 식 7–9]:

```
dn3/dt = (I/hν) σ_ESA n2 − (A3 + k3) n3 − s31 n1 n3 + Q23 n2²
dn2/dt = (I/hν)(σ_GSA n1 − σ_ESA n2) − (A2 + k2) n2 + (β32 A3 + k3) n3
         + 2 s31 n1 n3 − (Q22 + 2 Q23) n2²
dn1/dt = −dn2/dt − dn3/dt        (n1+n2+n3 = 1)
```

**Lee 2021 공개 코드의 rate tensor**(`Time-resolved photoluminescence calculation.py`) [verified, 코드 직접 확인]:

```
A_fl = [[0,  w2,  (1−b32) w3_r        ],
        [0, −w2,   b32 w3_r + w3_nr   ],
        [0,  0,   −w3_r − w3_nr       ]]
A_cr(n) = [[0, (Q22+Q23) n2,      −s3 n1],
           [0, −(Q22+2Q23) n2,   2 s3 n1],
           [0,  Q23 n2,           −s3 n1]]
pump = φ · [[−σG, 0, 0],[σG, −σE, 0],[0, σE, 0]],   φ = I/(hν)
dn/dt = (A_fl + A_cr(n) + pump) · n
```

해석: ³H₄의 비방사 감쇠(w3_nr)는 전부 ³F₄로 가고(³H₄→³H₅→³F₄), 방사 감쇠 중 b32 = 0.144 만 ³F₄로, 나머지는 ground로. CR: (³H₄, ³H₆)→(³F₄, ³F₄). Q23: (³F₄, ³F₄)→(³H₄, ³H₆) ETU. Q22: (³F₄, ³F₄)→(³H₆, ³H₅→³F₄), 즉 순수 손실형 ETU.

**Gaussian beam 평균** [verified, 코드]: 앙상블 신호 = ∫₀^{I₀} n₃(I)/I dI (`cumtrapz(n3/P, P)`). 이는 2D Gaussian 강도분포 I(r) = I₀exp(−2r²/w²)에서 dA ∝ dI/I 이므로 정확한 공간 적분이다 [derived]. 단일입자 측정(Kwock 2021)에서는 peak intensity를 그대로 사용 — **측정 형태(단일입자 vs 필름/앙상블)에 따라 forward model이 달라져야 함**.

**수치 적분**: 고전 RK4 + 적응 step dt = 0.4/max|A_total|, 정상상태는 t = 0.5 s 까지 적분한 값으로 취함 [verified, 코드]. Lee 2021 본문은 임계 부근 rise time이 최대 **~608 ms**(³F₄ 수명의 ~400배)라고 보고 [verified, eScholarship 원고]. → **0.5 s 적분은 임계 근방에서 정상상태에 도달하지 못할 수 있음** [derived]. 코드에는 `time_target_2`가 함수 내부에서 전역변수로 참조되는 등 재사용 시 주의할 점도 있음.

**정의** (Lee 2021) [verified]: slope s는 I_E ∝ I_P^s; rise time은 점근값의 95 % 도달 시간; threshold 이상에서 2배 파워 증가 시 ~500–10 000배 증가; 8 % Tm ANP 피팅에서 ESA/GSA(R₂/R₁) 비 ≈ 13 000 (SI Table S4); 나노입자는 Tm³⁺ 도핑 코어/비활성 셸, 총 직경 16–33 nm, Tm 1–100 %. "up/down scan에서 photobleaching·hysteresis 없음 → 열적 avalanche 기여 없음".

### 1.2 PA의 고전 이론 (bulk)

- **Chivian, Case, Eden (APL 1979)**: Pr³⁺:LaCl₃/LaBr₃ IR quantum counter에서 PA 최초 보고 [서지 verified(검색), DOI 10.1063/1.91044 unverified].
- **Goldner & Pellé (Opt. Mater. 5, 239, 1996)** [verified]: 3-level 정상상태 모델. "pump photon과 전자 전이 간 큰 mismatch가 필요", avalanche threshold 조건은 화합물에 거의 무관하며 최소 GSA에 해당, GSA/ESA 비는 multiphonon absorption 계산으로 예측 가능. PA 펌핑은 공진 펌핑보다 UC 세기가 낮다는 결론.
- **Joubert (Opt. Mater. 11, 181, 1999)** [verified]: ESA/ETU/PA 리뷰, PA를 임계 이상에서만 일어나는 bifurcation 현상으로 기술, rate-equation 모델링으로 bulk 결과 재현.
- **Guy, Joubert, Jacquier (PRB 55, 8240, 1997)** "Photon avalanche and the mean-field approximation" [서지 verified(APS 페이지 제목/권호), 내용 unverified]: PA에서 mean-field 가정의 타당성을 다룬 논문. 원문 접근 실패(403). 플랫폼에서 "ion-pair/비균질 CR이 mean-field와 어떻게 다른가"의 1차 문헌으로 확인 필요.
- **Pollnau et al. (PRB 61, 3337, 2000)** [verified]: n-photon 순차 UC의 파워 의존성은 UC rate와 선형 감쇠의 경쟁에 따라 Pⁿ ~ P¹ 사이 어디든 가능. ETU vs ESA를 slope로 구분하는 기준 제공. → **S-curve의 저파워 slope도 파라미터 조합에 의존**하므로 slope 자체를 "photon 수"로 해석하면 안 됨.

### 1.3 PA 최신 모델 변형

| 논문 | 시스템 | 모델 특징 | 비고 |
|---|---|---|---|
| Bednarkiewicz et al., Nanoscale Horiz. 4, 881 (2019) | Nd³⁺, 1064 nm | ⁴I₉/₂, ⁴I₁₁/₂, ⁴F₃/₂ 중심 DRE(SI S1–S5), W_CR 농도 의존, 분기비는 JO | GSA/ESA 비 R = 0.98×10⁻⁴ (Nd:NaYF₄ @1064 nm), time-gating 3–6 ms, N=80이면 ~41 nm 분해능 예측 [verified] |
| Szalkowski et al., Opt. Mater. X 12, 100102 (2021) | Tm:NaYF₄ | Lee 모델 + 온도 의존 MPR/PAET (Miyakawa–Dexter) | 300 K 상대 온도감도 최대 ~15 %/K, >400 K 에서 >40 %/K; 지배 인자는 ³F₄→³H₆ MPR 온도의존 [verified] |
| Kwock et al., JPCC 125 (2021) | Tm:NaYF₄ 단일입자 | Lee DRE + 셸두께 의존 W₂,NR | Γ₀ = 2 380 s⁻¹, λ = 0.9 nm⁻¹, W₂,R = 83.3 s⁻¹ (고정) [verified] |
| Dudek et al., Adv. Opt. Mater. 10, 2201052 (2022) | Tm:LiYF₄ nano/micro/bulk | 3-level DRE, 크기 효과를 k₂, k₃ 증가와 f_CR 감소로 표현 | threshold ~200 (bulk) → 300–400 (micro) → 700–800 kW/cm² (nano); S ≈ 6 → 8 → 10–12 (3 % Tm); t₅₀ 100 ms → 9 ms [verified] |
| Liang et al., Nat. Nanotechnol. 17, 524 (2022) | Yb/Pr 코어 + Ho/Tm 셸, 852 nm | migrating PA; Extended Data Fig. 1에 population 시뮬레이션 | 코어 26차, 최대 46차, threshold 60 kW/cm² [verified] |
| Dudek et al., Nanoscale 15, 18613 (2023) | Yb³⁺-sensitized Pr³⁺, 852 nm | 단순화 준위(Pr: ³H₄,³H₅,³H₆,¹G₄,³P₀; Yb 2준위), 파라미터 대부분 LiYF₄ 문헌값, MPR은 Kaminskii chart에서 | one-at-a-time 감도분석: σ_ESA, w65(ETU), τ_Yb, w64(back-transfer), τ₄가 지배적; NR₅₄, NR₄₃은 ±3 decade 변화에도 영향 미미; W_CR1은 w65에 가려 "PA 특성에 영향 없음" → **degeneracy의 직접 사례** [verified] |
| Skripka et al., Nano Lett. 23, 7100 (2023) | Tm 코어→Gd 버퍼→활성이온 | 메커니즘은 SI S9 시뮬레이션 | s = 10–17, Er 6 kW/cm², Ho 8 kW/cm² [verified 요약] |
| Skripka et al., Nat. Photon. 19, 212 (2025) | KPb₂Cl₅:Nd³⁺ | ~40 준위 전체 모델 + 7 준위 단순화; 열 항 없음 | IOB(intrinsic optical bistability) 비열적 기원, >200차 비선형, 77 K; KPb₂Cl₅ MPR W₀ = 0.5×10⁷ s⁻¹, ħω = 250 cm⁻¹; CR1 microparameter 363×10⁻⁴⁰ cm⁶/s; 77 K ⁴I₁₁/₂→⁴I₉/₂ MPR 185.85 s⁻¹ (NaYF₄에서는 654 905 s⁻¹로 비교) [verified, arXiv 2403.04098; MPR 식의 정확한 함수형은 추출 과정에서 왜곡됐을 가능성 있어 원문 확인 필요] |

### 1.4 Yb/Er, Yb/Tm UCNP의 표준 DRE 계열

- **Anderson, Smith, May, Berry, JPCL 5, 36 (2014)**; **Hossan et al., JPCC 121, 16592 (2017)** (DOI 10.1021/acs.jpcc.7b04567, verified): β-NaYF₄:Yb,Er 코어/코어-셸의 나노스케일 효과를 rate-equation으로 설명. 파라미터 표는 이번 조사에서 원문 접근 실패 → [unverified].
- **Würth et al., Nanoscale 9, 4283 (2017)**; **Kaiser et al., Nanoscale 9, 10051 (2017)**; **Kaiser et al., Nano Res. 12, 1871 (2019)** (verified 서지): 절대 양자수율(파워 의존)과 decay kinetics를 결합한 coupled rate equation 피팅. 파라미터 표는 SI(MOESM)에 있음 — 이번 조사에서 수치 추출 실패.
- **Villanueva-Delgado, Krämer, Valiente, JPCC 119, 23648 (2015)** [서지 verified]: β-NaYF₄:Yb,Tm ET/UC 시뮬레이션(microscopic 접근). 본문 수치 미확인.
- **Chan et al., JPCB (2012), DOI 10.1021/jp302401j** 및 **Chan, Chem. Soc. Rev. (2015) combinatorial review, DOI 10.1039/C4CS00205A**: JO 기반 radiative rate + Kushida 식 ET microparameter + Miyakawa–Dexter PAET로 수천 개 rate를 자동 생성하는 "Chan-type" DRE. SI에 식(S1–S10) [verified]:
  - A_ED = 64π⁴e²ν̃³ / [3h(2J+1)] · [n(n²+2)/9]² · S_ED, S_ED = Σ_λ Ω_λ |⟨U^(λ)⟩|²
  - W_ET = C_DA / R⁶,  C_DA ∝ s₀ · [(n²+2)/(3n)]⁴ · S_ij · S_kl /(g_i g_k)
  - C_DA^PAET = C_DA⁰ · exp(−β_MPR ΔE)
  - 이 계열이 Lee 2021 모델의 원형이며, Teitelboim 2019 KMC와 같은 rate 상수를 공유.

### 1.5 크기/표면 소광 모델

- **Johnson et al., JACS 139, 3275 (2017)** (주의: 요청서의 "ACS Nano"가 아니라 **JACS**; DOI 10.1021/jacs.7b00223, verified): 표면 소광과 농도 소광이 energy migration을 통해 결합되어 있다는 직접 증거.
- **Fischer, Bronstein, Swabeck, Chan, Alivisatos, Nano Lett. 16, 7241 (2016)** [verified]: 비활성 β-NaLuF₄ 셸 두께를 체계적으로 바꿔 준위별 표면 소광 rate를 측정; 셸 두께에 대한 효율 최적점.
- **Rabouw et al., ACS Nano 12, 4812 (2018)** (DOI 10.1021/acsnano.8b01545) [서지 verified, 권·쪽 unverified]: NaYF₄:Er,Yb 소광 경로; 리간드/용매 진동으로의 에너지 전달.
- **Wang & Meijerink, JPCC 122, 26298 (2018)** [verified]: 농도 소광 — Er ⁴S₃/₂, ⁴I₁₁/₂에서 강함, Yb는 코어-셸에서 60 %까지 효율 유지, 수계에서 O–H 진동에 의한 강한 소광.
- **Kwock 2021 셸 모델** [verified]:
  W₂,NR(d_shell) = Γ₀ · exp(−λ d_shell) · (R_core + d_shell)² / R_core²,  Γ₀ = 2 380 s⁻¹, λ = 0.9 nm⁻¹ (κ⁻¹ ≈ 1.1 nm).
  CR 지배 영역(Tm ≥ 8 %)에서 R₂,th ≈ σ_ESA I_th = W₂ → **I_th ∝ W₂**. 단일입자 threshold 분포는 높은 쪽으로 skew; 평균 셸 두께가 아니라 **가장 얇은 부분**이 I_th를 결정한다는 해석.
- **Dudek 2022**: 크기가 작아질수록 k₂, k₃ 증가, f_CR 감소로 표현. 나노에서 threshold↑ 이면서 S↑ — 표면 손실이 오히려 pre-avalanche를 억눌러 비선형성을 키우는 경향.
- **Cation intermixing**: 코어/셸 계면에서 양이온 상호확산(Hudry et al., Chem. Mater. 29, 9238 (2017)) [unverified DOI] — "셸 두께 = 거리" 가정의 오차원.

크기가 rate equation에 들어가는 방식(문헌에서 쓰인 세 가지):
1. **현상론적 추가 비방사 rate** (k_surf,i): Dudek 2022, Kwock 2021. 가장 단순, 파라미터 1–2개.
2. **Core–shell 거리 의존 소광 + 표면 부피비**: Kwock 식, Fischer 2016.
3. **공간 분해(shell-resolved 또는 per-ion) microscopic 모델**: Villanueva-Delgado 2016, Teitelboim 2019 KMC. 에너지 migration을 명시적으로 표현.

### 1.6 열 효과와 bistability

- **열 없이도 bistability 가능**: Skripka 2025 — MPR이 dark state를 안정화(초기 소수 여기 이온을 조기 소광), ESA+CR 양의 피드백이 bright state를 안정화; switch-on에서 ⁴I₉/₂와 ⁴I₁₁/₂ population 역전. 77–150 K에서 IOB, 160 K에서 소멸(s = 70은 유지). dwell time 0.1–10 s 변화에도 hysteresis 폭 불변 → 동역학적 지연(rise time)이 아닌 진짜 다중 정상해 [verified].
- **열적 bistability**: Yb³⁺ bromide에서 intrinsic optical bistability가 처음 협동 효과로 보고(Hehlen et al., PRL 73, 1103 (1994)) [서지 verified]; 이후 Gamelin, Lüthi, Güdel (JPCB 104, 11045 (2000))이 **레이저 가열의 역할**을 분석 [서지 verified, 결론 세부 unverified]. → 플랫폼은 "광학 bistability"와 "열 bistability"를 구분하는 옵션(온도 ODE 결합)을 가져야 함.
- **온도 의존 PA**: Szalkowski 2021(예측), Korczak et al. Low Temp. Phys. 49, 322 (2023)(Yb/Pr 실험: −175~175 °C에서 S 4.5–9, threshold 100–700 kW/cm², 0 °C에서 7.5 %/°C) [secondary, Nanoscale 2025 리뷰 경유]. NdAl₃(BO₃)₄에서 온도가 PA-like 메커니즘을 촉발하는 보고(J. Lumin., 2021, S0022231321007614) [서지 unverified].
- 함의: RT에서 PA 곡선의 상대 온도감도가 ~10 %/K 수준이면, 수 K의 레이저 가열만으로도 threshold가 눈에 띄게 이동한다 [derived]. 피팅 전 up/down sweep hysteresis 체크가 필수(Lee 2021의 절차).

---

## 2. 파라미터 표 (출처 포함)

### 2.1 Lee et al. 2021 공개 코드의 Tm³⁺:NaYF₄ PA 파라미터 (1064 nm = 9398 cm⁻¹)

코드 URL: https://github.com/nawhgnahc/Photon_Avalanche_DRE_calculation (파일 `Time-resolved photoluminescence calculation.py`, 2020-06-17 작성). 전부 [verified, 코드 직접 확인]; 파생값은 [derived].

| 기호 | 코드 식 | 값 | 단위 | 의미 |
|---|---|---|---|---|
| w2_r | 상수 | 162.60 | s⁻¹ | ³F₄ radiative |
| w2_nr | 162.60·f_w2_nr + 0.00239 | 샘플별(아래) | s⁻¹ | ³F₄ nonradiative |
| w3_r | 상수 | 636.01 | s⁻¹ | ³H₄ radiative |
| w3_nr | 6.02 + 162.60·f_w3nr | 샘플별 | s⁻¹ | ³H₄ nonradiative (→³F₄) |
| b32 | 상수 | 0.144 | – | ³H₄→³F₄ 방사 분기비 |
| σ_GSA | f_ab_1·0.12×10⁻²⁵, f_ab_1 = 0.005 | 6.0×10⁻²⁹ m² (= 6.0×10⁻²⁵ cm²) | m² | 비공진 GSA (코드 주석은 cm²라 적었으나 photon flux가 m⁻²s⁻¹ 이므로 실제 단위는 m²) |
| σ_ESA | f_ab_2·3.2×10⁻²⁵, f_ab_2 = 2.0 | 6.4×10⁻²⁵ m² (= 6.4×10⁻²¹ cm²) | m² | ³F₄→³F₂,₃ ESA |
| σ_ESA/σ_GSA | | 10 667 | – | [derived]; 논문 본문 "≈13 000"(8 % 피팅, SI Table S4)과 다름 — 버전 차이 가능 |
| s31 | f_s3·1.6·c²·1000, f_s3 = 0.1 | 10 240 @ c = 8 | s⁻¹ | CR (³H₄,³H₆)→(³F₄,³F₄); **c(%)의 2승 스케일** |
| Q22 | 0.08·320·c³/(c² + 4.3²) | 158.9 @ c = 8 | s⁻¹ | ETU (³F₄,³F₄)→(³H₆,³H₅→³F₄) |
| Q23 | 0.1·90·c³/(c² + 4.3²) | 55.9 @ c = 8 | s⁻¹ | ETU (³F₄,³F₄)→(³H₆,³H₄) |
| 적분시간 | time_target | 0.5 | s | 정상상태 근사 |
| 파워 범위 | 10^7 – 10^9.1 | W/m² | = 1–126 kW/cm² |

샘플별 계수 (코드의 `c_Tm_array`, `f_w2_nr_array`, `f_w3_nr_array`) [verified]; w_nr 값은 [derived]:

| 라벨(코드) | c_Tm (코드 입력값) | f_w2_nr | f_w3nr | w2_nr (s⁻¹) | w3_nr (s⁻¹) |
|---|---|---|---|---|---|
| "4 %, Sample No. 2" | 2.25 | 0.35 | 0.6 | 56.9 | 103.6 |
| "8 %, Sample No. 3" | 8.0 | 3.25 | 6.3 | 528.5 | 1030.4 |
| "8 %, Sample No. 4" | 8.0 | 0.25 | 0.5 | 40.65 | 87.3 |
| "8 %, Sample No. 5" | 8.0 | 0.0 | 0.0 | 0.0024 | 6.02 |
| "20 %, Sample No. 6" | 20.0 | 6.0 | 14.0 | 975.6 | 2282.4 |
| "20 %, Sample No. 7" | 20.0 | 3.6 | 7.2 | 585.4 | 1176.7 |

관찰:
- "4 %" 샘플에 c = 2.25를 입력 — 명목 농도와 모델 유효 농도가 다름(실측 조성 보정일 가능성; 코드에 설명 없음).
- 같은 명목 농도(8 %)에서 w_nr이 샘플 간 최대 ~10²배 차이 → **샘플별로 자유롭게 두는 파라미터는 사실상 비방사 손실(셸 품질)뿐**이고 σ, s31, Q는 공유(global) 파라미터로 취급됨. 플랫폼의 multi-sample global fit 설계의 선례.
- Szalkowski 2021 Table 1의 300 K 값 k₂ = 40.7 s⁻¹, k₃ = 87.3 s⁻¹, s₃₁ = 10 240 s⁻¹는 "Sample No. 4"와 정확히 일치 [derived 교차검증].
- **Q22/Q23 라벨 불일치**: 코드에서는 Q22 = 158.9 (손실형), Q23 = 55.9 (→³H₄). Szalkowski 표(추출본)는 Q₂₂ = 55.8, Q₂₃ = 158.9로 표기 — 라벨이 뒤바뀌었거나 추출 오류 가능. 원문 PDF로 최종 확인 필요.
- 파생 스케일 [derived]: Sample 4의 W₂ = 162.6 + 40.65 = 203 s⁻¹ (τ(³F₄) ≈ 4.9 ms). 특성 강도 I* = W₂hν/σ_ESA ≈ 5.9 kW/cm². Kwock의 R₂,th ≈ W₂ 관계와 차원적으로 일관.
- Kwock 2021은 W₂,R = 83.3 s⁻¹ (τ_R = 12 ms)를 사용 — Lee 코드(162.6 s⁻¹)와 2배 차이. **³F₄ radiative rate조차 논문마다 2배 차이**가 있으므로 플랫폼 기본값은 반드시 출처와 함께 저장해야 함.

### 2.2 Szalkowski et al. 2021 (Opt. Mater. X 12, 100102) — 온도 의존 파라미터 [verified]

| 과정 | 300 K rate | 전이 | ΔE (cm⁻¹) | 선택 ħω (cm⁻¹) | P (phonon 수) |
|---|---|---|---|---|---|
| k₂ (MPR) | 40.7 s⁻¹ | ³F₄→³H₆ | 5180 | 305 | 17 |
| k₃ (MPR) | 87.3 s⁻¹ | ³H₄→³H₅(→³F₄) | 4064 | 365 | 11 |
| s₃₁ (CR, PAET) | 10 240 s⁻¹ | (³H₄;³H₆)→(³F₄;³F₄) | 1401 | – | 4 |
| Q₂₂ | 55.8 s⁻¹ (라벨 주의) | (³F₄;³F₄)→³H₅→(³F₄;³H₆) | 741 | – | 2 |
| Q₂₃ | 158.9 s⁻¹ (라벨 주의) | (³F₄;³F₄)→(³H₄;³H₆) | 1028 | – | 3 |
| σ_ESA | 6.41×10⁻²⁵ m² | ³F₄→³F₂,₃ | | | |
| σ_GSA | 6.00×10⁻²⁹ m² | ³H₆→³H₅(→³F₄) | | | |
| β_PA | 10 666 | σ_ESA/σ_GSA | | | |

NaYF₄ 활성 phonon: 250, 305, 365 cm⁻¹. 온도 의존식 (원문 식 6):
W(ΔE, T) = W(ΔE = 0, T = 0) · exp(−α ΔE) · [±1 ∓ exp(∓ħω/k_BT)]^(−P)
(위 부호: phonon 방출, 아래 부호: phonon 흡수), n(T) = [exp(ħω/k_BT) − 1]⁻¹. α(MPR)와 β(PAET)는 서로 관련된다고만 서술; NaYF₄용 α, β, W(0) 수치는 추출본에 없음.

### 2.3 Tm³⁺:β-NaYF₄ 분광 파라미터 — Villanueva-Delgado et al., PCCP 18, 27396 (2016) [verified]

β-NaYF₄:0.01 % Tm (저농도, ET 무시 가능) 실측 수명·분기비(Table 1)와 ET microparameter:

| 준위 | 수명 | 주요 분기비 β_i,f |
|---|---|---|
| ¹D₂ | 67.5 μs | 0.40 (→³F₄), 0.09 (→³F₃) |
| ¹G₄ | 758 ± 1 μs | 0.07 (→³F₄), 0.22 (→³H₅), 0.07 (→³H₄) |
| ³F₃ | 20 μs | 0.999 (→³H₄) |
| ³H₄ | 1.94 ± 0.01 ms | 0.3 (→³F₄), 0.06 (→³H₅) |
| ³H₅ | 29 μs | 0.8 (→³F₄) |
| ³F₄ | 12 ms | – |

(주의: 위 분기비는 추출된 일부 값이며 나머지는 ground로 가는 것으로 해석해야 함. 원문 Table 1 전체 확인 권장.)

| 파라미터 | 값 | 비고 |
|---|---|---|
| CR microparameter W₆,₁ (d–d) | 4.205×10⁻³⁹ cm⁶ s⁻¹, R_c,CR = 12.1 Å | Scheme I (채택) |
| ETU microparameter W₄,₆ (d–d) | 1.258×10⁻³⁸ cm⁶ s⁻¹, R_c,ETU = 16.9 Å | Scheme I |
| (기각) ETU W₆,₆ | 1.60×10⁻³⁵ cm⁶ s⁻¹, R_c = 47.9 Å | 비물리적 → 기각 |
| Inokuti–Hirayama (n = 6) | C_SA = 3.2×10⁹ Å⁶ s⁻¹, R_c = 11.80 Å | n = 8: 5.5×10¹¹ Å⁸s⁻¹, 11.97 Å; n = 10: 2.0×10¹⁴ Å¹⁰s⁻¹, 13.10 Å |
| σ(³H₆→¹G₄) @ 21 140 cm⁻¹ | 7.77×10⁻²² cm² | 편광 평균 |
| 격자상수 β-NaYF₄ | a = 5.9734 Å, c = 3.5296 Å | |
| MPR ³F₃→³H₄ (Scheme II 피팅값) | 3.26×10⁸ s⁻¹ | 에너지갭 법칙 추정치(ΔE = 1800 cm⁻¹에서 4.61×10⁴ s⁻¹, Fischer et al. 파라미터 사용)보다 수 자리 큼 → Scheme II 기각 근거 |

모델 비교 결론: 571 이온 × 7 준위 = 3 997 변수의 microscopic RE가 decay를 재현(normalized RMSE 0.49), **average(mean-field) RE는 decay 곡선을 재현하지 못함**; energy migration ≈ 0. → 저농도·이 준위들에서는 **거리 분포가 핵심**이며 평균거리 mean-field는 비지수 감쇠를 설명 못 함.

**교차 비교 (³F₄ radiative/total rate)** [derived]: Villanueva(0.01 %) τ = 12 ms ↔ Kwock W₂,R = 83.3 s⁻¹ (12 ms, 일치) ↔ Lee 코드 w2_r = 162.6 s⁻¹ (6.15 ms). ³H₄: Villanueva 1.94 ms ↔ Lee w3_r = 636 s⁻¹ (1.57 ms).

### 2.4 KPb₂Cl₅:Nd³⁺ (Skripka 2025, arXiv 2403.04098) [verified, 일부 함수형 확인 필요]

| 항목 | 값 | 단위 |
|---|---|---|
| ħω | 250 | cm⁻¹ |
| MPR prefactor W₀ | 0.5×10⁷ | s⁻¹ |
| CR1 microparameter | 363×10⁻⁴⁰ | cm⁶ s⁻¹ |
| ⁴I₁₁/₂→⁴I₉/₂ MPR @77 K | 185.85 | s⁻¹ |
| 동일 전이 NaYF₄ 비교치 | 654 905 | s⁻¹ |
| 평균 switch-on power | 5.2 | kW/cm² |
| hysteresis 폭 | 1 (CW), 4.5–24 (pulsed) | kW/cm² |
| 시뮬레이션 rise time | 305 μs – 20 ms (저온), RT 50–70 ms | |

### 2.5 기타 PA 수치

| 항목 | 값 | 출처 |
|---|---|---|
| Nd:NaYF₄ @1064 nm GSA/ESA 비 | 0.98×10⁻⁴ | Bednarkiewicz 2019 [verified] |
| PA 기준 σ_ESA/σ_GSA | ≥10⁴ (단, cross-section 비인지 rate 비인지 문헌이 모호) | Skripka & Chan, Mater. Horiz. 2025 [verified] |
| Tm:LiYF₄ threshold | ~200 (bulk), 300–400 (μ), 700–800 (nano) kW/cm² | Dudek 2022 [verified] |
| Tm:NaYF₄ 단일입자 평균 I_th | 23 (20 % Tm, 2.6 nm 셸), 14.8 (8 %, 5.6 nm), 6 (8 %, 8.5 nm) kW/cm² | Kwock 2021 [verified] |
| 단일입자 예 | W₂ = 383.2 ± 3.7 s⁻¹ → I_th = 12.3 ± 0.1 kW/cm² | Kwock 2021 [verified] |
| Tm:NaYF₄ ANP 최대 rise time | ~608 ms (95 % 기준) | Lee 2021 [verified] |
| Dudek 2022 파라미터 스캔 격자 | A₂ ∈ {250…20}, A₃ ∈ {1000…25}, k₂ ∈ {1000…25}, k₃ ∈ {2000…25} s⁻¹; f_ESA ∈ {5…0.01}, f_CR ∈ {0.5…0.01} | [verified] — 파라미터 prior 범위 참고용 |

### 2.6 에너지갭 법칙 / host 파라미터

| 호스트 | 파라미터 | 출처/상태 |
|---|---|---|
| NaYF₄ | 활성 phonon 250/305/365 cm⁻¹ | Szalkowski 2021 [verified] |
| NaYF₄ | "fluoride 300–500 cm⁻¹" | Aggarwal et al. Nanoscale 2025 리뷰 [secondary] |
| KPb₂Br₅ | ~128 cm⁻¹ | 同 리뷰 [secondary] |
| KPb₂Cl₅ | W₀ = 0.5×10⁷ s⁻¹, ħω = 250 cm⁻¹ | Skripka 2025 [verified] |
| LaF₃ | B ≈ 6.6×10⁸ s⁻¹, α ≈ 5.6×10⁻³ cm, ħω_max ≈ 350 cm⁻¹ | Riseberg & Moos 계열로 널리 인용 **[unverified — 기억 기반]** |
| LiYF₄ | B ≈ 3.5×10⁷ s⁻¹, α ≈ 3.8×10⁻³ cm, ħω ≈ 400 cm⁻¹ | **[unverified — 기억 기반]** |
| LiYF₄ | MPR은 Kaminskii의 ΔE–rate chart에서 읽음 | Dudek 2023 Nanoscale [verified 방법론] |

**식의 형태**:
- 기본 에너지갭 법칙: W_MPR = B exp(−α ΔE) · [n(T)+1]^p, p = ΔE/ħω (Riseberg & Moos, Phys. Rev. 174, 429 (1968) [서지 verified]; Lehigh 강의노트는 "C, α, p는 host 의존 경험 상수이며 RE 이온에는 둔감"이라고 정리 [verified]).
- van Dijk & Schuurmans (J. Chem. Phys. 78, 5317 (1983)) [서지 verified via Crossref]: "modified exponential energy gap law". 흔히 W = β_el exp[−α(ΔE − 2ħω_max)] 형태로 인용됨 **[함수형 unverified — 원문 접근 실패]**.
- Miyakawa & Dexter (PRB 1, 2961 (1970)) [서지 verified]: 다중 phonon 이완과 phonon-assisted ET의 에너지 불일치 의존성을 같은 틀에서 유도; W_PAET(ΔE) = W_PAET(0) exp(−β ΔE)이고 β는 α와 관련(흔히 β = α − γ, γ = ln(1 + g_A/g_S 류)/ħω 형태로 인용) **[정확한 관계식 unverified]**.

플랫폼 제안 [derived]: MPR/PAET rate를 자유 파라미터로 두지 말고, **host당 (B, α, ħω) 3개 + 준위별 ΔE**로 파라미터화하면 수십 개의 MPR rate를 3개로 줄일 수 있다(Chan-type 모델이 실제로 이렇게 함). 온도 의존도 자동으로 따라온다.

### 2.7 농도 스케일링

- **Mean-field 기본**: 도너 1개당 ET rate는 Σ_j C_DA/R_j⁶ over acceptor sites → 희박 극한에서 acceptor 농도에 **선형**. 정규화 population(n_i = 준위 i의 분율)으로 쓴 이차항 s·n₁n₃에서 s가 농도에 선형 [derived, 표준 결과].
- **Lee 코드**: s31 ∝ c² (경험적), Q ∝ c³/(c² + 4.3²) (저농도 c³, 고농도 c¹로 포화) [verified 코드; 물리적 근거는 코드에 명시 없음]. Dudek 2022는 CR을 f_CR 스케일 인자로 조절.
- **Inokuti–Hirayama (J. Chem. Phys. 43, 1978 (1965))** [서지 verified]: 도너 감쇠 I(t) = exp[−t/τ₀ − Γ(1 − 3/s)(c_A/c₀)(t/τ₀)^{3/s}] (s = 6, 8, 10), c₀ = 3/(4πR₀³) — 표준식 [식 형태는 교과서 수준 지식]. Villanueva 2016이 이 모델로 R_c 추출(위 표).
- **Kushida/Dexter microparameter**: C_DA ∝ s₀ · S_ij · S_kl (Chan SI 식 S9) [verified].
- 농도 의존을 "계수 × 함수형"으로 고정하고 여러 농도 샘플을 동시에 피팅하는 것이 식별성을 크게 높임(Lee 코드의 방식) [derived].

### 2.8 Detailed balance / McCumber

McCumber (Phys. Rev. 136, A954 (1964)) [서지 verified]: σ_em(ν) = σ_abs(ν) · (Z_l/Z_u) · exp[(E_zl − hν)/k_BT]. PA 모델에서 ESA cross-section을 emission 스펙트럼으로부터 추정하거나, 역방향 ET(back-transfer)를 detailed balance로 묶는 데 쓸 수 있다. 여기서 ET의 역방향 rate W_back = W_fwd · (g_i g_k / g_j g_l) · exp(−ΔE/k_BT) 형태로 묶으면 자유 파라미터가 줄어든다 [derived, 표준 열역학 관계].

---

## 3. 식별성(identifiability)과 전역 최적화

### 3.1 PA S-curve에서 구조적 비식별 방향 [derived]

위 3-level 모델의 정상상태 조건 0 = [A_fl(θ) + A_cr(n;θ) + φ P(σ)] n 은 **모든 rate(A, k, s, Q)와 σ를 같은 인자 c로 곱해도 불변**이다(우변이 전체 c배가 될 뿐). 따라서:

1. 정상상태 데이터만으로는 **절대 시간 스케일이 식별 불가**. 식별 가능한 것은 무차원 조합뿐:
   - σ_ESA/σ_GSA
   - W₃/W₂, s₃₁/W₂, Q₂₂/W₂, Q₂₃/W₂, β₃₂
   - 강도 스케일 I* = W₂ hν/σ_ESA (x축이 절대 강도로 보정되어 있을 때만)
2. 측정 신호 = η_ch · A₃,rad · n₃ (η_ch: 채널 검출 효율×gain) → η_ch와 A₃,rad는 곱으로만 등장.
3. 따라서 **최소 한 개의 절대 rate(예: 저농도 수명, JO A_rad)를 고정**하거나 **시간분해 데이터(rise/decay)를 동시 피팅**해야 절대값이 나온다.
4. x축(여기 강도) 교정 오차는 σ_ESA(및 σ_GSA)에 그대로 흡수된다. Kwock 2021은 단일입자는 peak intensity, 필름은 FWHM intensity를 쓴다고 명시 — **강도 정의 차이만으로 σ가 수 배 달라질 수 있음**.

### 3.2 실제 문헌의 degeneracy 증거

- **Dudek et al. 2023 (Yb/Pr)**: W_CR1은 "PA에 핵심처럼 보이지만" 더 강한 w65(ETU)가 ¹G₄를 채우므로 결과에 영향 없음; NR₅₄, NR₄₃은 ±3 decade 변화에도 영향 미미 → **sloppy/비식별 파라미터**의 전형 [verified].
- **Skripka & Chan 2025 (Mater. Horiz.)**: DRE는 "수천 개 rate 상수의 추정치"를 쓰며 대부분 정밀 측정이 어렵다; rate 상수가 도핑·population·파워에 독립이라는 가정은 과단순화; 과-파라미터화 위험. 실무적으로 **pre-avalanche 영역은 극히 어두워 검출 노이즈에 묻히고**, log-log slope 추출은 **배경 차감에 매우 민감**, 인접 점 몇 개로 slope를 구하면 N을 과대추정 [verified].
- **Tessitore et al. Nanoscale 2019 리뷰**: "포함할 과정 수와 피팅할 데이터 수를 균형 있게 — rate equation의 과-파라미터화를 피하라" [verified].
- **Teitelboim et al. JPCC 123, 2678 (2019)**: KMC에서 ⁴S₃/₂를 탈여기시키는 전이가 1700가지 이상, 상위 5개가 53 %만 차지; 주요 전이 4종을 knockout 해도 τ_eff가 크게 버팀(679 → 285 μs) — 네트워크의 **집단적·강건(robust)·이력 의존** 동역학. 측정 수명은 개별 준위의 고유 수명이 아님 [verified]. → **수명을 "고정 파라미터"로 바로 넣는 것도 위험**; 파워 의존 수명은 모델 출력으로 피팅해야 함.
- Lee 2021은 여러 농도/셸 샘플을 σ, s31, Q 공유 + w_nr 샘플별 자유로 두어 degeneracy를 줄임 [verified 코드].

### 3.3 Sloppiness와 profile likelihood

- **Gutenkunst et al., PLoS Comput. Biol. 3, e189 (2007)** [verified 서지]: 다수의 systems biology 모델에서 FIM 고유값이 로그 스케일로 고르게 퍼져 있는 "universal sloppiness" — 개별 파라미터는 불확정이어도 예측은 잘 제약될 수 있음.
- **Transtrum, Machta, Sethna, J. Chem. Phys. 143, 010901 (2015)** [verified 서지] 및 **Transtrum & Qiu, PRL 113, 098701 (2014)** (Manifold Boundary Approximation Method, MBAM) [verified 서지]: sloppy 방향을 따라가 모델 경계(파라미터 → 0/∞ 극한)로 **체계적 모델 축약**. PA 모델에서 "빠른 MPR로 중간 준위 소거", "Q22→0" 같은 축약을 데이터 기반으로 정당화하는 도구.
- **Raue et al., Bioinformatics 25, 1923 (2009)** [verified]: profile likelihood로 structural(평평한 profile) vs practical(한쪽만 열린 profile) 비식별성 구분, 1차원 신뢰구간. **Wieland et al., Curr. Opin. Syst. Biol. 25, 60 (2021)** "On structural and practical identifiability" [서지 부분 verified].

### 3.4 최적화 알고리즘 — 벤치마크 결론

- **Raue et al., PLoS ONE 8, e74335 (2013)** [verified]: 15개 알고리즘 비교에서 **sensitivity equation 기반 gradient + Latin hypercube multistart deterministic local 최적화가 확률적 방법보다 정확도·속도에서 수 자리 우수**; finite difference 도함수는 ~10배 느리고 수치적으로 불안정하여 조기 종료·가짜 local optimum을 유발; **log-scale 파라미터**와 무차원화가 성능 향상; 노이즈 파라미터를 동시 추정하는 것이 전처리보다 나음.
- **Villaverde et al., Bioinformatics 35, 830 (2019)**, DOI 10.1093/bioinformatics/bty736 (요청서의 btz020 아님) [verified]: 36–383 파라미터 7개 벤치마크. 최고 성능은 **eSS(enhanced scatter search) + interior-point local + adjoint gradient 하이브리드**; gradient-based multistart도 경쟁력 있음; PSO 등 순수 확률적은 열세; log 파라미터화가 전 방법에서 유리.
- **Egea et al., BMC Bioinformatics 15, 136 (2014)** MEIGO [verified 서지]: eSS/VNS 구현.
- **CMA-ES** (Hansen & Ostermeier, Evol. Comput. 9, 159 (2001)) [verified 서지], **Differential Evolution** (Storn & Price, J. Glob. Optim. 11, 341 (1997)) [verified 서지]: gradient가 불안정(bistable 분기 점프 등)할 때의 대안. 위 벤치마크에 비추면 **전역 탐색 단계에만 쓰고 반드시 gradient local refinement와 결합**하는 것이 합리적 [derived].
- **Hierarchical optimization**: Weber et al. (IFAC Proc. 2011) [서지 verified], **Loos, Krause, Hasenauer, Bioinformatics 34, 4266 (2018)** [verified]: relative data에서 scale(및 noise) 파라미터를 inner problem에서 **해석적으로** 풀어 outer 차원을 줄임 — Gaussian noise면 가중 최소제곱 closed form, Laplace noise도 처리; 수렴성·계산시간 모두 개선. **Schmiester et al., Bioinformatics 36, 594 (2020)** (btz581) [서지 verified]: offset까지 포함한 대규모 확장.
  - PA 적용 [derived]: 채널 k의 측정 y_k(P_j) = s_k · m_k(P_j; θ) + b_k 에서 (s_k, b_k)는 θ가 주어지면 선형 최소제곱으로 닫힌해. log-space 잔차를 쓰면 s_k는 log-offset이 되어 평균 잔차로 닫힌해. **채널 gain이 미지인 다채널(³H₄ 800 nm, ¹G₄ 475 nm 등) 동시 피팅에 그대로 적용 가능**.
- 소프트웨어: **pyPESTO** (Schälte et al., Bioinformatics 39, btad711 (2023)) [verified 서지] — multistart, hierarchical, profile, sampling 통합; **AMICI**(CVODES 기반 forward/adjoint sensitivity) [서지 unverified].

### 3.5 정상상태 계산: 적분 vs root finding vs continuation

- **ODE 적분으로 정상상태**: Lee 코드 방식. 임계 부근 critical slowing down(rise time 수백 ms)으로 **수렴 불확실 + 매우 비쌈**. 단, rise time 데이터를 피팅할 때는 필요.
- **Root finding** [derived]: 3-level 모델의 정상상태는 n₂, n₃에 대한 **2차 다항식 2개**(Bezout 상한 4개 해)이므로, 제거법(resultant)으로 단변수 4차 다항식으로 환원하거나 Newton + 해석적 Jacobian으로 풀 수 있다. 물리해 조건(0 ≤ n_i ≤ 1, Jacobian 고유값 실부 음수=안정)으로 선별. 일반 N-level Chan형 모델도 질량보존 식 하나를 대체한 Newton 해법이 표준.
- **Numerical continuation** [derived + 도구 문헌]: 파워 P를 continuation 파라미터로 한 pseudo-arclength 추적으로 **fold(turning point) = switch-on/off threshold**를 직접 찾고, 불안정 가지까지 계산. MATCONT (Dhooge, Govaerts, Kuznetsov, ACM TOMS 29 (2003)) [verified 서지]; Python에서는 자체 구현 또는 PyDSTool/BifurcationKit.jl 등 [unverified]. Skripka 2025처럼 bistable한 경우 **실험의 up-sweep ↔ 하위 안정가지, down-sweep ↔ 상위 안정가지**로 매핑해야 하며, 정상상태 root를 "초기값 의존"으로 아무거나 고르면 피팅 목적함수가 불연속이 된다.
- 권장 [derived]: (1) 기본 forward model = continuation으로 P-grid 전체 가지 추적(첫 점은 Newton, 이후 predictor–corrector), (2) 데이터에 sweep 방향 메타데이터 저장, (3) rise/decay 피팅이 필요할 때만 stiff solver(BDF/Radau, 해석 Jacobian)로 적분.

### 3.6 Time-resolved 데이터로 degeneracy 깨기

- Critical slowing down: rise time은 임계 근방에서 최대(Lee 2021: ~608 ms; Skripka & Chan 2025: "bucket 비유" — 내부 rate는 빨라지지만 채워야 할 population이 커져 rise time이 늘어남) [verified].
- Dudek 2022: 크기에 따라 t₅₀ 100 ms(bulk) → 9 ms(core NP) — rise time이 k₂, k₃(표면 손실)에 민감 [verified].
- Bednarkiewicz 2019: rise time과 time-gated detection(3–6 ms)으로 PA 비선형성을 더 키우는 시뮬레이션 [verified].
- [derived] 3.1의 스케일 대칭은 동역학 dn/dt = c·f(n)에서 시간축을 c배 압축하므로, **rise/decay 한 곡선만 있어도 절대 스케일 c가 고정**된다. 또 rise time의 파워 의존(피크 위치 ≈ threshold, 피크 높이 ∝ 1/W₂ 스케일)은 정상상태와 다른 조합에 감도를 가져 profile likelihood를 좁힌다. 플랫폼은 **steady-state + rise + decay 공동 목적함수**를 기본 설계로 할 것.

---

## 4. Mean-field를 넘어서: microscopic / KMC

| 모델 | 방법 | 대상/규모 | 결론 | 출처 |
|---|---|---|---|---|
| Microscopic RE | 격자 위 이온별 RE, 실제 거리로 ET | 571 이온 × 7 준위 (3 997 변수), β-NaYF₄:Tm | average RE는 decay 재현 실패; microscopic은 재현; migration ≈ 0 | Villanueva-Delgado 2016 [verified] |
| KMC (Gillespie, rejection-free) | 8 nm 구형 UCNP, 무작위 도핑; JO radiative + Kushida ET, DRE와 같은 rate 상수 | Yb/Er 20/20, 20/60 % | 파워 의존 수명(20/60 %에서 10²–10⁶ W/cm² 범위 6배 감소) 재현; **고농도(>20 %)에서 DRE 적분이 time-step 축소에도 수렴하지 않음** | Teitelboim 2019 [verified] |
| Microscopic/KMC 기반 설계 | 다중 셸 이종구조; KMC가 >30 nm 다중셸에서 수주–수개월 소요 → GNN surrogate(SUNSET 데이터셋 6 000–30 000 구조) | – | Skripka & Chan 2025 리뷰 [verified]; Nature Comput. Sci. 2025 GNN 논문들 [서지 unverified] |

소형 입자에서의 이온 수 [derived]: β-NaYF₄ 단위셀 부피 = (√3/2)a²c = 109.1 Å³, RE 자리 1.5개/셀 가정 → RE 밀도 ≈ 1.38×10²² cm⁻³. 직경 10 nm 코어 ≈ 7 200 RE 자리 → 8 % Tm ≈ 580개, 1 % Tm ≈ 70개. Skripka 2023 Nano Lett.도 "<500 Tm per NP" 입자를 다룸 [verified]. 이 규모에서:
- 셸 근처(표면 수 nm 이내) 이온 비율이 크고, Kwock 식에서 λ⁻¹ ≈ 1.1 nm 이므로 **표면 소광은 공간적으로 매우 비균질** [derived].
- CR 이차항의 평균장 근사(⟨n₁n₃⟩ ≈ ⟨n₁⟩⟨n₃⟩)는 이온 수가 적고 상관(correlation)이 있을 때 틀릴 수 있음 — Guy/Joubert/Jacquier 1997이 PA에서 바로 이 문제를 다룬 논문 [내용 unverified].
- 단일입자 간 threshold 이질성(Kwock 2021)은 **앙상블 S-curve를 단일 모델로 피팅하면 slope가 무뎌진다**는 뜻 → 앙상블 forward model에 입자간 분포(예: W₂,NR의 lognormal)를 적분하는 옵션 필요 [derived].

**Mean-field 사용 가능 범위(문헌 근거로 본 판단)** [derived]: (i) 중간 농도(대략 ≤20 %), (ii) 정상상태 또는 느린 동역학(ms) 위주, (iii) 표면 손실을 균일 k_surf로 흡수 가능한 두꺼운 셸, (iv) 앙상블이 아니라 단일입자이거나 이질성 분포를 따로 모델링할 때. 이를 벗어나면(고농도, 얇은 셸, μs 이하 decay 형태 피팅) microscopic/KMC로 교차검증 권장.

---

## 5. 플랫폼 설계에 대한 구체적 권고 (문헌 근거 → 설계 결정)

1. **모델 층위**: (L1) Lee/Szalkowski 3-level Tm 모델을 reference로 내장하고 공개 코드 값으로 회귀 테스트(Sample No. 4의 s31 = 10 240 s⁻¹ 등 재현). (L2) Chan형 N-level 자동 생성(JO A_ij + Kushida C_DA + 에너지갭/PAET + detailed balance). (L3) 셸/표면: Kwock형 거리 의존 k_surf 또는 shell-resolved compartments. (L4) 선택적 KMC 교차검증.
2. **파라미터화**: 모든 양의 파라미터는 log10. MPR/PAET는 (B, α, ħω) + ΔE로 생성, rate를 개별 자유 파라미터로 두지 않음. 농도 의존은 함수형 고정 + 계수 공유(여러 농도 샘플 global fit). 샘플별 자유 파라미터는 비방사 손실(k₂, k₃ 또는 Γ₀, d_shell)로 한정(Lee 2021 관행).
3. **관측 모델**: y_ch = s_ch · A_rad · n_emit(P) (+ b_ch); (s_ch, b_ch)는 hierarchical 해석해(Loos 2018). 앙상블이면 Gaussian beam 적분 ∫n(I)/I dI, 단일입자면 peak intensity; 강도 정의(peak/FWHM/평균)를 메타데이터로 강제.
4. **정상상태 solver**: Newton + 해석 Jacobian + pseudo-arclength continuation, 안정성 판정, sweep 방향별 가지 선택. ODE 적분(stiff, BDF/Radau)은 동역학 피팅용.
5. **목적함수**: log-intensity 잔차(배경 차감 오차 모델 포함; Skripka & Chan 2025 경고) + rise time(95 % 또는 50 %) + decay; 노이즈 파라미터 동시 추정(Raue 2013).
6. **최적화**: Latin hypercube multistart(수백 회) + gradient local(forward sensitivity, 파라미터 수 ≤ ~30이면 forward로 충분 [derived]) → 필요 시 eSS/CMA-ES 전역 단계. 결과는 waterfall plot(시작점별 최종 목적값 정렬)로 수렴 진단(Raue 2013 관행).
7. **식별성 진단**: 매 피팅 후 profile likelihood; 평평한 profile은 해당 파라미터 고정 또는 MBAM식 축약. 구조적 스케일 대칭(3.1)은 사전에 하나의 절대 rate 고정으로 제거.
8. **열 효과**: 선택적 온도 ODE(C dT/dt = P_abs·η_heat − G(T − T₀))와 Szalkowski형 온도의존 rate 결합 — 실험의 up/down hysteresis가 있을 때만 활성화하고, Skripka 2025처럼 비열적 bistability도 가능함을 분리 판정.

---

## 6. 검증된 참고문헌 목록 (Bibliography)

상태: [V] = DOI/서지 원문·출판사 페이지·Crossref·공식 인덱스로 확인; [V-partial] = 제목/학술지 확인, 권·쪽 일부 미확인; [U] = 미확인.

### PA 모델 및 실험
1. Lee, C.; Xu, E. Z.; Liu, Y.; Teitelboim, A.; Yao, K.; Fernandez-Bravo, A.; Kotulska, A. M.; Nam, S. H.; Suh, Y. D.; Bednarkiewicz, A.; Cohen, B. E.; Chan, E. M.; Schuck, P. J. "Giant nonlinear optical responses from photon-avalanching nanoparticles." *Nature* **589**, 230–235 (2021). DOI 10.1038/s41586-020-03092-9 [V]. 코드: github.com/nawhgnahc/Photon_Avalanche_DRE_calculation [V].
2. Szalkowski, M.; Dudek, M.; Korczak, Z.; Lee, C.; Marciniak, Ł.; Chan, E. M.; Schuck, P. J.; Bednarkiewicz, A. "Predicting the impact of temperature dependent multi-phonon relaxation processes on the photon avalanche behavior in Tm³⁺:NaYF₄ nanoparticles." *Opt. Mater. X* **12**, 100102 (2021). DOI 10.1016/j.omx.2021.100102 [V].
3. Bednarkiewicz, A.; Chan, E. M.; Kotulska, A.; Marciniak, L.; Prorok, K. "Photon avalanche in lanthanide doped nanoparticles for biomedical applications: super-resolution imaging." *Nanoscale Horiz.* **4**, 881 (2019). DOI 10.1039/c9nh00089e [V].
4. Dudek, M.; Szalkowski, M.; Misiak, M.; Ćwierzona, M.; Skripka, A.; Korczak, Z.; Piątkowski, D.; Woźniak, P.; Lisiecki, R.; Goldner, P.; Maćkowski, S.; Chan, E. M.; Schuck, P. J.; Bednarkiewicz, A. "Size-dependent photon avalanching in Tm³⁺ doped LiYF₄ nano, micro, and bulk crystals." *Adv. Opt. Mater.* **10**, 2201052 (2022). DOI 10.1002/adom.202201052 [V].
5. Kwock, K. W. C.; Lee, C.; Teitelboim, A.; Liu, Y.; Yao, K.; Alam, S. B.; Cohen, B. E.; Chan, E. M.; Schuck, P. J. "Surface-sensitive photon avalanche behavior revealed by single-avalanching-nanoparticle imaging." *J. Phys. Chem. C* **125** (2021). DOI 10.1021/acs.jpcc.1c07721 [V-partial: 쪽 미확인].
6. Liang, Y.; Zhu, Z.; Qiao, S.; et al. "Migrating photon avalanche in different emitters at the nanoscale enables 46th-order optical nonlinearity." *Nat. Nanotechnol.* **17**, 524–530 (2022). DOI 10.1038/s41565-022-01101-8 [V]. (요청서의 "Nature 2022"가 아니라 Nat. Nanotechnol.)
7. Skripka, A.; Zhang, Z.; Qi, X.; Ursprung, B.; Ercius, P.; Cohen, B. E.; Schuck, P. J.; Jaque, D.; Chan, E. M. "Intrinsic optical bistability of photon avalanching nanocrystals." *Nat. Photon.* **19**, 212–218 (2025). DOI 10.1038/s41566-024-01577-x [V]; arXiv:2403.04098 [V].
8. Skripka, A.; et al. "A generalized approach to photon avalanche upconversion in luminescent nanocrystals." *Nano Lett.* **23**, 7100 (2023). DOI 10.1021/acs.nanolett.3c01955 [V]; arXiv:2306.06200.
9. Skripka, A.; Chan, E. M. "Unraveling the myths and mysteries of photon avalanching nanoparticles." *Mater. Horiz.* **12**, 3575 (2025). DOI 10.1039/d4mh01798f [V].
10. Dudek, M.; Korczak, Z.; Prorok, K.; Bezkrovnyi, O.; Sun, L.; Szalkowski, M.; Bednarkiewicz, A. "Understanding Yb³⁺-sensitized photon avalanche in Pr³⁺ co-doped nanocrystals: modelling and optimization." *Nanoscale* **15**, 18613 (2023). DOI 10.1039/d3nr04409b [V]. 데이터: Zenodo 10654092 [V].
11. Aggarwal, S.; et al. "Recent advances in fundamental research on photon avalanches on the nanometre scale." *Nanoscale* **17**, 6329 (2025). DOI 10.1039/d4nr03493g [V-partial] (리뷰; 일부 서술 오류 관찰 — 예: Nd ground를 ⁴I₁₅/₂로 표기 — 2차 인용 시 주의).
12. Korczak, Z.; et al. "Sensitized photon avalanche nanothermometry in Pr³⁺ and Yb³⁺ co-doped NaYF₄ colloidal nanoparticles." *Low Temp. Phys.* **49**, 322 (2023) [V-partial, DOI 미확인].
13. Bednarkiewicz, A.; et al. "All-optical data processing with photon-avalanching nanocrystalline photonic synapse." *Adv. Mater.* (2023). DOI 10.1002/adma.202304390 [V-partial].
14. Lee, C.; et al. "Indefinite and bidirectional near-infrared nanocrystal photoswitching." *Nature* (2023). DOI 10.1038/s41586-023-06076-7 [V-partial; 권·쪽 미확인].
15. Goldner, Ph.; Pellé, F. "Photon avalanche fluorescence and lasers." *Opt. Mater.* **5**, 239–249 (1996). DOI 10.1016/0925-3467(96)00003-1 [V].
16. Joubert, M.-F. "Photon avalanche upconversion in rare earth laser materials." *Opt. Mater.* **11**, 181–203 (1999). DOI 10.1016/S0925-3467(98)00043-3 [V].
17. Guy, S.; Joubert, M.-F.; Jacquier, B. "Photon avalanche and the mean-field approximation." *Phys. Rev. B* **55**, 8240 (1997). DOI 10.1103/PhysRevB.55.8240 [V-partial: 제목·권·쪽 확인, 저자 이니셜·내용 미확인].
18. Pollnau, M.; Gamelin, D. R.; Lüthi, S. R.; Güdel, H. U.; Hehlen, M. P. "Power dependence of upconversion luminescence in lanthanide and transition-metal-ion systems." *Phys. Rev. B* **61**, 3337–3346 (2000). DOI 10.1103/PhysRevB.61.3337 [V].
19. Chivian, J. S.; Case, W. E.; Eden, D. D. "The photon avalanche: a new phenomenon in Pr³⁺-based infrared quantum counters." *Appl. Phys. Lett.* **35**, 124 (1979) [V-partial; DOI 10.1063/1.91044 U].
20. Hehlen, M. P.; Güdel, H. U.; et al. "Cooperative bistability in dense, excited atomic systems." *Phys. Rev. Lett.* **73**, 1103 (1994). DOI 10.1103/PhysRevLett.73.1103 [V-partial].
21. Gamelin, D. R.; Lüthi, S. R.; Güdel, H. U. "The role of laser heating in the intrinsic optical bistability of Yb³⁺-doped bromide lattices." *J. Phys. Chem. B* **104**, 11045 (2000) [V-partial; DOI U].
22. "Kolesov" 관련 PA 리뷰 — 이번 조사에서 해당 문헌을 특정하지 못함 [U].

### UCNP DRE, 분광 파라미터, 표면 소광
23. Villanueva-Delgado, P.; Krämer, K. W.; Valiente, R. "Simulating energy transfer and upconversion in β-NaYF₄:Yb³⁺,Tm³⁺." *J. Phys. Chem. C* **119**, 23648–23657 (2015). DOI 10.1021/acs.jpcc.5b06770 [V].
24. Villanueva-Delgado, P.; Krämer, K. W.; Valiente, R.; de Jong, M.; Meijerink, A. "Modeling blue to UV upconversion in β-NaYF₄:Tm³⁺." *Phys. Chem. Chem. Phys.* **18**, 27396 (2016). DOI 10.1039/c6cp04347j [V].
25. Teitelboim, A.; Tian, B.; Garfield, D. J.; Fernandez-Bravo, A.; Gotlin, A. C.; Schuck, P. J.; Cohen, B. E.; Chan, E. M. "Energy transfer networks within upconverting nanoparticles are complex systems with collective, robust, and history-dependent dynamics." *J. Phys. Chem. C* **123**, 2678 (2019). DOI 10.1021/acs.jpcc.9b00161 [V].
26. Chan, E. M.; Gargas, D. J.; Schuck, P. J.; Milliron, D. J. "Concentrating and recycling energy in lanthanide codopants for efficient and spectrally pure emission: the case of NaYF₄:Er³⁺/Tm³⁺ upconverting nanocrystals." *J. Phys. Chem. B* (2012). DOI 10.1021/jp302401j [V-partial].
27. Chan, E. M. "Combinatorial approaches for developing upconverting nanomaterials: high-throughput screening, modeling, and applications." *Chem. Soc. Rev.* (2015). DOI 10.1039/C4CS00205A [V-partial]; ESI 식 S1–S10 [V].
28. Tessitore, G.; Mandl, G. A.; Brik, M. G.; Park, W.; Capobianco, J. A. "Recent insights into upconverting nanoparticles: spectroscopy, modeling, and routes to improved luminescence." *Nanoscale* (2019). DOI 10.1039/C9NR02291K [V-partial].
29. Johnson, N. J. J.; He, S.; Diao, S.; Chan, E. M.; Dai, H.; Almutairi, A. "Direct evidence for coupled surface and concentration quenching dynamics in lanthanide-doped nanocrystals." *J. Am. Chem. Soc.* **139**, 3275 (2017). DOI 10.1021/jacs.7b00223 [V-partial: DOI·제목 확인, 저자·쪽은 기억 기반].
30. Fischer, S.; Bronstein, N. D.; Swabeck, J. K.; Chan, E. M.; Alivisatos, A. P. "Precise tuning of surface quenching for luminescence enhancement in core–shell lanthanide-doped nanocrystals." *Nano Lett.* **16**, 7241–7247 (2016). DOI 10.1021/acs.nanolett.6b03683 [V].
31. Rabouw, F. T.; et al. "Quenching pathways in NaYF₄:Er³⁺,Yb³⁺ upconversion nanocrystals." *ACS Nano* (2018). DOI 10.1021/acsnano.8b01545 [V-partial].
32. Wang, Z.; Meijerink, A. "Concentration quenching in upconversion nanocrystals." *J. Phys. Chem. C* **122**, 26298–26306 (2018). DOI 10.1021/acs.jpcc.8b09371 [V].
33. Hossan, M. Y.; et al. "Explaining the nanoscale effect in the upconversion dynamics of β-NaYF₄:Yb³⁺,Er³⁺ core and core–shell nanocrystals." *J. Phys. Chem. C* **121**, 16592 (2017). DOI 10.1021/acs.jpcc.7b04567 [V-partial].
34. Anderson, R. B.; Smith, S. J.; May, P. S.; Berry, M. T. "Revisiting the NIR-to-visible upconversion mechanism in β-NaYF₄:Yb³⁺,Er³⁺." *J. Phys. Chem. Lett.* **5**, 36 (2014) [V-partial; DOI U].
35. Kaiser, M.; Würth, C.; Kraft, M.; Hyppänen, I.; Soukka, T.; Resch-Genger, U. "Power-dependent upconversion quantum yield of NaYF₄:Yb³⁺,Er³⁺ nano- and micrometer-sized particles – measurements and simulations." *Nanoscale* **9**, 10051 (2017). DOI 10.1039/C7NR02449E [V-partial].
36. Kaiser, M.; Würth, C.; Kraft, M.; Soukka, T.; Resch-Genger, U. "Explaining the influence of dopant concentration and excitation power density on the luminescence and brightness of β-NaYF₄:Yb³⁺,Er³⁺ nanoparticles: measurements and simulations." *Nano Res.* **12**, 1871–1879 (2019). DOI 10.1007/s12274-019-2450-4 [V].
37. Würth, C.; Kaiser, M.; Wilhelm, S.; et al. "Excitation power dependent population pathways and absolute quantum yields of upconversion nanoparticles in different solvents." *Nanoscale* **9**, 4283–4294 (2017) [V-partial].
38. Hudry, D.; et al. "Direct evidence of significant cation intermixing in upconverting core@shell nanocrystals." *Chem. Mater.* **29**, 9238 (2017) [V-partial; DOI U].

### 물리 제약 (MPR, ET, 분광)
39. van Dijk, J. M. F.; Schuurmans, M. F. H. "On the nonradiative and radiative decay rates and a modified exponential energy gap law for 4f–4f transitions in rare-earth ions." *J. Chem. Phys.* **78**, 5317–5323 (1983). DOI 10.1063/1.445485 [V, Crossref].
40. Riseberg, L. A.; Moos, H. W. "Multiphonon orbit-lattice relaxation of excited states of rare-earth ions in crystals." *Phys. Rev.* **174**, 429 (1968). DOI 10.1103/PhysRev.174.429 [V].
41. Miyakawa, T.; Dexter, D. L. "Phonon sidebands, multiphonon relaxation of excited states, and phonon-assisted energy transfer between ions in solids." *Phys. Rev. B* **1**, 2961 (1970). DOI 10.1103/PhysRevB.1.2961 [V].
42. Inokuti, M.; Hirayama, F. "Influence of energy transfer by the exchange mechanism on donor luminescence." *J. Chem. Phys.* **43**, 1978 (1965). DOI 10.1063/1.1697063 [V].
43. McCumber, D. E. "Einstein relations connecting broadband emission and absorption spectra." *Phys. Rev.* **136**, A954 (1964). DOI 10.1103/PhysRev.136.A954 [V].
44. Judd, B. R. *Phys. Rev.* **127**, 750 (1962); Ofelt, G. S. *J. Chem. Phys.* **37**, 511 (1962) [U — 이번 조사에서 미확인, 표준 인용].
45. Tm³⁺:NaYF₄/LiYF₄의 JO Ω_λ 및 A_ij 표 — 이번 조사에서 원문 표를 확보하지 못함 [U]. 대안: Villanueva-Delgado 2016의 실측 수명·분기비(위 2.3) 사용, 또는 β-NaGdF₄:Tm 단결정 JO 분석(Villanueva-Delgado et al., J. Lumin. 2016) [U].

### 식별성·최적화
46. Gutenkunst, R. N.; Waterfall, J. J.; Casey, F. P.; Brown, K. S.; Myers, C. R.; Sethna, J. P. "Universally sloppy parameter sensitivities in systems biology models." *PLoS Comput. Biol.* **3**, e189 (2007). DOI 10.1371/journal.pcbi.0030189 [V].
47. Raue, A.; Kreutz, C.; Maiwald, T.; Bachmann, J.; Schilling, M.; Klingmüller, U.; Timmer, J. "Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood." *Bioinformatics* **25**, 1923 (2009). DOI 10.1093/bioinformatics/btp358 [V].
48. Raue, A.; Schilling, M.; Bachmann, J.; Matteson, A.; Schelke, M.; Kaschek, D.; Hug, S.; Kreutz, C.; Harms, B. D.; Theis, F. J.; Klingmüller, U.; Timmer, J. "Lessons learned from quantitative dynamical modeling in systems biology." *PLoS ONE* **8**, e74335 (2013). DOI 10.1371/journal.pone.0074335 [V].
49. Villaverde, A. F.; Fröhlich, F.; Weindl, D.; Hasenauer, J.; Banga, J. R. "Benchmarking optimization methods for parameter estimation in large kinetic models." *Bioinformatics* **35**, 830–838 (2019). DOI 10.1093/bioinformatics/bty736 [V].
50. Weber, P.; Hasenauer, J.; Allgöwer, F.; Radde, N. "Parameter estimation and identifiability of biological networks using relative data." *IFAC Proc. Vol.* **44** (2011) [V-partial].
51. Loos, C.; Krause, S.; Hasenauer, J. "Hierarchical optimization for the efficient parametrization of ODE models." *Bioinformatics* **34**, 4266–4273 (2018). DOI 10.1093/bioinformatics/bty514 [V].
52. Schmiester, L.; Weindl, D.; Hasenauer, J. "Efficient parameterization of large-scale dynamic models based on relative measurements." *Bioinformatics* **36**, 594 (2020). DOI 10.1093/bioinformatics/btz581 [V-partial].
53. Transtrum, M. K.; Machta, B. B.; Brown, K. S.; Daniels, B. C.; Myers, C. R.; Sethna, J. P. "Perspective: Sloppiness and emergent theories in physics, biology, and beyond." *J. Chem. Phys.* **143**, 010901 (2015). DOI 10.1063/1.4923066 [V-partial: 저자 목록 기억 기반].
54. Transtrum, M. K.; Qiu, P. "Model reduction by manifold boundaries." *Phys. Rev. Lett.* **113**, 098701 (2014). DOI 10.1103/PhysRevLett.113.098701 [V].
55. Egea, J. A.; Henriques, D.; Cokelaer, T.; Villaverde, A. F.; MacNamara, A.; Danciu, D.-P.; Banga, J. R.; Saez-Rodriguez, J. "MEIGO: an open-source software suite based on metaheuristics for global optimization in systems biology and bioinformatics." *BMC Bioinformatics* **15**, 136 (2014). DOI 10.1186/1471-2105-15-136 [V-partial: 저자 목록 기억 기반].
56. Hansen, N.; Ostermeier, A. "Completely derandomized self-adaptation in evolution strategies." *Evol. Comput.* **9**, 159 (2001). DOI 10.1162/106365601750190398 [V].
57. Storn, R.; Price, K. "Differential evolution – a simple and efficient heuristic for global optimization over continuous spaces." *J. Glob. Optim.* **11**, 341–359 (1997). DOI 10.1023/A:1008202821328 [V].
58. Dhooge, A.; Govaerts, W.; Kuznetsov, Yu. A. "MATCONT: a MATLAB package for numerical bifurcation analysis of ODEs." *ACM Trans. Math. Softw.* **29**, 141 (2003). DOI 10.1145/779359.779362 [V-partial: 쪽 미확인].
59. Schälte, Y.; et al. "pyPESTO: a modular and scalable tool for parameter estimation for dynamic models." *Bioinformatics* **39**, btad711 (2023) [V-partial].
60. Wieland, F.-G.; et al. "On structural and practical identifiability." *Curr. Opin. Syst. Biol.* **25**, 60 (2021) [V-partial].

### 조사 한계
- 접근 실패(403/robots/rate-limit): APS 원문(Guy 1997, Pollnau 2000 본문), ACS 원문(Johnson 2017, Hossan 2017, Teitelboim 2019 본문 일부), PMC(에너지갭 법칙 온도의존 리뷰), OSTI 전반, Crossref 대량 조회. 따라서 **Yb/Er 표준 파라미터 표, Tm³⁺ JO Ω_λ·A_ij 표, 호스트별 에너지갭 법칙 상수(B, α)는 원문 수치로 검증하지 못했다.** 해당 값은 위에 [unverified]로 표시했으며, 플랫폼 기본값으로 쓰기 전 원문(PDF) 대조가 필요하다.
- Lee 2021 SI Tables S3–S5 원문은 확보하지 못했고, 수치는 공개 코드와 Szalkowski 2021 재현표로 대체 검증했다.
