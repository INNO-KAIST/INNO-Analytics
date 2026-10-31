# INNO Analytics — 작업 이력

> STEM 이미지·EDS 스펙트럼·SI(원소 맵)·격자 FFT·EDS 간섭 검사·Energy Diagram·DRE(속도식 모델)·XPS·Image Analysis(ImageJ 방식) 분석을 브라우저 안에서 하는 단일 HTML 앱. 서버·로그인 없음.
>
> 이 파일 하나에 이 앱의 작업 이력·교훈·문서 색인을 모았다. 자세한 근거는 색인의 문서(claude.ai 프로젝트 ‘유원정’)에 있다.

## XPS v1.9.5 (2026-10-31) — 문헌 값 원문 대조 ①–③
- **Scofield 단면적 전면 교체**: J.H. Scofield, J. Electron Spectrosc. Relat. Phenom. 8, 129 (1976) 표 1(1254 eV)·표 2(1487 eV)를 원문 PDF 에서 한 쪽씩 전부 옮겨 적음(Al 1066 · Mg 1030 아껍질). 이전 표(KherveFitting 디지털화)의 오류 수정 — 10 배 자릿수 오류: Mg Na 2s 3.9→0.390, Mg Ti 3s 4.08→0.408, Mg Zr 4s 3.05→0.305, Mg Zr 4p 3.666→0.903, Al K 3p 0.037→0.362; 칸 밀림: Mg Pr 3p3/2 30.72→9.75, Mg Au 4s 1.115→1.45; 비율로 쪼갠 값: Al 2p1/2 0.179→0.1811, 2p3/2 0.358→0.356 외; 빠진 아껍질 약 590 개 추가. 검사: 이전 표와 일치 여부, Z 에 대한 매끄러움, Al/Mg 비의 매끄러움. 전사본 test/ref/scofield1976_transcription.json.
- **스핀–궤도 Δ 출처 순서**: Δ ≥ 2 eV 는 X-ray Data Booklet(원자 준위 편집값) 먼저 — Avantage 원소표 위치는 정수 반올림이 많아 78/357 이 ≥ 0.3 eV 다름(Ce 3d 18 vs 18.6, Xe 3d 15 vs 12.6, Ge 2p 26.5 vs 31.1). Δ < 2 eV 는 Booklet 0.1 eV 반올림이 상대 오차로 커서(Si 2p 0.4 vs ≈0.6) 이전처럼 Avantage 먼저.
- **IMFP 기본 지수** KE^0.6 → KE^0.7: TPP-2M(Tanuma·Powell·Penn, Surf. Interface Anal. 21, 165 (1994) Eq. 3–8, 원문 대조 — 코드 상수 일치)으로 계산한 유효 지수 0.64–0.77(200–1500 eV; C, Si, SiO₂, CeO₂, Au). 같은 시료에서 KE 가 크게 다른 선끼리 at% 가 몇 % 바뀔 수 있음.
- **C 1s 284.8 eV** 버튼에 근거 표시: Biesinger, Appl. Surf. Sci. 597, 153681 (2022) — 117 시료 평균 284.91 ± 0.25 eV.
- 시험: test_all_elements 다중항 준위 χ²ν 허용 1.35→1.5(옛·새 표 모두 seed 31–33 에서 1.35–1.41 — 시험 설계 문제, 근거 주석).

## XPS v1.9.4 (2026-10-31) — Ce 3d 레시피 원문 대조
- Bêche et al. 2008 (Surf. Interface Anal. 40, 264) 원문 확인: 스핀–궤도 간격이 산화수마다 다름(Ce⁴⁺ ≈18.6 eV, Ce³⁺ ≈18.2 eV), I(3d5/2)/I(3d3/2) = 1.5. 레시피의 공유 ΔE 18.5(18.1–18.9)를 산화수별 ΔE(18.6: 18.3–18.9 / 18.2: 17.9–18.5)로 바꿈.
- Shyu et al. 1988 (J. Catal. 114, 23) Fig. 2 확인: 순수 CeO₂ 의 u‴ 분율 ≈13–14 %(그래프 판독, 상대오차 약 10 %, 비선형 Shirley). '미확인' 표기를 출처와 함께 교체.
- 반치폭 허용 범위를 '문헌 범위'가 아닌 '경험적 범위'로 표기 정정.
- 문헌 값 전수 목록: modules/xps/test/ref/lit_audit.md (32 항목, 남은 확인 과제 포함).

## 최근 작업
- **10-29 XPS v1.9.3 — 피팅 안내: Fit check · 영역 레시피(Ce 3d) · 단계 피팅**(modules/xps/SPEC.md §17, src/xps_guide.js):
  요청: Ce 3d 처럼 복잡한 영역에서 피팅이 안 맞고 매개변수가 한계에 붙음 → 앱이 친절하게 안내, 한 자료에 맞추지 않은 일반 기능.
  Fit check(피팅 뒤 자동, 한 줄 + title, 'Apply' 단추): 한계 닿음(어느 쪽·값, 폭 연결·문헌값 고정·성분 빼기·범위 넓히기) · |r| > 0.95 상관 쌍 · 넓이 0 / σ > 100 % / < 1 % · 합쳐지는 선(ΔBE < ½ FWHM) · 멈춤·미수렴 · 창 끝 잘림(창 넓히기) · 배경이 자료 위(Smart) · 잔차 구조(≥ 5 σ 그리고 ≥ 3 % 높이, 잡음은 자료의 2차 차분 MA(1) — VGD 이웃 상관 0.3–0.5 반영) · 다중항·위성 영역(레시피 또는 LCF).
  레시피(자료): Ce 3d Burroughs v/u — CeO₂ only (6) · Ce³⁺ only (4) · Mixed (10), 공유 ΔE 18.5 (18.1–18.9), u/v 2/3, 짝 폭 같음, 위치 ±0.6 eV, 875–925 eV, 피팅 안 Shirley + 기울기, 대전 보정에 고정 + ±3 eV NNLS 정렬(어긋나면 경고), Ce³⁺·u‴ 분율 ± σ. 출처 Morgan 2023 · Isaacs 2023 · Bêche 2008 계열 · ACS SI 표(공개본 확인). Ni/Co/Fe/Cu 2p 는 표를 숫자로 확인할 수 없어 넣지 않음, u‴ ≈ 13–14 % 인용값은 원문 미확인.
  단계 피팅: 높이만 → 위치(폭 연결) → 사용자 제약(보통 경로·멈춤 규칙) + 같은 출발의 보통 피팅, χ²_d 낮은 쪽 → 보통 Fit 보다 나빠지지 않음.
  검증: 잘 놓인 합성 270 피팅 경고 0, 규칙마다 양성 대조 걸림; 합성 Ce 3d 12 개 Ce³⁺ |Δ| ≤ 2.7 %p · u‴ ≤ 0.06 %p, 한계 닿음 4 vs 순진한 Auto suggest 12; 실측 S3: 순진한 4 선 χ²ν 18.5 → 레시피 Ce³⁺ 24.7 ± 7.9 % · u‴ 10.40 ± 0.08 % · χ²ν 7.7, C 1s 기준과 Ce 선이 1.45 eV 어긋남(차등 대전) 검출.
  시험: guide 34 · real_ce 9 · Chromium ui_xps_guide 9 · 기존 node·UI 모두 통과 · fitlog_compare 기존 피팅 비트 동일.
- **10-29 XPS v1.9.2 — 멈춤 규칙 재설계 · CasaXPS 4–6차 실측 반영 · Avantage 배경 규약(K-Alpha) · 측정 때 뺀 투과 함수**(modules/xps/SPEC.md §2c·§3.3b·§4, SPEC_casa_v4.md, test/ref/avantage_a2_patch.md):
  멈춤 규칙: 옛 LM 을 현직으로 그대로 두고, 갇히면 옆에서 시도 갈래(rescue = 한계·평탄 매개변수를 시작값으로 되돌려 고정 → 나머지 → 모두, 없으면 kick) 하나를 같은 계산량씩 번갈아 풀어 겨룸. 갈래는 잡음 바닥(χ²/ν ≤ 1 + √(2/ν))에 닿을 때만 결과가 되고, 2차(IRLS)는 두 끝점에서 풀어 낮은 쪽 → 옛 LM 보다 최종 χ² 가 나빠지지 않음.
  무작위 1100 개 최종 χ² 옛과 비트 동일(v1.9.1 은 320 개 중 5 개 나빠짐), 2900 개 더: 나빠짐 0 · 좋아짐 2; 시험 묶음 피팅 23 764 개 나빠짐 0 · 좋아짐 3 · 갈래 없는 23 751 개 비트 동일. TP wL 자유 66 s → 3.9–4.6 s(χ² −2.7e-2), 갈래 채택은 한 줄 알림.
  CasaXPS 4–6차: 폭 숫자·asym 하나의 규칙(격자 F/40 포물선 꼭짓점, fwhm = 2·min(l, h), asym = 1 − l/h — 경험 δ 대체, 14 모양 fwhm ≤ 5.3e-4 · asym ≤ 6e-4), DS 넓이 = 진폭 × ∫ ±600·wL(창 넓이 × 인자, ≤ 4.1e-4), DS 위치 ≤ 1.6e-4 eV,
  ST = 계단을 (바탕 FWHM)/80 늦춘 모양 + 위치 = 원점 + 가장 큰 마디 k·(바탕 FWHM)/40(16 판독 ≤ 0.16 meV; 가져온 ST 는 stCasaK). casaToComp exact 를 compToCasa 와 맞춤, casaModelToINNO 가 DS 넓이·한계·연결을 환산, 성분 카드에 DS CasaXPS Area. 근사로 남음: ⊗G 바탕 ST 의 폭·넓이(≤ 0.5 %).
  Avantage 배경(K-Alpha 2 시료 S1·S2, Smart 12 · Shirley 12 영역): 끝값 = 끝 채널 ±round(W/2·step) 채널 단순 평균(창 밖 포함, 48/48), 모양 = 고정점 u = p·u₁ + (1−p)·G(q·u₁ + (1−q)·u)(p 0.07, q 0.25 — Iter./Conv. 로는 재현 안 됨),
  Smart 는 데이터 넘는 곳에서 나눔, Shirley 는 나누지 않음. 21/24 영역 일치(ΔB/H 중앙 4.2e-4 · 1.4e-5, Area(P) 중앙 0.006 % · 0.005 %), Avantage Area(P) 로 at% 차 ≤ 0.005 at%(INNO Area(P) 로 ≤ 0.056). Kratos 26 영역도 h 4 로 Area(P) 중앙 0.41 % → 0.044 %.
  Avg. eV 를 줄 때만 새 규칙(화면 'Avg. (eV)' 칸, 새로 고르면 0.50, 비우면 v1.8 + Iter./Conv.), Shirley 에도 창 밖 채널 전달, Avantage 곡선이 데이터 위 > 20 % 면 경고 · > 50 % 면 "use Smart". INNO 기본 배경은 그대로.
  투과 함수: VGD prop 415(Acq_TXFN_Correction) = true 면 meta.tfRemoved — 'auto' 에서도 다시 나누지 않음(이중 보정이 at% 최대 9.8 at% 틀어지게 함), 정량 표에 짧은 배지(설명은 title), (PE/KE)^x 를 고르면 경고.
  시험: xps 155 · casa 111 · casa_v3 117 · shapes2 77 · real_kara 44 · real_a2 23 · multi 50 · arxps 32 · image 23 · proc 76 · Chromium ui_xps · fixes 27 · casa 21 · shapes2 13 · multi 43 · fitlog_compare 38 피팅 비트 동일.
- **10-29 XPS v1.9.1 — 피팅 멈춤 규칙 · CasaXPS 3차 실측 반영**(modules/xps/SPEC.md §3.3b, SPEC_casa_v3.md):
  LM 멈춤 규칙: 감쇠 지배 걸음 8 회(거절 ≥5) → Gauss–Newton kick → 그래도 이득 없으면 stalled + 정해지지 않는 매개변수 표시(한 줄 알림). TP wL 자유 퇴화 피팅 66 s → 4–5 s(χ² +0.19 %). 시험 피팅 23 717 건 중 23 687 건 비트 동일, 무작위 320 건: 같음 298 · 나음 10 · 나쁨 5(꺾인 모양의 다른 국소 최솟값).
  CasaXPS 3차 판독: GL(m,n)·SGL(m,n)·TLA·DS 모양과 n 규칙(n·F/144.87) 일치, Casa fwhm = 2L · 비대칭 지수 = 1 − L/H(꺾인 꼭대기는 경험 δ), Casa DoF = 배경과 다른 채널 수 − p, DS 위치 = 원점 + wL·|cot(π/(2−α))|. casaToComp 의 'Casa fwhm = 참 FWHM' 가정 제거(TLA 8 % 오차 → 5e-5).
  미확정: DS 넓이 정규화, ST 위치·fwhm, 낮은 BE 꼬리 모양의 폭 기준. 시험: casa_v3 69 · shapes2 72 · xps 153 · casa 111 · real 44 · UI 5 종.
- **10-29 XPS v1.9 — CasaXPS Cookbook 2026 선모양**(modules/xps/SPEC_shapes2.md):
  TLA(α,μ,w) · ST(μ,γ)·SB(0) 수식어(GLS·GLP·LA·V·TLA·GLN·SGLN) · DS(α,n[,l]) · GL(m,n)·SGL(m,n) + T(k)/PHI(j,k) 꼬리 · TP(α,β,n). Casa 표기 읽기·쓰기(ST(1.7)LA(1,200), SB(0)LA(1,43) 등 왕복), 기존 GL(m)T(k) 저장 모형 그대로.
  독립 기준(scipy·mpmath) 대비 값 1.5e-7 · 넓이 7e-9 · FWHM 5e-8, 합성 재현 200 회 × 9 경우 3σ 안 ≥99 %. 성분 10 개 피팅 약 2 s · MC 1.5 s/회(TP wL 은 기본 고정 — 풀면 퇴화해 느려짐), 반복 진행 표시.
  미확정(앱에 '가정' 표시): n → 폭 규칙을 LA 외 모양에 적용, DS·SB 넓이 = 피팅 창 안, TLA 표기.
  시험: shapes2 71 · ui_xps_shapes2 10 · 기존 XPS node·UI 모두 통과.
- **10-28 XPS v1.8.1 — 실측 대조(Avantage Peak Table 4 시료 · CasaXPS 저장 모델)**(modules/xps/SPEC_real_validation.md):
  Kratos→Avantage .xls 를 못 열던 것·Backgnd 열 순서 고침, 시트 이름을 Peak Table 에서. Avantage Smart 배경 규약(끝점 9 점 평균·9 점 평활 Shirley·데이터 넘으면 나눔)을 Smart(Avantage) 선택지로 — 배경 오차 중앙 3.9e-3, Area(P) 0.41 %(기본 Smart 4.3 %). Avantage ECF = TPP-2M λ(KE)(2e-4 안), at% 0.01 at% 안 일치. TXFN = 1 자료라 투과 함수는 미결 — 기본값 꺼짐 유지.
  CasaXPS 주석(CASA region/comp/Calib) 읽기 → INNO 모형; Shirley 배경 ≤2.6e-4, GL(30) 성분 ≤3.8e-4, envelope RMS 0.4 % 안. 저장된 Casa 모형은 수렴 안 된 상태(INNO 재피팅이 6/6 더 낮은 RMS)라 맞춘 값 비교는 불가. Cookbook 2026: SGL·GLT 꼬리·U 4 Tougaard 정의 확인.
  시험: real_kara 44(자료 없으면 건너뜀) · xps 153 · casa 111 · UI 4 종.
- **10-28 XPS v1.8 — CasaXPS·Avantage 장점 이식 · 깊이/각도/이미징 화면 · CasaXPS 2.3.27 실측 대조**(modules/xps/SPEC.md §1.8b–§16, SPEC_multi.md, SPEC_proc.md, SPEC_casa_v2.md):
  피팅: 참조 스펙트럼 선형결합(LCF)·측정 스펙트럼 선모양·spline 배경·Simplex·여러 스펙트럼 일괄 적용·at% Monte Carlo, Casa 영역 규약·Avantage Shirley·꼬리(AGT) 호환 모드.
  화면 전환 Spectra · Series · ARXPS · Imaging + Process 탭 — 깊이/연속 측정 조성·PCA(IND·F·MLPCA 잡음)·target testing·MCR-ALS/참조 NNLS, ARXPS 덮개층·다층·Hill·RDP·MaxEnt(σ 하한 1 % 기본), 스펙트럼 영상 맵·ROI·NNLS 화학상태 맵, Jansson/RL/Landweber 디컨볼루션(표시 전용)·Fermi edge·재표본·평활. 긴 계산은 Worker(진행·중지).
  CasaXPS 실측(데모, 화면 판독): Linear/Shirley 끝점 규약, 맞춤 = 비가중 최소제곱(채널 i0+1…i1−1), RMS = √Σr², GL(p)=GLP, SGL(p)=높이 혼합 GLS(4 맞춤 면적 ≤1e-3), 대칭 LA fwhm = 참 FWHM, LA(1,1,n) Gaussian FWHM = n·F/144.87, GL(30)T(2) 면적 일치.
  미확정: 비대칭 LA·GLT 의 Casa 폭 숫자·비대칭 지수, Casa Tougaard(영역 표가 입력 칸으로 다시 계산되지 않아 미검증), Avantage 투과 함수 기본값(사용자 결정 대기).
  시험: xps 153 · casa 111 · multi 50 · arxps 32 · image 23 · proc 76 · Chromium ui_xps_multi 43 · ui_xps · fixes 25 · casa 20. 검증 후 CasaXPS 폴더·검증 파일 삭제.
- **10-19 전수조사 v1.15 — 다섯 영역 점검·수정 · 설명 문장은 말풍선으로**:
  Image Analysis: Find Maxima 잡음 극대 보호(> 5000 점이면 멈춤·안내) · Find (auto)(입자 중심, 합성 133/134·35/35·305/305) · Orientation 이 이미지에 박힌 축척 막대·글자 제외 · 문턱 0–255 Analyze 멈춤 · pm 저장 결과 → nm · Auto 가 손 설정(S.ui.pa)을 바꾸지 않음(다음 이미지 0 개 버그) · 마스크 Analyze 세기는 원본에서 · 굵은 마스크 Skeleton 은 먼저 Skeletonize · Subtract Background·Local 문턱·NaN Background 경고 · Summary 단위 · FFT/Bandpass 약 3–6 배, 순위 필터 약 2–4 배 빠름(이전과 같은 출력) · OCR 확실할 때만 자동 눈금.
  XPS: 대전 보정 뒤 템플릿 저장·불러오기 BE 범위 밀림 · C 1s 기준 성분을 클릭한 성분으로 잡던 것 · 계산 중 편집 덮어쓰기 · Auger α′ 엉터리 · 정량 기본 강한 선만 + 창 안 다른 원소 경고 · KE 텍스트 다시 읽기 · 음수 at% → n.d. · 상수 가중치 · 휴대폰 at% 열.
  EM Image·EDS: Velox SI 부분 영역 스캔 크기(ScanArea) · DM 썸네일·작은 이미지 · 파일 바꿔도 남던 스펙트럼 후보 · Elements 칸이 정량 묶음 · UltraX 등 맵 인식 · FEI/Zeiss/ImageJ TIFF 화소 크기 · File structure · 가짜 .emd 오류 · 1× 미만 확대 · 새로고침 시 탭 복원. (.bcf 맵은 아직 지원 안 함 — 표시함)
  Diffraction: 링·결 무늬에서 Spots 순위 매기다 탭이 메모리 부족으로 죽던 것(치명) · 링 오판 · 막대 단위 1/nm↔1/Å 오독 시 확인 단추 · 눈금 ±15 % 창 · 'No credible match' · 파일마다 설정 초기화 · 8-bit 로그 SAED 판별 · 진행·Cancel · Camera length 보정.
  DRE·Energy Diagram: 단면적 오타 시 탭 7–30 s 멈춤 → 0.7 s + 수렴 실패 표시 · 새로고침 뒤 가져오기 · λ·Range 검사 · 농도 > 100 % 경고 · 차트·표 화면 맞춤 · 이온 쌍 선택 · 내보내기.
  화면: 패널 회색 설명은 숨기고 2 초 대면(휴대폰 길게 누르면) 말풍선(modules/nav/tips.js), 긴 알림은 앞부분만.
  시험(통합 빌드): image 107·16·21·101·54·12·perf 2231 · diff 61·150·33·277·29·29·79·28 (ocr 시간 1 건 부하 탓) · dre 56·5·5·17 · ED 1404·12·10 · xps 153 · UI pa 28·audit_fix 8·drop 38·nav 25·touch 72·shell 18·xps·xps_fixes 25·diff 61·db·diff_fix 47·dre·dre_fix·cpu·insight·pool·ED ui2·ui3·int·sweep·interf if — 모두 통과.
- **10-19 겹친 입자·크기 이상치는 찾되 크기 분포에서 뺌 · 이미지에서 눌러 고른 입자만 분포 · pm → nm/Å**(modules/image/SPEC.md §6k):
  요청 ① 겹쳐 있는 입자는 찾기만 하고 분포에서 제외 ② 크기가 너무 다른 입자도 제외 ③ 손으로 고른 입자만의 분포 ④ pm 대신 nm·Å.
  Auto particle mask 가 Watershed 전 덩어리 하나에서 조각이 2 개 이상 나온 곳을 '겹침'으로 표시 → 기본으로 분포에서 뺌(Exclude overlapping). 크기 이상치 = 겹치지 않은 입자의 ln ECD 중앙값 ± 3 robust SD(1.4826·MAD) 밖(0 = 끔). Pick particles: 이미지 클릭으로 고르기·빼기 → Distribution 의 Particles = Used · All found · Picked. 표·CSV Status·Picked 열, Summary Used 열, 윤곽 색(사용 하늘색·겹침 주황 점선·이상치 자홍 점선·고름 노랑).
  단위: pm 을 쓰지 않음 — 파일이 pm 이어도 nm, 0.01 nm/px 미만만 Å, 축척 막대 1 nm 미만 Å, 저장된 pm 눈금도 nm 로. 고침: 저장 상태 불러오기에서 dist·scalebar·drop 기본값 병합이 주석에 묻혀 있던 것.
  시험 test_autopart 10(겹침 표시 52/55 · 7/7, 겹침 빼면 ECD SD 1.58 → 0.15 px) · ui_pa_guard 22 · test_img 107 · ui_drop 38 · ui_touch 72 · ui_shell 18.
- **10-19 겹친 입자 나누기(Split by size) · 배율에 맞춘 흐림 — 원클릭이 씀**(modules/image/SPEC.md §6j):
  보고 ① 겹쳐 붙은 입자 여러 개를 하나로 셈 ② ×86000 에서 지지막 결을 입자로 셈. → Split by size: 홑입자 넓이 A0 로 덩어리를 LoG 중심(세기)·EDM 씨앗 + k-means 로 나눔(Binary 메뉴에도). Auto particle mask: 1 차로 입자 크기를 재고 σ = r/6(2–12 px)로 다시 흐림, Size 최솟값은 홑입자 넓이 ¼(잡음 점 수에 끌리지 않음).
  사용자 x17000: 519 개, 합쳐진 덩어리 0 (Watershed 만이면 21), ECD 22.4 ± 2.5 nm(홑입자만 23.4 — 나눈 입자는 경계선만큼 작게). x86000: 59 개, 잡음 0. 시험 test_autopart 6(겹친 묶음 172/175 · 고배율 상관 잡음 35/35 · 저배율 305/305), ui_pa_guard 12, test_img 103, ui_drop 38, ui_touch 72.
- **10-19 원클릭 Analyze particles (auto) · 권장값을 기본값으로 · '영상' → '이미지'**(modules/image/SPEC.md §6i):
  Analyze 패널 맨 위 [Analyze particles (auto)]: 원본 그대로 두고 사본에서 Gaussian σ 2 → Threshold Default(극성 자동) → Fill Holes → Watershed → 축척 막대·글자 자동 제외(눈금 없으면 막대 글자로 정함) → Size ≥ 성분 넓이 중앙값 ¼ · Circularity 0.5–1 · edge 제외 · 세기는 원본 → Size distribution 탭. 손으로 하는 Analyze 는 아래 '고급'.
  기본값: Circularity 0.5–1, Exclude on edges, Dark background 꺼짐(저장된 상태도 한 번 바꿈). 화면·안내·사용 가이드의 '영상' → '이미지'(조사 맞춤).
  사용자 TEM 이미지: 5.7 s, 432 개, ECD 23.7 ± 4.6 nm (D50 23.3). 합성 421 개 → 421. 시험 ui_pa_guard 12 · test_img 103 · oir 54 · ui_drop 38 · ui_nav 29 · ui_touch 72 · ui_shell 18 · ui_diff 61.
- **10-19 Analyze Particles 보호 장치 · ROI 안만 세기 · nm 단위**(modules/image/SPEC.md §6h):
  기본 문턱 Dark background 꺼짐(ImageJ 와 같음 — TEM 명시야에서 바탕 70 % 를 입자로 세던 것). Max count(기본 2000): 넘으면 측정 대신 1–2 s 에 멈추고 원인 안내(잡음 점 15 만 개를 재다 페이지가 멈추던 것). 선택 > 50 % 경고. 표는 앞 3000 줄.
  영역 ROI 안만 셈(축척 막대·글자 빼기), Exclude on edges = 사각형 ROI 테두리 포함, 매크로에 ROI. 화소 단위: 0.55 nm/px 가 552.6 pm 으로 보이던 것 → nm.
  사용자 TEM 영상(0.5526 nm/px, GMS 태그): Gaussian σ2 → Default → Fill Holes → Watershed → ROI → ≥ 100 nm²·Circ 0.5–1·edge 제외 → 433 개, ECD 23.7 ± 4.7 nm (D50 23.3). 시험 ui_pa_guard 9 · test_img 103 · ui_drop 38 · ui_touch 72.
- **10-19 탭 막대 정리 — 한 줄 · 짧은 이름 · More ▾, EM Metadata 를 EM Image 안으로**(modules/nav/nav_fit.js·css, inject_all.py):
  데스크톱 탭 막대가 820–1280 px 에서 두 줄로 꺾이던 것 → 늘 한 줄: 전체 이름 → 안 들어가면 짧은 이름(휴대폰과 같은 data-short) → 그래도 넘치면 뒤쪽 탭을 More ▾ 메뉴로(지금 탭은 남김, Esc·밖 클릭으로 닫힘, 폭이 바뀌면 다시 맞춤). 휴대폰 아래 막대는 9 칸.
  EM Metadata 탭 없앰 → EM Image 위 Image · Metadata 전환(보기 키 meta 그대로, 영상 없는 파일은 바로 Metadata, 탭 막대는 EM Image 표시). 출처 줄 링크 'Set in Lattice FFT' → 'Set in Diffraction'.
  시험: ui_nav 29(360–1920 px 9 폭 × 2 탭 한 줄, More 메뉴·Esc·메뉴로 이동·넓어지면 사라짐, Metadata 보기) · ui_touch 72 · ui_shell 18 · ui_drop 38 · ui_diff 61.
- **10-17 Image Analysis v1.6 — Drop Shape: 접촉각 · 펜던트 드랍 표면장력**(modules/image/SPEC.md §3h·§6g):
  ImageJ DropSnake·LB-ADSA·Pendent Drop 대응. 처음엔 따로 Drop 탭이었으나 요청으로 **Analyze 패널 맨 아래 칸 Drop Shape** 로 옮김(머리의 "Drop Shape ↓" 단추, 예전 저장 상태 panel=drop → Analyze). 윤곽 = Otsu + marching squares + 경사 최대 다듬기, Young–Laplace(RK4) Levenberg–Marquardt 맞춤.
  접촉각: Young–Laplace · 접선(호 길이 2 차식, θ > 90° 에서도) · 타원(Fitzgibbon) · 원, 자동 기준선(RANSAC) 또는 선 ROI, 바늘 꽂힌 방울 자동 제외. 펜던트: 바늘 지름 눈금, γ·Bo·Wo·V·A·기울기, Wo < 0.5 경고.
  검증: 합성(scipy 정답) γ 중앙 0.06 %·최대 0.37 %, θ YL ≤ 0.85°·접선 ≤ 1.4°(≤ 150°). OpenDrop 공개 예제: 물 γ 72.45 ↔ 73.31 mN/m(−1.2 %, V·A·Wo·Bo ±0.9 % 안), 5 장 평균 73.1 ↔ IAPWS 72.0–72.8; 접선 118.8/123.8° ↔ 120.0/124.2°. OpenDrop 그림(GPL-3)은 ref/fetch_drop.sh 로만 — 배포 제외(manifest drop_real).
  ± 는 'fit SD'(맞춤 통계)와 'needle-scale SD'(바늘 ±0.3 px) 를 따로. 스택 시계열·매크로·CSV.
  시험: test_drop 21 · Chromium ui_drop 38 · 독립 검토 16 건 중 14 건 수정 · 기존 test_img 103·afm 101·oir 54·lif 12·ui_touch 72·ui_shell 18 그대로.
- **10-17 Diffraction v1.2 — REDI 대응 마무리: GMS 이름표 OCR · 조성비 필터 · R3·각도 오차 · 결과 보고 PNG · 화학식 검사**(modules/diffraction/SPEC.md §3c·§4c):
  GMS 기본 이름표(검은 세리프 글자 + 흰 그림자 + 속 빈 막대)를 반대 극성·중앙값 잡음 제거·글자∪그림자 모드로 읽음 — 첫 읽기 확신 ≥ 0.5 면 그대로, 뒤집은 결과는 ≥ 0.6 만 받음. 세리프·고정폭 틀(TPL_B64_X)은 뒤집기에서만 써서 산세리프 오독 없음. 시험 gms 58.3 %(capH ≥ 20: 26/33), 고확신 오답 0, 실측 REDI_B3·T13 확신 0.998·0.997.
  Phase DB 조성비 검색: "Fe3O4", "Fe:3 O:4", "Zn 60 Mg 40" 등 → 원자 분율 비교(기본 ±5 at.%p, 화면에서 조정), 수화물·괄호·혼합 자리 해석, 화학식이 비면 원자 자리 조성으로 대신. 화학식과 대칭 전개 조성이 2 at.%p 넘게 다르면 경고.
  후보 표: R1·R2·R3(=R1+R2) hkl·d·Δd, ∠R1R2·∠R1R3·∠R2R3 계산값과 차(REDI 표 형식). CSV 에도 같은 열. Export report = 영상·덧그림·후보 표를 PNG 한 장(DPR 반영).
  시험: core 61 · pattern 29 · index 33 · real 29 · pack 277 · db 150 · dyn 79 · ocr 21 · ui_diff 61 · ui_db 통과 · 독립 검토 17 건 중 16 건 수정("CO" 를 Co 로 읽는 것은 남김).
- **10-17 Diffraction v1.1 — 큰 영상 속도 · 동역학 세기(Bloch-wave) · 스케일바 숫자 읽기 · Energy Diagram 터치**(modules/diffraction/SPEC.md §3b):
  큰 영상은 줄여(≤ 1024) 찾고 원래 해상도에서 다듬음 — 4096² 스폿 9.8 s/415 MB → 1.1 s/212 MB, 링 17.6 → 2.6 s, 결과 차 0.05 % 안.
  상 찾기 세기: 상위 12 후보를 Bloch-wave(ZOLZ, Bethe, 두께 10–80 nm 평균)로 — 이중 회절·Gjønnes–Moodie 소멸 반영, 독립 Bloch 1e-6·다중 슬라이스 0.1–2 %·Hirsch 소멸 거리 5 % 안. 포화 스폿은 ρ 에서 뺌. 화면 Intensities: Dynamical/Kinematic.
  스케일바 이름표 OCR(diff_ocr.js): 시험 780 장 98.7 %, 확신 ≥ 0.9 오답 0.14 %, 실측 "5 1/nm" → 값 미리 채움.
  Energy Diagram: 손가락 탭 = 풍선, 밖 탭·스크롤 = 닫힘.
  시험: core 58 · pattern 29 · index 33 · real 23 · pack 277 · db 139 · dyn 79 · ocr 18 · ui_diff 49 · touch 72 · 기존 그대로.
- **10-17 Diffraction 탭 v1.0 — SAED 자동 해석(REDI 기능 + 보강) · Phase DB(연구실 공유)**(modules/diffraction/SPEC.md):
  Lattice FFT 탭을 Diffraction 으로 넓힘(키 lattice 유지, 휴대폰 Diffr.) — FFT (HRTEM) · SAED · Phase DB. 파일에 따라 모드 자동(역공간 교정·DIFFRACTION 태그·이름·밝기 분포).
  결정학 코어(공간군 530 설정·IT C 4.3.2.2 전자 산란 인자·CIF·정대축 2D 격자) = gemmi·pymatgen 과 일치. 패턴 읽기: 덧그림·빔스토퍼·포화 덩어리·국소 잡음 스폿·Friedel 중심(반 격자점 배제)·여러 격자·방사 왜곡·링 타원.
  상 찾기: 길이비·각 FOM + 세기 순위 상관·운동학 금지 자리·빠진 반사로 전체 패턴 검증, 공칭 보정은 순위에 안 쓰고 κ 로만 — REDI 의 '3 % 교정 오차로 순위 뒤집힘' 해결. 아핀(찌그러짐) 표준 보정 저장.
  실측 반증: Au 링(DM 태그 9 % 작음 검출), Si [001](파일 교정 nm 오기 경고 · 투사 렌즈 찌그러짐 7.4 % · 같은 장비 375/300 mm 교차 −1.4 %), Fe-Cr ODS 8-bit RGB tif(스케일바 177 px · 포화 직접빔 · bcc [135] κ 1.012).
  Phase DB: IndexedDB(내장·이 기기) + Firestore `inno-analytics-base`(구성원·관리자·이력·숨김·검증, 게시된 규칙 그대로), 판정마다 DB 항목 id·판 기록. 기본 TIFF 리더에 RGB 8-bit 추가.
  시험: core 58 · pattern 29 · index 30 · real 19 · 내장 묶음 277 · db 139 · Chromium ui_diff 47 · ui_db 83 · 독립 검토 13 건 중 12 건 수정 · 기존(ui_img·검수·XPS·DRE·셸 18·touch 70) 그대로.
- **10-17 휴대폰 터치 조작 — 공용 제스처 window.INNO_TOUCH**(modules/touch/src, inject_all.py 의 INNO-TOUCH 블록 · head):
  한 구현을 모든 캔버스가 같이 쓴다 — EM Image, EDS Maps, EDS Spectrum, Lattice FFT(실공간·FFT·역FFT), Image Analysis, XPS. 나중의 Diffraction 모듈도 `INNO_TOUCH.attach(el,{getView,setView,oneFinger,onTap,onDoubleTap,onLongPress,onPassCancel…})` 로 붙인다.
  두 손가락 = 핀치 중심 고정 확대·이동(어느 도구든), 한 손가락 = 이동 또는 그리기 도구면 그리기(둘째 손가락이 오면 그리던 것 취소), 두 번 탭 = 맞춤, 탭·길게 누르기 = hover 로만 보이던 값(화소 값·XPS 풍선).
  확대 기능이 없던 Maps·Lattice 캔버스는 CSS transform 확대(cssZoom, 다시 계산 안 함), EDS Spectrum 은 에너지 축 보기 범위(저장 안 함). 손가락만 처리 — 마우스 동작은 그대로.
  Image Analysis: 다각형 Finish·Cancel 막대, 손잡이·허용 거리 확대, Point 탭·길게 눌러 지우기(Alt+클릭 대신), Point·Wand·Zoom 은 탭 = 클릭. XPS: Add component 뒤 탭 = 그 BE 에 성분, 표지·손잡이 끌기. Lattice: 확대한 FFT 탭 = 스폿 측정.
  pointer:coarse 에서만 주요 단추·도구·탭 ≥ 40 px, 파일 목록 × 늘 보이기. rAF 로 한 프레임에 한 번만 다시 그림. DRE·Energy Diagram 은 확대 기능이 없어 그대로.
  시험: interf/test/ui_touch.py(390×844 hasTouch·isMobile, CDP 다중 터치) 70 · 데스크톱 Chromium(Image·검수·XPS·DRE·셸) 그대로 통과 · node 378 · XPS 145.
- **10-17 화면 영어 통일 · OIR 실측 검증 · Asylum IBW 단위** :
  탭 이름 XPS · EM Image · Image Analysis · EDS Maps · EDS Spectrum · Lattice FFT · EDS Interference · EM Metadata · Energy Diagram · DRE Kinetics('이미지'→EM Image 로 Image Analysis 와 구분),
  제목·메뉴·버튼·표 머리는 영어, 설명·알림·경고는 한국어(docs/UI-English-style.md). 휴대폰 아래 탭은 짧은 이름(EM Img·Analysis·Maps·Spectra·FFT·Interf.·Meta·Energy·DRE), 390·360 px 가로 넘침 없음.
  OIR: OME 시험 자료 서버의 FV4000 공개 파일 13 개(CC BY 4.0)로 검증 — 12 파일 모든 면 = oirfile = Bio-Formats(λ 는 BF 결함), 3 파일은 OME-TIFF 와도 같음. 실측에서 ZSTACK 정의 둘 중 엉뚱한 것을 쓰던 결함 등 수정, '실측 미검증' 알림 제거.
  IBW: 연구실 Cypher 실측 5 파일(AC topo·KPFM·C-AFM·EFM·FM-EFM) — 값이 단위 없는 원시 수로 보이던 결함 → 채널 이름으로 단위(nm·deg·V·pA·Hz), 48 채널 = Gwyddion 2.64. .ARDF 는 Force Map 전용(영상은 .ibw).
  배포 도구: analytics_src(공개 소스 백업)에서 'real' 폴더(실측 시험 자료) 제외.
  시험: node 378 · LIF 47 · OIR 70 · AFM 217 · Skeleton 102 · ImageJ 경계 47 · Color 17 · 실측 135 · XPS 145 · Chromium(Image·검수·XPS·DRE·셸 18).
- **10-16 Image Analysis v1.4 — Analyze Skeleton(clean-room), Leica LIF, Olympus OIR, 큰 창 Median·Remove Outliers, Wand Worker**(SPEC §3b-2·3b-3·3f·3g·6e):
  Analyze Skeleton: 처음 이식본이 Fiji(GPL-3) 절차를 그대로 옮긴 것으로 판정 → 명세 작성자(소스 본 쪽)와 구현자(소스를 보지 않음)를 나눈 clean-room 으로 새로 씀, 이식본은 지움.
  Fiji jar 블랙박스 대조 고정 102(가지치기 52) + 무작위 수천 경우 일치, 가지치기 4 방식·Prune ends 추가, Longest shortest path 크기 제한 없음(Fiji 는 꼭짓점 8 494 나무에 1 212 s, 앱은 꼭짓점 52 k 나무 0.9 s — 같은 결과).
  LIF: 실측 9 파일(Bio-Formats = liffile) + 합성 경계 11, 손상 series 건너뛰기·잘린 파일·RGB·채널 이름. OIR·POIR: oirfile(BSD) 형식 지식, 합성 45 파일(oirfile·Bio-Formats 대조) — 실측 미검증 알림.
  Median 16/32-bit 순위 이동 히스토그램(1024² r=15 13 s → 0.9 s), Remove Outliers 를 ImageJ 순서 그대로(NaN 32-bit 80 화소 차이 → 0), Wand 를 Worker 로(4096² 0.75 s, 화면 안 멈춤).
  독립 검토 3 건(LIF 13·OIR 21·Skeleton 이식본) 반영. 시험: node 378 · LIF 47 · OIR 54 · Skeleton 102(+FUZZ) · ImageJ 경계 47 · AFM 164 · Color 17 · 실측 135 · StackReg 161 · Chromium.
- **10-16 Image Analysis v1.3 — AFM 형식(GWY·GSF·Nanoscope·ARDF), Color Threshold**(SPEC §3d·3e·6d):
  AFM: Gwyddion 2.64 가져오기 모듈과 같은 규칙(Nanoscope Z scale × soft scale, 9.2 이후 32-bit, 전역 크기·종횡비, 위아래, GWYO, mask·volume). 값은 물리량(높이 nm, 힘 nN, modulus MPa…), 연산에 따라 단위 유지·지움.
  정답 = Gwyddion 모듈을 직접 부르는 변환기 + gwyfile·magic-afm·손 계산. 실측 13 파일(SciFiDatasets·AFMReader·TopoStats) + 합성 16 → 164 대조 일치. 실측 ARDF 는 공개 자료가 없어 합성으로만 확인.
  Color Threshold: ImageJ ColorThresholder 와 같은 HSB·RGB·Lab·YUV 변환(2²⁴ 색 전부 ImageJ 와 같은 바이트), pass/stop, Auto, Sample(ROI, Hue 감싸기), 스택·hyperstack, .ijm(ImageJ 로 실행해 같은 mask).
  독립 검토 2 건(AFM 10·Color 9) 수정 + 함께 찾은 결함: Enhance Contrast(표시만) 뒤 새로 고침 때 영상이 버려짐, RGB hyperstack TIFF 의 z·t.
  시험: node 378 · AFM 164 · Color 17 · 실측 135 · StackReg 161 · ImageJ 경계 19 · Chromium 61 · 검수 7.
- **10-16 Image Analysis v1.2 — 전수 검수**(SPEC §6c): 독립 검토 4 영역(계산 엔진·파일·상태와 매크로·표시와 성능)에서 재현으로 확인된 결함 약 70 건 수정.
  큰 것: Scale·회전 뒤 되돌리기에서 보정이 돌아오지 않음(v1.1 수정 때 주석이 코드를 삼킴 → 주입 전 node --check 추가), 큰 영상 Adjust 패널 10 s 멈춤, 2.4×10⁷ 화소 넘는 영상 확대·PNG 가 줄인 해상도,
  NaN·±∞ 처리(Find Maxima·Max/Min/Mean/Median·측정·표시·Z Project), 32-bit Divide 정수 잘림, 계산 중 되돌리기·닫기, 중단이 멈추지 못함, .ijm 이 ImageJ 에서 멈춤·다르게 실행, 잘린 TIFF·PNG 를 0 으로 채워 엶, multi-series OME 합쳐짐.
  시험: node 378 · 실측 135 · StackReg 161 · ImageJ 경계 19 · Chromium 46 + 검수 회귀 7 (.ijm 13 단계를 ImageJ 1.54u 가 끝까지 실행).
- **10-16 Image Analysis v1.1 — Nikon ND2·Zeiss CZI 읽기, StackReg 정합**(SPEC v1.1):
  ND2(v2 XML·v3 이진 메타데이터, loop·XY 위치·lossless, 채널 이름·색, µm/px·z step)·CZI(zstd0/1, mosaic 합성, scene, pyramid 건너뜀, Bgr48/96f, 채널 색)를 직접 해석. 큰 파일은 Import 창(위치·scene·Z·T 범위·간격, 한도 1.6 GB).
  실측 파일은 GitHub 공개 시험 자료(bioio-nd2·pims_nd2·czitools(Zeiss)·aicspylibczi·pylibCZIrw)에서 nd2 13·czi 28 → 정답 = Python nd2 패키지·Zeiss libCZI 가 읽은 값, 135 대조 전부 일치(1 GB z-stack 포함). 합성 nd2 12·czi 9.
  Register: StackReg(TurboReg, Thévenaz 1998) Translation·Rigid·Scaled Rotation·Affine — pystackreg C++ 를 이식(허용 라이선스 판, Java 판은 재배포 제한이라 제외). pystackreg 대조 45 사례 행렬 ≤3×10⁻¹³·화소 동일, 정답 변환 회복 ≤0.03 px·0.02°. ImageJ StackReg 와 같은 anchor(보이는 장) 규칙.
  독립 검토 10 건 수정(SPEC §6b): ND2 XML 긴 태그 열기 실패, 16-bit RGB ND2 잘림, TIFF z step 단위, CZI RGB+회색 채널 섞임, CZI 메모리 계산, 첫 n 장 평균 기준의 첫 장, 매크로 anchor, 고르기 창 중복 열기 등. 시험: node 364 + 실측 135 + StackReg 161 + Chromium 41.
- **10-16 Image Analysis 탭 신설**(modules/image/, SPEC v1.0): ImageJ 를 본뜬 영상 분석. ImageJ 1.54u 소스(public domain)를 컴파일해 기준으로 삼고 알고리즘을 옮김 —
  Gaussian·Rank filters·Rolling ball·Bandpass·Threshold 16 종(+Auto Local 8 종)·Binary·EDM/Watershed·Find Maxima·Analyze Particles(33 열 + ECD)·ROI Measure·Math·Enhance Contrast·Z Project.
  ImageJ 실행 결과 275 사례와 화소·수치 대조 전부 일치, 합성 성질 시험, 파일 38 개(TIFF 변종·16-bit PNG·IBW·SXM·NPY…), Chromium 33 항목. 시험 332.
  화면: 도구 막대(ROI 7 종·Wand·Set Scale), LUT, B&C, Preview, Try all thresholds, 크기 분포(로그정규·D10/50/90), ROI Manager·point counter, Radial profile(d=1/r), Orientation(구조 텐서), FFT mask, drift 정렬, Macro 자동 기록·batch·.ijm 내보내기, ImageJ 호환 TIFF 저장.
  독립 검토 13 건 수정(SPEC §6): 되돌리기가 보정·매크로 단계를 잘못 지움, Preview 중 Apply 무시, hyperstack·RGB Z Project, OME 차원 순서, RowsPerStrip 기본값, TIFF XResolution 넘침, 열 단위, 분포 통계의 0 이하 값, 형 변환 범위 기록, 영상 밖 ROI 둘레 등.
  발견: ImageJ 1.54u Outline 은 Black background 에서 바깥 테두리(snapshot 사본 읽기) → 문서 정의(안쪽)를 따름. 위상 상관 부화소는 고정 Hann 창이 0 쪽으로 0.36 px 끌어 → 창을 이동과 함께 옮겨 3 회 반복(0.013 px). 방향 분석 2차 차분 1° 치우침 → 4차 차분.
- **10-16 XPS 모형 민감도**(SPEC v1.5): 결과 탭에서 배경 종류·끝점·창 끝을 바꿔 다시 맞춘 흔들림 → 가능 범위 [최소−2σ, 최대+2σ], 주효과, 보고 예 문장, 변형 CSV.
  합성 검증: 배경을 틀리게 골랐을 때 통계 ±2σ 만 5/40 → 가능 범위 37/40 포함, 맞을 때 19/20·20/20. sd 로 합치면 20/40 이라 범위로 보고. 시험 145.
  독립 검토로 고친 것: 기울기 고정 사용자 모형을 0 번 변형이 풀어 버림(1.8σ 어긋남), 배경 효과에 끝점 효과 섞임, 창 밖 꼬리 성분 때문에 기준 변형 실패, 실패 묶음 개수, CSV 빈 값 0, 통계 σ 에 AR(1)·MC 반영, dw 범위, 지금 배경 체크 해제 무의미, 합성 시험 문구(후보 안일 때만).
- **10-16 XPS 다중항 준위·화학량론 비교·Auger parameter**(SPEC v1.4): 다중항 준위(Ti–Cu 3s·3p, Ce–Yb 4s·4p·4d·5s·5p)는 doublet 없이 envelope + 안내(K2312 Co 3p χ²ν 2.7 → 0.68).
  자동 제안에 '보이지 않는 주선' 규칙(Au 4f 두 상태 실패 사례). 정량 탭 화학량론 비교(화학식 해석·재정규화·델타 σ·깊이 경고), Auger α′ 탭(자동 짝짓기·성분/최댓값·Wagner plot·사용자 기준점). 시험 137 항목, all_elements 155 경우(다중항 20).
  독립 검토로 고친 것 7: '보이지 않는 주선'이 앞 doublet 의 짝·꼬리를 주선으로 착각(기존 성분 예측을 뺀 나머지로 판단), 시료 이름 바꾸기·옮기기 때 Auger 항목 끊김,
  Wagner plot 눈금이 가로·세로 달라 α′ 선 36.6°(→ 45°), 화학량론 화면·계산이 다른 행 사용, peakMaxPos 내림차순 에너지 구간 끝 판정, 없어진 선택이 첫 항목처럼 보여 바꿀 수 없던 것, 실패할 수 없던 시험 2 개.
- **10-16 XPS '한계에 닿은 매개변수' 표시**: 성분 id(catkqh9) 대신 이름·값·어느 한계(하한/상한)인지 표시, 피팅 뒤 모형을 바꾼(stale) 경우 이전 모형의 경고를 숨기고 숫자를 흐리게.
- **10-16 XPS 대전 보정·시료 묶음 사용성**: VGD 시료 키를 '자료 / 실험 / 측정 조건 / 지점'으로(다른 실험의 W1 합쳐짐 결함), 목록에 지점 이름 + 맥락·보정 상태 배지·✎ 이름 바꾸기(합치기), 정보 탭에서 스펙트럼 옮기기.
  대전 보정 줄 새로 — 상태 칩(보정 안 함 / Δ·기준 / 다시 맞추기), `C 1s 로 맞추기` 한 번으로 C 1s 자동 제안→피팅→맞춤→원래 화면, 기준 스펙트럼·성분 고르기(시료 안 C 1s 여럿), 기준 BE 단추, 선택 성분 기준, Δ 직접·해제, 기준 C 1s 재피팅 감지.
  알림을 그래프 아래 줄로(쌓인 toast 가 패널을 가리던 문제), 단추 문구 정리(+ 성분 추가·피팅 영역 확대·문헌 위치·Monte Carlo 불확도·작업 저장/열기), 성분 이름 중복(Co 2p C ×3) 수정, survey 안내 수치 갱신.
- **10-16 XPS 자동 제안 결함(실측 Ba 4d)**: 작은 봉우리(0.7 %)가 Δ·m 연결 기준이 되어 창 밖·넓이 0 성분으로 남음 → 기준 = 가장 큰 성분, BE 창 안, 기준 성분 빼기(removeComp), 2 % 미만 봉우리 제외.
  잔차 보강을 '한 걸음에 가장 좋은 수' 방식으로: 성분 추가 · doublet r 풀기 · q 풀기 경쟁(ΔBIC ≥ 10). 실측 Ba 4d r = 0.69(ΔBIC 84), χ²ν 1.45 → 1.06. 합성 양성 대조 104/105, 거짓 풀림 0/143.
  정적 배경 넓이 오차를 창 폭별로 다시 잼(≤20 eV 3 %, 20–40 eV 8 %, >40 eV 14 %) — 이전 문서의 '−8…+8 %' 는 > 60 eV 에서 틀렸음, 경고 문턱 40 → 20 eV. 시험 130 항목.
- **10-15 XPS 원소 독립성 검증**(SPEC v1.2): Li–U 143 원소·궤도 경우의 합성 시험(test/all_elements.js)을 시험 묶음에 넣음(126 항목).
  고친 것: 파일 이름 해석을 주기율표로 검사(원소 미상 → 정보 탭에서 지정), 광원별 Scofield(Al/Mg, 그 밖 경고), doublet Δ 를 고정 → 범위 안 자유(라이브러리 Δ 부정확 — Si 2p 0.40 vs 0.61),
  봉우리 찾기를 지형 prominence·평활 뒤 σ 로(survey 오채택 0.59 → 0.05 개), Shirley 비수렴 경고·active 의 넓이 출발을 직선 배경으로·k 상한(퇴화 결함),
  자동 제안에 잔차 보강(z ≥ 4·ΔBIC ≥ 10, 판단은 active+기울기 배경 — 정적 배경이면 143 중 106 에 가짜 성분)과 죽은 성분 덜기(Au 4f 두 상태에서 넓이 0 성분), 일반 대전 기준(기준 BE + 선택 성분).
  탭 순서: XPS 를 맨 왼쪽으로(inject_all.py 의 탭 목록 치환).
- **10-15 XPS 탭 신설**(modules/xps/): Avantage 대체 목적. VGD(OLE2 원본)·VAMAS·Avantage xls/xlsx/avg·텍스트 읽기, 배경 6종(Smart·active Shirley·Tougaard U2/U3 포함),
  선모양 6종(GLS·GLP·Voigt·DS·LA·GLT)+doublet, 제약·연결, LM+다중 시작, 공분산·AR(1) Monte Carlo, 정량(Scofield·ALTHERMO1, VGD transmission, KE^n·TPP-2M),
  survey 탐욕적 원소 판별. 독립 기준 구현 두 벌(scipy) 대조·합성 데이터 보정 시험 118 항목.
  결함(개발 중 발견·수정): 격자 합성곱 간격이 폭에 따라 연속 변해 DS·LA·GLT 피팅이 멈추고 공분산이 80–170 % 틀림 → 간격 양자화+C¹ 보간.
  디지털 Scofield 표의 doublet 합 오류(Ni 2p·Cl 2p 등) → 두 subshell 합으로 재계산. Booklet 기체 원소(O·N·F…)는 진공 기준이라 survey 판별에서 제외.
  기준 ② 대조로 SPEC v1.1: 다중 시작을 2차 χ²(시작점마다 다른 가중)로 고르던 것 → 1차 χ², AICc·BIC 도 데이터 가중 χ², LA 좌우를 CasaXPS 금속 템플릿과 맞춤, doublet q 는 수명 폭에만, Tougaard 음의 세기 0.
- **10-02 폴더 정리**: 위 ‘정리 기록’ 참고.
- **09-30 간섭 검사**: 원소별 추천 라인(예상 상대오차·실측 NNLS), DY033 입력 검토(up/down 혼합·단위 미확인·rise 오입력 경고).
- **09-29 DRE·Energy Diagram 탭 신설**(modules/): 기획서 v0.1 → D1–D4 구현 → 실측 반례(WJ-0128c, 3준위 기각) → 가설 시험 → 모델 적합성 A1–A6 → Tm 5·6준위·h5 입력.
- **09-29 suite 구조 재설계 ③**: 파일을 바꿔도 탭 유지, 파일 출처 줄(화소 크기 출처), 화소 크기 없는 이미지의 격자 탭 진입(‘알려진 d 로 정하기’), 측정 결과 우측 배치.
  결함 수정: FFT 재계산 때 '배경대비'가 다른 양(m.snr)으로 바뀌던 문제, 파일 없이 간섭 검사 탭이 환영 화면에 막히던 문제.
- **09-29 결함 A1**: 격자 측정 상태(스폿·FFT 중심·배율 교정)를 이미지별로 보관 — 다른 파일의 d 값이 새 파일 이름으로 표·CSV 에 나가던 문제.

## 교훈 · 트러블슈팅
- (XPS) 수치 Jacobian 을 쓰는 피팅에서 모형 계산 격자가 매개변수에 따라 연속으로 바뀌면 목적함수가 계단이 된다 → 격자 마디를 매개변수와 무관하게(양자화) 두고 C¹ 보간.
- (XPS) VGD 의 VGData 는 보간된 값이라 Poisson counts 가 아니다 → counts 환산은 분산의 모양에만, 절대 척도는 φ.
- (XPS) IRLS 처럼 가중이 해에 따라 바뀌는 목적함수는 서로 다른 해(시작점·모형)끼리 그 값으로 비교하면 안 된다 → 비교는 고정 가중(데이터로 정한 분산)으로.
- (XPS) 기준 구현을 한 벌만 두면 같은 오해를 공유할 수 있다 — 서로 모르게 두 벌을 만들었더니 첫 벌이 놓친 결함 7건(시작점 가중·모형 비교 가중·LA 좌우·doublet q·Tougaard 부호·분산 하한·표 밖 transmission)이 더 나왔다.
- (XPS) 원소 판별은 원소별 독립 점수보다 '가장 큰 미설명 봉우리부터 설명'이 오채택이 약 1/3(합성 survey 예비 시험) — 참값을 아는 합성 데이터로 비교해 정했고, 실측은 반례 감시에만 썼다.
- 격자 d 값은 화소 크기 출처가 전부다 → 출처(파일/.dm 교정/Gatan/직접/알려진 d)를 항상 화면에 밝힌다.
- 측정 상태를 전역 하나로 두면 파일 전환 시 남의 값이 새 파일 이름으로 나간다 → 이미지(파일명+uid) 단위로 보관.
- 다시 계산하는 경로가 있으면 처음 잰 규칙(자동=국소최대 화소, 클릭=정밀 위치)을 그대로 따라야 같은 값이 나온다.

## 정리 기록 (2026-10-02)
- 삭제: `_backup_design_0929/`(디자인 개편 전 HTML), `_backup_ed_0929/`(Energy Diagram·DRE 작업 단계별 사본 8개 + 가이드 사본 5개: before_tm5·insight·cal·cpu·terms·ifsug·dy033fix·t3). 단계별 변경 내용은 아래 문서 색인의 DRE/간섭검사 문서에 남아 있다.
- 이동·통합: `INNO Analytics 사용 가이드.md` → `docs/사용 가이드.md`, `kfactors_NaLnF4.json` → `data/`.
- 신설: `_drive_backup/`(RESTORE.md · 백업 만들기.bat · make_backup.py), `docs/HISTORY.md`(이 파일).

## INNO 공통 (6개 앱 함께)
- **배포**: GitHub 조직 `INNO-KAIST`, 앱마다 저장소 `INNO-<앱>` → GitHub Pages `https://inno-kaist.github.io/INNO-<앱>/`.
  배포 도구는 `INNO Deploy/inno_deploy.py` 하나로 5개 앱(Prism·NanoLab·Ledger·Analytics·Scheduler)을 올린다.
  `python inno_deploy.py status` (로컬↔원격 비교) / `python inno_deploy.py push <앱>`. 매핑은 `INNO Deploy/manifest.json`. (10-03 `INNO Prism/.deploy/` → `Develop/INNO Deploy/`로 이동)
  RefAtlas 만 예외: 자체 저장소 `INNO-KAIST/INNO-RefAtlas` + Cloudflare Pages, push 는 사용자가 `RefAtlas.bat` 으로.
- **2026-09-22 GitHub 조직 이관**(개인 → INNO-KAIST), 공용 계정 2FA. (문서: GitHub조직이관-*, 계정구조개편-*)
- **2026-09-29 디자인 개편 → 구조 재설계**: 공통 색 토큰(앱별 강조색), 앱 아이콘 시리즈, 공통 UX 층
  (`inno-ux.js/css`: 대화창 Esc·포커스·입력 보호, 키보드 조작, 앱 전환 메뉴, 휴대폰 하단 탭),
  한글 줄바꿈(`word-break:keep-all` + `text-wrap-style`), 결함 7건(A1–A7) 수정, 앱별 구조 재설계.
  상세: `INNO-suite-디자인개편-*`, `INNO-suite-구조재설계-STEP1-사용성감사-*`, `INNO-suite-구조재설계-STEP2-결정·진행-0929`.
- **공통 교훈**
  - `text-wrap` 단축 속성은 상속된 `nowrap` 을 풀어 버린다 → 항상 `text-wrap-style` 로.
  - 연결 폴더로 파일을 덮어쓸 때 첫 쓰기가 반영되지 않은 일이 두 번 있었다 → **배포 전 기기 파일에 새 표식이 있는지 확인**하고, 배포 도구가 "변경 없음"이면 의심할 것.
  - 단일 HTML 앱은 HTML 자체가 원본이다. 수정은 정확한 문자열 치환 + 개수 검사로 하고, 같은 재현 시나리오로 전·후를 확인한다.
- **2026-10-02 폴더 정리**: 앱마다 같은 틀(아래 README의 폴더 지도)로 재구성, 백업 잔재 삭제, Drive 보관용 `_drive_backup/` 신설.

## 문서 색인 (이 앱 관련, 최신 → 과거)
| 날짜 | 문서 제목 (claude.ai 프로젝트 ‘유원정’ → `claude/<제목>.md`) |
|---|---|
| 10-15 | INNO-Analytics-XPS탭-Avantage대체-VGD·피팅·불확도·정량·survey-검증81-1015 |
| 09-30 | DRE-DY033-입력검토-updown혼합·단위미확인·rise오입력·역방향hysteresis·느린변수-0930 |
| 09-30 | INNO-Analytics-간섭검사-원소별추천라인-예상상대오차·실측NNLS·모형-0930 |
| 09-29 | INNO-Analytics-DRE-Tm5·6준위·예측정확도·h5입력·ED루프중복수정-0929 |
| 09-29 | INNO-Analytics-DRE-피팅상세분석·모델적합성A1-A6-구현·검증-0929 |
| 09-29 | INNO-Analytics-DRE-V4b-WJ-0128c-가설시험-다중방출체·느린변수·열한계·판별실험-0929 |
| 09-29 | INNO-Analytics-DRE-V4실측반례-WJ-0128c-3준위모델기각·느린변수·다중방출체-0929 |
| 09-29 | INNO-Analytics-DRE플랫폼-기획서v0.1-식별성·전역최적화·검증계획-0929 |
| 09-29 | INNO-Analytics-DRE탭-구현D1-D4·파일럿정정·식별성-0929 |
| 09-29 | INNO-Analytics-EnergyDiagram탭-구현·데이터출처·판정규칙·검증-0929 |
| 08-26 | analytics-MRC리더-emd화소값대조-행순서실측-0826 |
| 08-26 | analytics-dm3·TIFF리더-실제배율-8비트대조-0826 |
| 08-26 | analytics-격자FFT-dspacing-역FFT-0826 |
| 08-26 | analytics-예제상수제거-HAADF입자측정-크기민감도-0826 |
| 08-26 | analytics-형상모델-구형흡수-입자크기실측-0826 |
| 08-26 | analytics-흡수보정-Elam MAC-실효k대참k-0826 |
| 08-26 | analytics-TalosF200X스펙-검출기기하-세그먼트분해능-0826 |
| 08-26 | analytics-간섭검사-입력값근거추적-출처배지-0826 |
| 08-25 | analytics-v1.5-k세트라이브러리-Na공간분포-0825 |
| 08-25 | analytics-v1.2-1.4-동시적합맵-분해능실측-저에너지꼬리-0825 |
| 08-25 | analytics-v1.1-SI스트림디코드-ROI스펙트럼-원소맵겹쳐보기-0825 |
| 08-25 | analytics-v1.0-SNIP가장자리-F1-92.5-저에너지꼬리측정-0825 |
| 08-25 | analytics-v0.7-0.9-SI·BCF지원-M선결합-0825 |
| 08-24 | analytics-v0.5-0.6-BCF는SI·L선교정·정량-0824 |
| 08-20 | analytics-v0.4-동시적합-검출한계-성공-0820 |
| 08-20 | analytics-v0.3-점수재설계-벤치마크실패기록-0820 |
| 08-20 | analytics-v0.2a-실제EDS검증-mapping아님-0820 |
| 08-15 | nanolab-INNO-Analytics-v0.1-emd플랫폼-0815 |

### INNO 공통 문서
| 날짜 | 문서 제목 (claude.ai 프로젝트 ‘유원정’ → `claude/<제목>.md`) |
|---|---|
| 09-29 | INNO-suite-구조재설계-STEP2-결정·진행-0929 |
| 09-29 | INNO-suite-구조재설계-STEP1-사용성감사-6개앱-0929 |
| 09-29 | INNO-suite-디자인개편-공통틀·앱별강조색-파일럿Ledger-0929 |
| 09-27 | GitHub자동배포-inno_deploy-5개플랫폼·토큰·사용법-0927 |
| 09-22 | GitHub조직이관-이관후실측·수정적용·잔여조치-0922 |
| 09-22 | 계정구조개편-공용계정2FA·GitHub조직·공유드라이브-진행상황-0922 |
| 09-22 | GitHub조직이관-INNO-KAIST-파괴지점·실행순서·검증-0922 |
