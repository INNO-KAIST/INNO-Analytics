# INNO Analytics — 작업 이력

> STEM 이미지·EDS 스펙트럼·SI(원소 맵)·격자 FFT·EDS 간섭 검사·Energy Diagram·DRE(속도식 모델)·XPS·Image Analysis(ImageJ 방식) 분석을 브라우저 안에서 하는 단일 HTML 앱. 서버·로그인 없음.
>
> 이 파일 하나에 이 앱의 작업 이력·교훈·문서 색인을 모았다. 자세한 근거는 색인의 문서(claude.ai 프로젝트 ‘유원정’)에 있다.

## 최근 작업
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
