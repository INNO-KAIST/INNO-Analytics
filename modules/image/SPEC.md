# Image Analysis 탭 — SPEC v1.3 (2026-10-16)

ImageJ 를 본뜬 영상 분석 탭. 계산 결과가 ImageJ 와 같도록 ImageJ 1.54u 의 알고리즘을 옮겼다.
ImageJ 1.x 는 미국 NIH 직원이 직무로 작성해 public domain 이다(`LICENSE.txt`). 옮긴 원본 파일은 각 함수 머리말에 적었다.

## 1. 구조
| 파일 | 내용 |
|---|---|
| `src/core/*.js` | 계산 엔진 `IMG_CORE`. 이름 순서로 이어 붙여 Web Worker 원본(`<script id="img-core-src">`)이 된다 |
| `src/img_io.js` | 파일 읽기·TIFF 쓰기 `IMG_IO` |
| `src/img_micro.js` | Nikon ND2 · Zeiss CZI 읽기 (`IMG_IO.readND2/readCZI`) |
| `src/img_afm.js` | AFM·SPM 읽기: Gwyddion GWY(GWYP·GWYO)·GSF, Bruker Nanoscope, Asylum ARDF (`IMG_IO.readGWY/readGSF/readNanoscope/readARDF`) |
| `src/core/80_color.js` | Color Threshold (HSB·RGB·Lab·YUV, pass/stop, Auto, Sample) |
| `ref/gen_afm.py`, `ref/gwyconv/gwyconv.c` | AFM 합성 파일 + 정답(Gwyddion 2.64 가져오기 모듈·gwyfile·magic-afm·손 계산). gwyconv = Gwyddion 모듈을 직접 불러 GWY 로 쓰는 변환기 |
| `test/test_afm.js` | AFM 시험 164 (합성 16 파일·실측 13 파일) |
| `test/test_color.js` | Color Threshold 시험 17 (ImageJ 실행 대조: 2²⁴ 색 전부, RGB hyperstack TIFF) |
| `src/zstd_vendor.js` | fzstd 0.1.1 (MIT, 원본 그대로) — CZI zstd 해독 |
| `src/ui/*.js` | 화면. `inject_all.py` 가 하나의 IIFE 로 묶는다 |
| `src/img.css` | `.imx` 범위 스타일 |
| `ref/gen_ref.py` | ImageJ 실행으로 기준 결과(`ref/data`) 생성 |
| `ref/gen_io.py` | 파일 읽기 정답 파일(`test/io`) 생성 |
| `ref/gen_micro.py` | 합성 nd2·czi + 정답(nd2 패키지·czifile·libCZI 가 읽은 값) |
| `ref/fetch_real.sh`, `ref/gen_real.py` | 실측 nd2 13·czi 28 내려받기(GitHub 공개 시험 자료) + 정답(nd2 패키지·libCZI) — 파일은 저장소에 넣지 않음 |
| `test/test_real.js` | 실측 nd2·czi 대조 |
| `test/test_img.js` | node 시험 332 항목 |
| `test/ui_img.py` | Chromium 시험 61 항목 |
| `test/ui_audit.py` | 전수 검수 회귀(계산 중 되돌리기·닫기, 중단, .ijm 을 ImageJ 가 끝까지 실행 등) 7 항목 |
| `test/test_ij_edge.js` | NaN·±∞·일정한 영상 경계 사례를 ImageJ 를 직접 실행해 대조 19 항목 |
| `ref/gen_audit.py` | 파일 읽기 회귀 시험 파일(`test/io/audit`) |
| `src/core/75_stackreg.js` | StackReg/TurboReg 엔진 (pystackreg C++ 이식) |
| `test/stackreg_parity.py`, `test/test_stackreg.js` | pystackreg 대조 45 사례 |

화소 표현: plane = `{w, h, bit: 8|16|32, data: Float32Array}`. 8/16-bit 도 정수를 Float32 에 담고, ImageJ 와 같은 규칙 `(int)(v+0.5)` 후 범위 자르기로 다시 정수화한다.
영상 = planes 배열, 순서 c + nc·(z + nz·t) (ImageJ hyperstack 과 같음). RGB 는 8-bit plane 3 장 + `rgb` 표시.

## 2. 연산과 원본
| 범주 | 연산 | ImageJ 원본 |
|---|---|---|
| Filters | Gaussian (kernel accuracy 0.002/0.0002, σ>8.5 축소 경로), Mean/Median/Min/Max/Variance/Remove Outliers/Despeckle/Top Hat (원형 kernel), Smooth/Sharpen/Find Edges, Convolve, Unsharp Mask | GaussianBlur, RankFilters, Byte/Short/FloatProcessor.filter3x3, Convolver, UnsharpMask |
| Background | Subtract Background (rolling ball, shrink 2/4/8, arc trim 24/32/40 %, 3×3 pre-smoothing) | BackgroundSubtracter |
| FFT | Bandpass (mirror tiling, Gaussian band, stripe suppression, autoscale/saturate), power spectrum, ROI mask inverse FFT | FFTFilter |
| Contrast | Enhance Contrast (saturated %, normalize, √/classic equalize), Apply LUT | ContrastEnhancer |
| Math | Add…Reciprocal, NaN background, Invert | Byte/Short/FloatProcessor.process |
| Threshold | 16 global methods, 16/32-bit 은 8-bit 비례 히스토그램 → 원래 단위로 되돌림, Convert to Mask | AutoThresholder, ImageProcessor.setAutoThreshold |
| Local threshold | Bernsen·Contrast·Mean·Median·MidGrey·Niblack·Phansalkar·Sauvola (원형 창) | Fiji Auto Local Threshold 의 문서 정의 |
| Binary | Erode/Dilate/Open/Close (count, iterations), Fill Holes, Outline, Skeletonize, Distance Map, Watershed, Ultimate Points | ByteProcessor.filter, BinaryProcessor, Binary, EDM, MaximumFinder |
| Maxima | Find Maxima (prominence, strict, exclude edges, light) + Single/Within tolerance/Segmented | MaximumFinder |
| Measure | Area·Mean·StdDev·Mode·Min·Max·Centroid·Center of mass·Perimeter·Bounding·Ellipse·Circ.·AR·Round·Solidity·Feret 5종·IntDen·Median·Skew·Kurt + ECD | ImageStatistics, EllipseFitter, Roi.getFeretValues, PolygonRoi.getTracedPerimeter/getSmoothedPerimeter, OvalRoi(Ramanujan), PolygonFiller |
| Particles | 8-연결, exclude edges, include holes(안쪽 입자를 삼킴 — ImageJ 와 같음), size·circularity 필터(화소 수 기준), redirect | ParticleAnalyzer, Wand |
| Stack | Z Project (1.54u: avg·max·min 은 원래 형식), Z-axis profile | ZProjector |
| Registration | StackReg: Translation · Rigid Body · Scaled Rotation · Affine (B-spline 피라미드, 세기 차 제곱합 최소화), 기준 = 이전 장(anchor = 보이는 장에서 앞뒤로 누적)·첫 장(첫 n 장 평균)·평균; Translation phase correlation | TurboReg/StackReg (Thévenaz 1998) — pystackreg 0.2.8 의 C++ 이식 |
| 그 밖 | Radial profile, Orientation (구조 텐서, OrientationJ 정의), size distribution (log-normal MLE) | — |

## 3. ImageJ 와 일부러 다르게 한 것
1. **Outline**: ImageJ 1.54u 는 Black background 에서 snapshot 사본을 읽어 물체 *바깥* 테두리를 만든다(`Binary.run` 의 `setSnapshotCopyMode(true)` → `BinaryProcessor.process` 가 반전 전 화소를 읽음). 문서 정의인 안쪽 1 px 테두리를 따르고, 흰 배경 경로의 ImageJ 결과와 같은지 시험한다.
2. **Phase correlation 부화소**: 이동하는 쪽 영상의 Hann 창을 추정 이동만큼 옮겨 다시 계산(3 회). 고정 창은 0 쪽으로 최대 0.36 px 끌린다(합성 시험). 결과 최대 오차 0.013 px.
3. **Orientation 기울기**: 4차 중앙 차분. 2차 중앙 차분은 주기 9 px 결에서 1° 가량 기운다.
4. ImageJ 의 Java `int` 넘침(Li·Huang 의 합, 2³¹ 초과)은 재현하지 않는다(double 로 셈).
5. **Mean/Variance + NaN**: ImageJ 는 이동 합이 NaN 을 지난 뒤 줄의 다음 화소 하나를 더 NaN 으로 남긴다(전체 재계산 전). 여기서는 커널 안에 NaN 이 있는 화소만 NaN — 유한값은 ImageJ 와 같음.
6. 16-bit 의 일부 결과(rolling ball shrink 경로·bandpass)는 ImageJ 의 float32 중간 계산과 ±1 count 차이가 0.5 % 미만 화소에서 난다(전체 눈금의 1.5×10⁻⁵).

## 3b. 현미경 형식 (ND2 · CZI)
- **ND2** (NIS-Elements): chunk map → `ImageAttributesLV!`·`ImageMetadataSeqLV|0!`·`ImageMetadataLV!`(v3, CLx Lite Variant 이진) 또는 같은 이름의 XML(v2.x). 해석 규칙은 Python `nd2` 패키지 0.12 의 modern reader 와 같다:
  실험 loop 은 바깥 loop 먼저(T·NETime·XY·Z, Spect loop 은 건너뜀, 같은 층의 같은 종류 loop 은 큰 쪽), 프레임 = `ImageDataSeq|i!` 자료의 8 바이트(시각) 뒤, 행 간격 = `uiWidthBytes`,
  component 가 채널당 3 이면 RGB, 아니면 모든 component 가 채널. lossless = zlib. 없는 프레임은 0(알림). XY 위치 = 영상 따로. 시간형 loop 이 둘 이상이면 하나의 t 로 잇는다(알림).
  v1.0(JPEG 2000)·lossy 는 분명한 오류로 알린다.
- **CZI** (ZEN): ZISRAWFILE → directory(DV 항목) → subblock. 층 0 만(저장 크기 = 논리 크기, pyramid type 0). 화소 형식 Gray8/16/32f, Gray32(int)·Gray64(→ float, 알림), Bgr24(→ RGB), Bgr48·Bgr96f(→ R·G·B 3 채널), Bgra32.
  압축: 없음·zstd0·zstd1(첫 바이트 = 머리 크기, 항목 1 = HiLoByteUnpack: 낮은 바이트 전부 → 높은 바이트 전부)·JPG(브라우저). JPEG XR·LZW 는 오류로 알림.
  mosaic: scene 의 층 0 subblock 경계 상자 위에 M 순서로 그려 뒤 tile 이 위(libCZI 와 같은 결과 — 실측 겹침 파일 포함 일치). scene·그 밖의 차원(H·R·I·V·B 중 값이 둘 이상) = 영상 따로.
  Scaling·채널 이름·색은 `ImageDocument/Metadata/Scaling`, `…/Information/Image/Dimensions/Channels`, `…/DisplaySetting/Channels` 자리만 읽는다(같은 이름의 다른 요소 무시 — LLS7 파일은 다른 곳에 다른 z 값이 있음).
- 채널 색 → 가까운 LUT(Red·Green·Blue·Cyan·Magenta·Yellow). 광학 영상 보정 단위는 µm(10 nm 미만이면 nm).
- 큰 파일: 위치/scene 이 6 개를 넘거나 32-bit 메모리로 600 MB 를 넘으면 고르기 창(위치·scene, Z·T 범위와 간격). 한도 1.6 GB. 간격을 두고 고른 Z 는 z step 도 그만큼 곱한다.

## 3c. Registration (StackReg)
- 엔진 `75_stackreg.js` = pystackreg 0.2.8 C++(TurboReg ANSI-C 이식, 허용 라이선스 — 머리말에 원문·논문 인용) 를 연산 순서까지 옮김. Java TurboReg/StackReg(재배포 제한 라이선스)는 옮기지 않았다.
- pystackreg 대조 45 사례(5 변환 × 크기 4, 음수·8-bit, 512², stack previous/first/mean): 행렬 최대 차 3×10⁻¹³(rigid — Node 와 glibc 의 atan2·sin 1 ulp 차이까지 추적), 나머지 0, 변환 영상 화소 동일.
- 정답 변환 회복(합성): translation 0.005 px, rigid 중심 0.007 px·각도 0.02°, scaled 0.004 px·0.009°, affine 모서리 0.03 px.
- pystackreg 의 특징을 그대로 따른 것: 맞춤 전 마지막 행·열을 뺌(`ref[:-1,:-1]`), Levenberg–Marquardt 단계에서 Hessian 대각만 갱신(원래 Java TurboReg 와 다름 — 위 정확도는 이 방식의 값).
- pystackreg 와 다르게 한 것 하나: 기준 '첫 n 장 평균'(n>1)이면 첫 장도 그 평균에 맞춘다(pystackreg 는 첫 장을 그대로 둬 첫 장만 어긋남). 대조 시험도 이 규칙으로.
- 매크로: '이전 장' 기준의 기준 장(anchor)을 기록해 재생·batch 에서 같은 결과.
- 앱 규칙: 기준 채널로 변환을 구해 모든 채널에 적용. z·t 가 모두 있으면 t 마다 z 최대 투영으로 구해 모든 z 에 같은 변환. 결과는 원래 bit 로 반올림(8/16-bit). 영상 밖은 0.
  보고값: drift = 영상 중심이 기준 대비 옮겨진 양(px), 회전 = 화면에서 반시계 +(Orientation 과 같은 규약), scale = √|det|. ImageJ 매크로로는 `run("StackReg ", "transformation=[Rigid Body]")`.
- 속도(512², node): rigid 0.17 s, affine 0.2 s / 장.

## 3d. AFM·SPM 형식 (v1.3)
- 해석 규칙의 기준은 Gwyddion 2.64 의 가져오기 모듈(modules/file/nanoscope.c·gwyfile.c·gsffile.c·ardf.c). 코드는 옮기지 않고 동작을 맞췄으며, 시험 정답은 그 모듈이 같은 파일을 읽은 값이다.
  `gwyddion --convert-to-gwy` 는 GTK 를 초기화하지 않아 ARDF 모듈이 실패하므로, 모듈을 직접 부르는 작은 변환기(`ref/gwyconv`)로 GWY 를 쓰고 gwyfile(독립 파이썬 파서)로 읽는다.
- **값 단위**: 영상 값은 물리량으로 바꾼다. 길이(m)는 늘 nm, N·Pa·A·F·S 는 중앙값 크기에 맞는 접두어(nN·MPa 등), 그 밖(V·deg·Hz·단위 없음)은 그대로. 배율은 정보 칸에 적고, 마우스 아래 값에 단위를 붙인다.
  같은 물리량을 남기는 연산(필터·회전·자르기·정합·평균/최대 투영 등, `KEEP_UNIT`)만 단위를 이어받고, 형 변환의 눈금 바꾸기·문턱·FFT·곱셈 등은 지운다. 되돌리면 단위도 돌아온다.
- **Nanoscope**: 머리 길이 = 앞 8 필드 안 `\Data length`. 값 = 원시 정수 / 256^(Bytes/pixel) × `@2:Z scale`(또는 `@4:`)의 hard value × soft scale(`@Sens. …`, 단위 곱 — 예 V × nm/V = nm).
  soft scale 에 단위가 없으면 V(Phase 는 deg). 9.2(0x09200000) 이후는 자료가 32-bit(Bytes/pixel 은 눈금만). 영상 절 크기가 틀리면 Ciao scan list 의 전역 크기, 종횡비는 Gwyddion 의 두 규칙.
  첫 줄이 아래(위아래 뒤집기). 채널 이름은 따옴표 안 화면 이름(Gwyddion 은 soft 이름). 힘 곡선(Force)·단면 묶음·텍스트형은 분명한 오류, Force Volume 은 영상 채널만(알림). `\*Force file list` 로 시작해도 Operating mode 가 정함.
- **GWY**: 13 성분 형식, GWYO(1.x) 의 옛 맨 위 컨테이너(GType 번호), DataField·mask(→ 8-bit 0/255 영상)·Brick(→ z-stack, z 가 길이면 z step), meta(GWYO 는 `/meta/…`). 단위 문자열의 접두어는 Gwyddion 처럼 버림(값은 늘 기본 SI).
  그래프·XYZ·curve map·spectra 는 알리고 읽지 않음. **GSF**: 머리 + NUL 1–4(4 바이트 경계) + float32, 크기가 정확히 맞아야 함, NaN 은 그대로(알림).
- **ARDF**: 블록 머리 CRC32 를 모두 확인(손상은 오류), FTOC → IMAG(TTOC·IDEF·IBOX·IDAT·GAMI). IDEF 의 간격·단위·이름을 쓴다(Gwyddion 은 FastScanSize 로 정사각 가정 — 직사각이면 틀림).
  영상은 첫 줄이 아래(IBW 와 같음; magic-afm 도 origin=lower — Gwyddion ARDF 모듈만 뒤집지 않음). 0 으로 채운 IDAT 줄(중단된 주사)은 NaN + 알림. 힘 곡선(VOLM)은 크기·채널만 알림.
  **실측 ARDF 는 구하지 못했다**(공개 저장소에 없음) — magic-afm·Gwyddion 이 함께 받아들이는 합성 파일 4 개로만 확인.
- 실측: SciFiDatasets(pycroscopy)·AFMReader·TopoStats 시험 자료 13 파일(Nanoscope 5.31·8.15·9.2·10, GWYP 3·GWYO 2, GSF, 힘 곡선 1) — 모두 Gwyddion 결과와 값·크기·단위 일치(상대 1.5×10⁻⁷, float32 반올림).

## 3e. Color Threshold (v1.3)
- ImageJ `ColorThresholder` 와 같은 색 변환: HSB(java.awt.Color.RGBtoHSB 의 float 연산 순서 그대로), Lab(ColorSpaceConverter D65, L×2.55·a/b×1.0625+128), YUV(Landini, byte 로 감음).
  ImageJ 를 실행해 **2²⁴ 색 전부** 같은 바이트임을 확인. 문턱 = 채널마다 pass(안)·stop(밖), 세 조건 AND → 8-bit mask(255). Auto = HSB 는 Brightness, RGB 는 세 채널, Lab·YUV 는 첫 채널(autoSetThreshold).
  Sample = ROI 안 범위, HSB Hue 는 쓰인 범위와 가장 긴 빈 구간 중 좁은 쪽(빨강처럼 0/255 를 감싸면 stop).
- 결과는 새 mask 영상(기본) 또는 제자리(스택은 '모든 장' 필요). .ijm 은 ImageJ 가 만드는 매크로와 같은 단계(변수 이름 ct_ — batch 고리의 i 를 덮지 않음); ImageJ 매크로는 한 장만.

## 4. 화소 크기 (calibration) 규칙
- 출처가 확실한 값만 자동 적용: Nikon ND2 `dCalibration`(bCalibrated 일 때), Zeiss CZI `Scaling/Distance`, ImageJ ImageDescription(unit)+XResolution, OME-XML PhysicalSize, FEI/Thermo 34680·34682 `[Scan] PixelWidth`, Zeiss 34118 `Image Pixel Size`, TESCAN 50431, Hitachi 옆 `.txt` `PixelSize`(nm), JEOL 옆 `.txt` 축척 막대, DM·MRC·EMD·BCF(앱 기존 파서), Igor wave scaling, Nanonis SCAN_RANGE, Gwyddion xreal/xres, Nanoscope Scan Size(전역 크기·종횡비 규칙), ARDF IDEF 간격.
- TIFF/PNG 해상도 태그만 있으면 '후보'로 보이고 적용하지 않는다(72 dpi 같은 기본값이 흔함).
- FEI 데이터 막대(`[Image] ResolutionY` < 영상 높이)는 경고와 한 번에 잘라내기.
- 크기 필터의 단위 변환: size(unit²)/(pw·ph) = 화소 수. Circularity 필터는 화소 단위 둘레·넓이(ImageJ 와 같음 — 보정과 무관).

## 5. 검증
- `node test/test_img.js` (378): 아래에 더해 파일 읽기 회귀 13(`test/io/audit`). ① ImageJ 1.54u 실행 결과 275 사례(8·16·32-bit 필터·문턱 16×2×3 histogram 모양 8 종·이진·maxima 위치·입자 33 열·ROI 8 종(영상 밖으로 나간 ROI 포함)·계산기·투영) ② 정답을 아는 합성 영상(줄무늬 방향 5, 부화소 이동 4, 고리 반지름, EDM 대 완전 탐색 거리, 원 ECD·원형도, 로그정규 회복, local Mean) ③ 파일 38(tifffile·pypng·numpy·igor2 로 만든 정답) + TIFF 쓰기→읽기(아주 작은 화소 크기, RowsPerStrip 기본값), 저장한 TIFF 를 ImageJ 가 보정·hyperstack 째 읽는 것 확인.
- 현미경 형식: 합성 nd2 12·czi 9 (`node test/test_img.js` 안, 정답 = nd2 패키지·czifile·libCZI), 실측 nd2 13·czi 28 → 135 대조 (`node test/test_real.js`, 정답 = nd2 패키지·libCZI, Bgr96Float 만 czifile).
- `python3 test/ui_img.py` (41): nd2·czi 열기(화소 = 정답, 고르기 창, 창이 떠 있을 때 또 열기), Register rigid(정답 회전·drift, 화소 = pystackreg), Worker 결과 = ImageJ, 되돌리기, Preview, Threshold 패널·Convert to Mask = ImageJ, Analyze Particles = ImageJ, 마우스 ROI·Ctrl+M = numpy, Find Maxima 개수, Macro 기록→다른 영상, .ijm, 새로 고침 복원, 모바일, 콘솔 오류 0.
- AFM: `node test/test_afm.js` (164) — 합성 16 파일(정답 = Gwyddion·gwyfile·magic-afm·손 계산) + 실측 13 파일(`ref/real/afm`, 정답 = Gwyddion) + 손상 파일 오류·단위 해석. 정답 다시 만들기: `python3 -I ref/gen_afm.py test/io/afm ref/real/afm` (gwyddion·libgwyddion20-dev·gwyfile·magic-afm 필요).
- Color Threshold: `node test/test_color.js` (17, ImageJ 실행) — 2²⁴ 색 변환, 문턱 7 경우, Auto·Sample, RGB hyperstack TIFF. Chromium: 패널·Sample·Auto, 스택·hyperstack .ijm 을 ImageJ 로 실행해 대조.
- ImageJ 기준 다시 만들기: ImageJ 소스(github.com/imagej/ImageJ) 를 javac 로 컴파일 → `python3 -I ref/gen_ref.py ij.jar ref/data`.

## 6. 독립 검토로 고친 것 (v1.0 안)
되돌리기가 보정을 지우고 엉뚱한 매크로 단계를 빼던 것(단계별 표지로) · Preview 계산 중 Apply 가 무시되던 것(전역 표시 → 매개변수) · 적용 뒤 Preview 가 한 번 더 겹쳐 보이던 것 ·
RGB Z Project 의 Sum/SD 가 깨진 영상 · hyperstack Z Project 가 시간까지 합치던 것 · OME DimensionOrder 무시(채널·z 뒤섞임) · RowsPerStrip 기본값(2³²−1) 에서 멈춤 ·
TIFF 쓰기 XResolution 넘침 · 열 단위가 지금 고른 영상을 따르던 것(측정 때 단위를 행에 저장) · 분포 통계가 0 이하 값을 조용히 빼던 것(선형 = 전부, 로그정규 = 양수, 알림) ·
열 때 표시 범위(ImageJ 와 같은 min–max) · 형 변환 범위를 매크로에 기록 · Duplicate 기록 · .ijm 에 문턱 명령 · FFT radial profile 중심 반 화소 · 영상 밖으로 나간 ROI 의 둘레·Feret(ImageJ 와 같음) ·
batch 의 redirect 를 원본 상대 지정(__parent/__root)으로 · RGB 입자 분석 막음 · Space 고정 · Flip·Register 뒤 남던 겹쳐 그림 · 보정 단위 Gaussian 의 세로 σ.

## 6b. 독립 검토로 고친 것 (v1.1)
ND2 v2 XML 의 4096 자 넘는 태그(긴 Description)에서 열기 실패 → 따옴표를 지키며 태그 끝까지 읽기, `&#..;` 해독 · 16-bit RGB ND2 를 8-bit RGB 로 잘라 버림 → R·G·B 3 채널 16-bit ·
TIFF spacing 단위(nm 보정에서 z step 1000 배) · CZI 의 RGB 채널 + 회색 채널 섞임 → 따로 영상 · CZI 고르기 창이 첫 scene 크기로 메모리 계산, 전체 한도 없음 → scene 별 크기·고른 것 합으로 한도 ·
'첫 n 장 평균' 기준에서 첫 장 미정렬 · Register 매크로 재생 시 기준 장 바뀜 · 위치 64 개 넘을 때 기본 간격 계산 · 고르기 창이 떠 있을 때 또 열면 앞 파일이 사라짐 → 차례로 ·
채널 수가 바뀐 뒤 남는 채널 이름 · 여러 scene 영상마다 메타데이터 원문 복사(저장 공간) → 첫 영상에만 · nd2 sub-loop 규칙 · CZI JPG 색 공간 변환 끔.

## 6c. 전수 검수 (v1.2) — 독립 검토 4 영역에서 확인된 결함 수정
- **계산 엔진(ImageJ 직접 실행으로 대조)**: NaN 이 있는 32-bit 에서 Find Maxima 봉우리 누락(정렬 키 (int)NaN), Max/Min·Open/Close·Top Hat 이 NaN 근처에서 틀림, Mean/Variance 의 NaN 이 줄 끝까지 번짐,
  Median 의 NaN 처리(ImageJ getNaNAwareMedian: 가운데 화소 추정·아래쪽 중앙값), 32-bit Image Calculator Divide 가 정수로 잘림, 측정·히스토그램·표시 범위에 ±∞ 포함(ImageJ 는 뺌),
  min=max 에서 8-bit 로 NaN, Z Project 의 NaN(평균·최대·최소는 뺌)·16-bit median 은 32-bit, Enhance Contrast 가 일정한 영상을 0 으로, 임의 각도 회전 가장자리 1 px, 정수·축 평행 다각형 둘레(traced),
  얇은 영상 Subtract Background → 분명한 오류, 8/16-bit ÷0 은 그대로·Reciprocal 1/0=NaN, Apply LUT 반올림, float 최빈값 칸 폭 1, 빈 histogram 문턱 범위, 일정한 16/32-bit 문턱 없음, 다각형 float32 좌표,
  NaN·∞ 가 있는 영상의 정렬 → 분명한 오류, StackReg 최소 16 px, 일정한 영상 Orientation = 방향 없음. 속도: Include holes 입자마다 전체 크기 배열(1024² 3.2 s → 상자 크기), redirect 없을 때 반복 복사.
- **파일**: TIFF·PNG·NPY 메모리 한도와 잘림·손상 → 분명한 오류(0 으로 채우지 않음), LZW 손상 검출, ImageJ 한 채널 표시 범위, multi-series OME → series 마다, OME PhysicalSizeZ, (H,W,3) uint16/float NPY,
  텍스트(';'·',' 소수점, NaN·inf, 줄 끝 구분자, 머리 줄), ImageJ signed 16-bit, 4-bit·1-bit 팔레트, ImageJ 이어진 RGB·부호 스택, int64·uint64, IBW 이름 위치(28)·음의 간격·4D·잘림,
  쓰기 해상도 유리수 정밀도·8/16-bit 반올림·자르기, StripByteCounts 없음, CMYK·YCbCr 알림, UTF-8 µ(태그·옆 파일), 경고 두 번, 파일이 밝힌 단위(micron → µm) 유지, CZI XML 경로 되돌림.
- **상태·매크로**: Scale 등 되돌리기에서 보정이 돌아오지 않던 것(이전 수정에서 주석이 코드를 삼킴 — 주입 전 node --check 추가), 계산 중 되돌리기·닫기·모두 지우기 막기, 계산 중 닫힌 영상의 화소가 되살아남·고아 화소 삭제,
  중단이 기다리던 일·batch 를 멈추지 못함(중단 세대 번호), .ijm 이 ImageJ 에서 멈추거나 다르게 실행(Flip·Rotate·Scale·Duplicate·Bandpass·Enhance·Convert to Mask·Conversions·redirect·문자열·EDM·Measure ROI·Apply LUT·대화 상자 없는 단일 장 처리·영상 파일만 열기),
  재생: Orientation 색 지도·Measure ROI·Apply LUT 범위·이진이 아닌 영상의 Analyze·없는 redirect, batch 파일 실패 시 정리·옆 .txt, Convert to Mask 가 스택을 한 장으로, 저장 실패 후 옛 화소(판 번호 확인),
  복사 TSV·plot CSV 머리, 되돌리기의 composite·databar, Summary 시간축 Slice, n=1 SD, FFT 영상 역공간 보정(1/unit), 되돌리기 메모리 계산(공유 배열), 매크로 덮어쓰기 확인, IndexedDB 연결 재사용.
- **표시·상호작용**: 큰 영상에서 Adjust 패널이 다시 그릴 때마다 10 s 멈춤(통계 캐시·형식 배열 정렬), 2.4×10⁷ 화소 넘는 영상 확대 시 원래 해상도(보이는 영역만)·PNG 내보내기 원래 해상도,
  ±∞ 영상이 검게, 끝난 일꾼이 큰 배열을 붙잡음, 닫은 영상의 표시 버퍼, 막대 끌기(한 프레임에 한 번)·16-bit 65536 칸 표, 터치로 다각형, min=max·max<min 표시(ImageJ 와 같은 이진),
  16-bit 표시 매핑 = ImageJ ShortProcessor, 재생 중 fps, 썸네일 갱신, 눈금 막대(보이는 영상 안·단위 자동·글자 크기 비례), 오른쪽 단추로 ROI 그리지 않음·ROI 좌표 내림(ImageJ), Esc 로 창 닫기, HiLo composite.

## 6d. 독립 검토로 고친 것 (v1.3)
- AFM: `\*Force file list` 로 시작하는 Force Volume 을 힘 곡선으로 거부, 0 으로 채운 ARDF IDAT 줄에서 파일 전체 실패, `Data offset` 이 수가 아니면 머리를 영상으로 읽음,
  Nanoscope E(`Samps/line: 512 512`), Scan Size 두 수 사이 빈칸 여럿, ARDF TEXT 가 24 바이트보다 짧을 때, `º`(0xBA) 단위, 파일 고르기 창의 .100–.999, GWYO meta. — 각 경우 합성 파일과 Gwyddion 대조 시험 추가.
- Color Threshold: 변환 중 Auto·Sample 이 이전 색 공간·다른 영상의 캐시로 계산(→ 준비될 때까지 끔·확인), 스택 .ijm 이 ImageJ 에서 멈춤·한 장만(→ 장마다 고리·Images to Stack·hyperstack 되살림, ImageJ 로 실행해 화소·차원 대조),
  큰 영상에서 끌 때 전체 세기(→ 4 M 화소 넘으면 듬성 '≈', 놓으면 정확히), 변환 일이 쌓임(→ 한 번에 하나, 마지막 요청만 이어서, Uint8 로 보냄, Lab 선형화 표), 지금 장 다시 변환(→ 캐시 사용),
  실패 뒤 '변환 중' 고정, 재생 중 겹쳐 그림 사라짐, 닫은 영상 캐시 남음, min > max 허용(→ ImageJ 처럼 min ≤ max).
- 함께 찾은 결함: Enhance Contrast(표시만) 가 영상 ver 를 올려 새로 고침 때 그 영상을 버림. RGB hyperstack TIFF 의 z·t 를 읽고 쓰지 못함(z 로 합침).

## 7. 한계 (다음 단계)
- lif·oir·vsi, ND2 v1(JPEG 2000), CZI JPEG XR, JPEG 압축 TIFF 는 아직 읽지 않는다. AFM 힘 곡선(Force·VOLM·curve map)은 읽지 않는다.
- Analyze Skeleton, 3D 연산, bilinear·비강체(elastic) 정합은 없다.
- 2.4×10⁷ 화소를 넘는 영상은 전체를 볼 때 줄여 보이고(확대하면 원래 해상도), 문턱 겹쳐 그림은 줄인 해상도.
- Wand 는 주 스레드(2048² 균일 영역 0.6 s). 16/32-bit Median·Remove Outliers 는 큰 반지름에서 느림(1024², r=15: 13 s). 터치 두 손가락 확대 없음.
