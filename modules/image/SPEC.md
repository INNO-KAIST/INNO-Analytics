# Image Analysis 탭 — SPEC v1.6 (2026-10-17)

ImageJ 를 본뜬 이미지 분석 탭. 계산 결과가 ImageJ 와 같도록 ImageJ 1.54u 의 알고리즘을 옮겼다.
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
| `src/core/85_skeleton.js` | Analyze Skeleton (2D/3D) — 동작 명세·블랙박스 대조로 새로 쓴 clean-room 구현 (가지치기 포함) |
| `cleanroom/` | 그 명세(SPEC_analyze_skeleton.md)·증거 기록(NOTES.md)·실험 스크립트 |
| `test/test_skeleton.js`, `ref/askel/SkelDriver.java` | Fiji AnalyzeSkeleton(GPL-3, 시험에서만 실행)과 고정 102 경우 + FUZZ 대조 |
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
| `test/test_ij_edge.js` | NaN·±∞·일정한 이미지 경계 사례를 ImageJ 를 직접 실행해 대조 19 항목 |
| `ref/gen_audit.py` | 파일 읽기 회귀 시험 파일(`test/io/audit`) |
| `src/core/75_stackreg.js` | StackReg/TurboReg 엔진 (pystackreg C++ 이식) |
| `test/stackreg_parity.py`, `test/test_stackreg.js` | pystackreg 대조 45 사례 |

화소 표현: plane = `{w, h, bit: 8|16|32, data: Float32Array}`. 8/16-bit 도 정수를 Float32 에 담고, ImageJ 와 같은 규칙 `(int)(v+0.5)` 후 범위 자르기로 다시 정수화한다.
이미지 = planes 배열, 순서 c + nc·(z + nz·t) (ImageJ hyperstack 과 같음). RGB 는 8-bit plane 3 장 + `rgb` 표시.

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
2. **Phase correlation 부화소**: 이동하는 쪽 이미지의 Hann 창을 추정 이동만큼 옮겨 다시 계산(3 회). 고정 창은 0 쪽으로 최대 0.36 px 끌린다(합성 시험). 결과 최대 오차 0.013 px.
3. **Orientation 기울기**: 4차 중앙 차분. 2차 중앙 차분은 주기 9 px 결에서 1° 가량 기운다.
4. ImageJ 의 Java `int` 넘침(Li·Huang 의 합, 2³¹ 초과)은 재현하지 않는다(double 로 셈).
5. **Mean/Variance + NaN**: ImageJ 는 이동 합이 NaN 을 지난 뒤 줄의 다음 화소 하나를 더 NaN 으로 남긴다(전체 재계산 전). 여기서는 커널 안에 NaN 이 있는 화소만 NaN — 유한값은 ImageJ 와 같음.
6. 16-bit 의 일부 결과(rolling ball shrink 경로·bandpass)는 ImageJ 의 float32 중간 계산과 ±1 count 차이가 0.5 % 미만 화소에서 난다(전체 눈금의 1.5×10⁻⁵).

## 3b. 현미경 형식 (ND2 · CZI)
- **ND2** (NIS-Elements): chunk map → `ImageAttributesLV!`·`ImageMetadataSeqLV|0!`·`ImageMetadataLV!`(v3, CLx Lite Variant 이진) 또는 같은 이름의 XML(v2.x). 해석 규칙은 Python `nd2` 패키지 0.12 의 modern reader 와 같다:
  실험 loop 은 바깥 loop 먼저(T·NETime·XY·Z, Spect loop 은 건너뜀, 같은 층의 같은 종류 loop 은 큰 쪽), 프레임 = `ImageDataSeq|i!` 자료의 8 바이트(시각) 뒤, 행 간격 = `uiWidthBytes`,
  component 가 채널당 3 이면 RGB, 아니면 모든 component 가 채널. lossless = zlib. 없는 프레임은 0(알림). XY 위치 = 이미지 따로. 시간형 loop 이 둘 이상이면 하나의 t 로 잇는다(알림).
  v1.0(JPEG 2000)·lossy 는 분명한 오류로 알린다.
- **CZI** (ZEN): ZISRAWFILE → directory(DV 항목) → subblock. 층 0 만(저장 크기 = 논리 크기, pyramid type 0). 화소 형식 Gray8/16/32f, Gray32(int)·Gray64(→ float, 알림), Bgr24(→ RGB), Bgr48·Bgr96f(→ R·G·B 3 채널), Bgra32.
  압축: 없음·zstd0·zstd1(첫 바이트 = 머리 크기, 항목 1 = HiLoByteUnpack: 낮은 바이트 전부 → 높은 바이트 전부)·JPG(브라우저). JPEG XR·LZW 는 오류로 알림.
  mosaic: scene 의 층 0 subblock 경계 상자 위에 M 순서로 그려 뒤 tile 이 위(libCZI 와 같은 결과 — 실측 겹침 파일 포함 일치). scene·그 밖의 차원(H·R·I·V·B 중 값이 둘 이상) = 이미지 따로.
  Scaling·채널 이름·색은 `ImageDocument/Metadata/Scaling`, `…/Information/Image/Dimensions/Channels`, `…/DisplaySetting/Channels` 자리만 읽는다(같은 이름의 다른 요소 무시 — LLS7 파일은 다른 곳에 다른 z 값이 있음).
- 채널 색 → 가까운 LUT(Red·Green·Blue·Cyan·Magenta·Yellow). 광학 이미지 보정 단위는 µm(10 nm 미만이면 nm).
- 큰 파일: 위치/scene 이 6 개를 넘거나 32-bit 메모리로 600 MB 를 넘으면 고르기 창(위치·scene, Z·T 범위와 간격). 한도 1.6 GB. 간격을 두고 고른 Z 는 z step 도 그만큼 곱한다.

## 3b-2. Leica LIF (v1.4)
- 구조: [0x70, 크기, 0x2A, 글자 수, XML(UTF-16LE)] 다음 memory block 들(v1 크기 int32 · v2 int64). XML 의 Element 마다 `Data/Image/ImageDescription` 이 이미지 하나(문서 순서, Memory Size > 0 인 것만).
- 차원 배치(Bio-Formats LIFReader 와 같음): 앞 두 차원 = 이미지 면(XZ 주사면이면 세로 = Z), 남은 공간 차원 하나 = Z, 시간·파장(T·λ·Λ·그 밖)은 하나의 T 로 이음(앞 차원이 빠르게), tile(M)마다 이미지 따로.
  Z 가 1 장이고 λ 등이 있으면 Bio-Formats 는 그것을 Z 로 두지만 여기서는 T(화소 같음, Image Info 에 적음).
- 화소 크기 = Length/(N−1) (µm, 10 nm 미만이면 nm). 세로 축이 공간이 아니면(XT·X-λ) 화소 높이 = 너비 + 알림, 세로 간격은 Image Info. z step = ScannerSettingRecord dblStepSize(있으면, Bio-Formats 와 같음).
- RGB(색 카메라, X BytesInc = 3·bps): Leica 는 B·G·R 로 저장 → ChannelTag(1 R·2 G·3 B, 없으면 BytesInc 거꾸로)로 R·G·B. 8-bit = RGB 이미지, 16-bit = R·G·B 3 채널(ND2·CZI 와 같은 모형).
- 채널 이름(liffile 과 같은 차례): ChannelProperty DyeName('Leica/' 뗌) → 모두 없으면 HardwareSetting 의 MultiBand DyeName(채널 수와 같을 때) → 순차 주사 단계마다 하나 → C1…. LUT 이름 → 색(모르는 이름은 그 채널만 비움).
- 손상: memory block 없음·채널 자료형 섞임·차원 값(N<1, BytesInc≤0)·크기 초과인 series 는 알림과 함께 건너뛰고 나머지를 읽는다. 파일 끝이 잘린 series 는 넘는 화소를 0 으로 채우고 알림(Bio-Formats 와 같은 방침).
  XML 머리·memory block 머리 안에서 잘림 → 분명한 오류/알림. 잘못된 문자 참조(&#x110000;)는 U+FFFD.
- 고르기 창: series 마다 크기가 다르면 Z·T 는 가장 큰 값으로 묻고 series 마다 그 안의 것만(메모리 계산도 series 별). 잎 이름이 겹치면 경로(맨 위 요소 뺌)로 이름. metaText = 첫 이미지 요소의 XML.

## 3b-3. Olympus/Evident OIR·POIR (v1.5 — 공개 실측 13 파일로 검증)
- 형식 지식: oirfile(BSD-3, © C. Gohlke — 머리말에 고지). Bio-Formats 7.3.1 은 블랙박스로만 실행.
- **실측 검증**(`ref/real/oir`, OME 시험 자료 서버의 FV4000 공개 파일 13 개, CC BY 4.0 — 저장소·배포에 올리지 않음, `ref/gen_oir_real.py`):
  화상이 있는 12 파일 모든 면이 oirfile = Bio-Formats(λ 파일 제외) = 앱 (면별 SHA-1), 그중 3 파일은 함께 배포된 OME-TIFF 와도 같음. 크기·채널 이름·화소 크기·z step(10⁻⁶)·frame interval(10⁻⁴)·참조 이미지 일치.
  λ 파일은 Bio-Formats 가 λ1 면을 모든 칸에 되풀이(결함) → oirfile 을 따르고 파일의 LAMBDA 위치(122/122 프레임)로 확인. Map(지도) 파일 = 화소 없음 → 타일 파일을 열라는 분명한 오류.
- 실측에서 고친 것: ZSTACK 축 정의가 둘(쓰인 것 paramEnable=true + 안 쓰인 'Range')인 파일에서 엉뚱한 쪽을 써 'Z 4/21' 거짓 알림 → paramEnable 우선(oirfile 은 마지막 것이라 틀림).
  frame 이름(l001z002…)으로 (t, λ) 안의 z 간격·t 첫 프레임 간격, LUT 앞의 채널 GUID 로 채널 색, UID 머리의 면 안 위치·길이로 블록 배치, 잘린 파일 차례 읽기(BMP·NULL 블록·UID 검사).
- 차원: C·Z·T·λ(λ 는 T 에 이음, λ 가 빠르게). 4 바이트 깊이 = float32, 값이 float 로 맞지 않으면 uint32 + 알림. 참조 이미지(REF_LSM0, 파일 맨 앞) = 둘째 이미지 — 보정은 본 이미지 값을 붙이고 알림.
- 손상·보조: 빠진 면 = 0 + 알림, 설정보다 면이 적음 → 알림(보조 파일 `_00001…` 은 읽지 않음), 보조 파일만 열면 분명한 오류.
- POIR(ZIP): stored·deflate, data descriptor, UTF-8 이름. ZIP64·암호화는 분명한 오류. 메모리 한도는 묶음 전체.
- 실측 미검증(합성만): 선 주사, RGB/COLOR2D, 4 바이트 깊이, 보조 파일, POIR, 참조 이미지 보정, 카메라 이미지. 시험: `node test/test_oir.js` (70, 실측 15 포함).

## 3c. Registration (StackReg)
- 엔진 `75_stackreg.js` = pystackreg 0.2.8 C++(TurboReg ANSI-C 이식, 허용 라이선스 — 머리말에 원문·논문 인용) 를 연산 순서까지 옮김. Java TurboReg/StackReg(재배포 제한 라이선스)는 옮기지 않았다.
- pystackreg 대조 45 사례(5 변환 × 크기 4, 음수·8-bit, 512², stack previous/first/mean): 행렬 최대 차 3×10⁻¹³(rigid — Node 와 glibc 의 atan2·sin 1 ulp 차이까지 추적), 나머지 0, 변환 이미지 화소 동일.
- 정답 변환 회복(합성): translation 0.005 px, rigid 중심 0.007 px·각도 0.02°, scaled 0.004 px·0.009°, affine 모서리 0.03 px.
- pystackreg 의 특징을 그대로 따른 것: 맞춤 전 마지막 행·열을 뺌(`ref[:-1,:-1]`), Levenberg–Marquardt 단계에서 Hessian 대각만 갱신(원래 Java TurboReg 와 다름 — 위 정확도는 이 방식의 값).
- pystackreg 와 다르게 한 것 하나: 기준 '첫 n 장 평균'(n>1)이면 첫 장도 그 평균에 맞춘다(pystackreg 는 첫 장을 그대로 둬 첫 장만 어긋남). 대조 시험도 이 규칙으로.
- 매크로: '이전 장' 기준의 기준 장(anchor)을 기록해 재생·batch 에서 같은 결과.
- 앱 규칙: 기준 채널로 변환을 구해 모든 채널에 적용. z·t 가 모두 있으면 t 마다 z 최대 투영으로 구해 모든 z 에 같은 변환. 결과는 원래 bit 로 반올림(8/16-bit). 이미지 밖은 0.
  보고값: drift = 이미지 중심이 기준 대비 옮겨진 양(px), 회전 = 화면에서 반시계 +(Orientation 과 같은 규약), scale = √|det|. ImageJ 매크로로는 `run("StackReg ", "transformation=[Rigid Body]")`.
- 속도(512², node): rigid 0.17 s, affine 0.2 s / 장.

## 3d. AFM·SPM 형식 (v1.3)
- 해석 규칙의 기준은 Gwyddion 2.64 의 가져오기 모듈(modules/file/nanoscope.c·gwyfile.c·gsffile.c·ardf.c). 코드는 옮기지 않고 동작을 맞췄으며, 시험 정답은 그 모듈이 같은 파일을 읽은 값이다.
  `gwyddion --convert-to-gwy` 는 GTK 를 초기화하지 않아 ARDF 모듈이 실패하므로, 모듈을 직접 부르는 작은 변환기(`ref/gwyconv`)로 GWY 를 쓰고 gwyfile(독립 파이썬 파서)로 읽는다.
- **Asylum IBW (Cypher·MFP-3D, AR 소프트웨어 Igor)**: 한 wave 에 여러 채널, wave 단위는 하나(m) → 채널 이름([Nap]<종류>(Trace|Retrace))으로 단위를 정한다:
  Height·ZSensor·Amplitude·Deflection = m(→ nm), Phase = deg, Potential·Voltage·Bias·UserIn = V, Current = A(→ pA·nA), Frequency = Hz, Capacitance = F. 모르는 종류는 단위 없이 + 알림.
  연구실 실측 Cypher 5 파일(AC topography·KPFM(SKPM)·C-AFM·EFM·FM-EFM, AR 16.23·18.40) 48 채널: 화소(방향 포함)·단위·화소 크기 = Gwyddion 2.64 (`ref/real/ibw`, 저장소·배포에 올리지 않음).
  .ARDF 는 AR 의 Force Map(힘 곡선 지도) 전용 — 이미지 주사는 .ibw 로 저장된다.
- **값 단위**: 이미지 값은 물리량으로 바꾼다. 길이(m)는 늘 nm, N·Pa·A·F·S 는 중앙값 크기에 맞는 접두어(nN·MPa 등), 그 밖(V·deg·Hz·단위 없음)은 그대로. 배율은 정보 칸에 적고, 마우스 아래 값에 단위를 붙인다.
  같은 물리량을 남기는 연산(필터·회전·자르기·정합·평균/최대 투영 등, `KEEP_UNIT`)만 단위를 이어받고, 형 변환의 눈금 바꾸기·문턱·FFT·곱셈 등은 지운다. 되돌리면 단위도 돌아온다.
- **Nanoscope**: 머리 길이 = 앞 8 필드 안 `\Data length`. 값 = 원시 정수 / 256^(Bytes/pixel) × `@2:Z scale`(또는 `@4:`)의 hard value × soft scale(`@Sens. …`, 단위 곱 — 예 V × nm/V = nm).
  soft scale 에 단위가 없으면 V(Phase 는 deg). 9.2(0x09200000) 이후는 자료가 32-bit(Bytes/pixel 은 눈금만). 이미지 절 크기가 틀리면 Ciao scan list 의 전역 크기, 종횡비는 Gwyddion 의 두 규칙.
  첫 줄이 아래(위아래 뒤집기). 채널 이름은 따옴표 안 화면 이름(Gwyddion 은 soft 이름). 힘 곡선(Force)·단면 묶음·텍스트형은 분명한 오류, Force Volume 은 이미지 채널만(알림). `\*Force file list` 로 시작해도 Operating mode 가 정함.
- **GWY**: 13 성분 형식, GWYO(1.x) 의 옛 맨 위 컨테이너(GType 번호), DataField·mask(→ 8-bit 0/255 이미지)·Brick(→ z-stack, z 가 길이면 z step), meta(GWYO 는 `/meta/…`). 단위 문자열의 접두어는 Gwyddion 처럼 버림(값은 늘 기본 SI).
  그래프·XYZ·curve map·spectra 는 알리고 읽지 않음. **GSF**: 머리 + NUL 1–4(4 바이트 경계) + float32, 크기가 정확히 맞아야 함, NaN 은 그대로(알림).
- **ARDF**: 블록 머리 CRC32 를 모두 확인(손상은 오류), FTOC → IMAG(TTOC·IDEF·IBOX·IDAT·GAMI). IDEF 의 간격·단위·이름을 쓴다(Gwyddion 은 FastScanSize 로 정사각 가정 — 직사각이면 틀림).
  이미지는 첫 줄이 아래(IBW 와 같음; magic-afm 도 origin=lower — Gwyddion ARDF 모듈만 뒤집지 않음). 0 으로 채운 IDAT 줄(중단된 주사)은 NaN + 알림. 힘 곡선(VOLM)은 크기·채널만 알림.
  **실측 ARDF 는 구하지 못했다**(공개 저장소에 없음) — magic-afm·Gwyddion 이 함께 받아들이는 합성 파일 4 개로만 확인.
- 실측: SciFiDatasets(pycroscopy)·AFMReader·TopoStats 시험 자료 13 파일(Nanoscope 5.31·8.15·9.2·10, GWYP 3·GWYO 2, GSF, 힘 곡선 1) — 모두 Gwyddion 결과와 값·크기·단위 일치(상대 1.5×10⁻⁷, float32 반올림).

## 3e. Color Threshold (v1.3)
- ImageJ `ColorThresholder` 와 같은 색 변환: HSB(java.awt.Color.RGBtoHSB 의 float 연산 순서 그대로), Lab(ColorSpaceConverter D65, L×2.55·a/b×1.0625+128), YUV(Landini, byte 로 감음).
  ImageJ 를 실행해 **2²⁴ 색 전부** 같은 바이트임을 확인. 문턱 = 채널마다 pass(안)·stop(밖), 세 조건 AND → 8-bit mask(255). Auto = HSB 는 Brightness, RGB 는 세 채널, Lab·YUV 는 첫 채널(autoSetThreshold).
  Sample = ROI 안 범위, HSB Hue 는 쓰인 범위와 가장 긴 빈 구간 중 좁은 쪽(빨강처럼 0/255 를 감싸면 stop).
- 결과는 새 mask 이미지(기본) 또는 제자리(스택은 '모든 장' 필요). .ijm 은 ImageJ 가 만드는 매크로와 같은 단계(변수 이름 ct_ — batch 고리의 i 를 덮지 않음); ImageJ 매크로는 한 장만.

## 3f. Analyze Skeleton (v1.4) — clean-room 독립 구현
- **출처·라이선스**: `src/core/85_skeleton.js` 는 Fiji `AnalyzeSkeleton_`(Arganda-Carreras 2010, GPL-3) 의 **소스를 보지 않은** 작업자가
  동작 명세(`cleanroom/SPEC_analyze_skeleton.md` — 출력이 무엇인지만 적음)와 jar 를 블랙박스로 실행한 결과만으로 새로 쓴 구현이다(clean-room).
  명세에 없던 규칙(가지치기 전부, LSP 대각 항목)은 실험으로 찾아 `cleanroom/NOTES.md` 에 입력·출력 증거를 남겼다. 명세 작성자(소스를 본 쪽)와 구현자를 나눴다.
  이전(v1.3 개발 중) 의 GPL 절차 이식본은 배포 전에 지웠다. 시험에서만 Fiji jar 를 실행한다(앱에 넣지 않음). ※ 법률 자문이 아님 — 외부 배포 전 기관 확인 권장.
- **결과 = Fiji**: 나무 수·모든 열·가지 목록(길이·두 꼭짓점·running average·줄기 수)·태그 이미지·Longest shortest path 이미지(경로 96) 이 같음(상대 10⁻¹²).
  `node test/test_skeleton.js` 고정 102 경우(가지치기 52 포함) + `FUZZ=3000` 씨앗 여럿(각 5852, 0 실패), 구현자의 별도 가지치기 무작위 9300 경우 0 불일치.
- **Fiji 보다 나은 점**: Longest shortest path 크기 제한 없음 — 꼭짓점 ≤ 600 은 Fiji 와 같은 행렬 절차, 그 위는 block-cut tree 이심률·이심률 경계 + Dijkstra 로 같은 쌍·경로·시작점
  (Fiji 실행과 꼭짓점 1232–8494 나무·고리 골격 대조 같음). 길이가 같은(10⁻¹² 안) 쌍·경로가 있어 Fiji 의 갱신 순서에 따라 갈릴 수 있으면 경고(길이는 같음).
  속도: 1500² 꽉 찬 덩어리 0.9 s, 2000² 뱀 경로 + LSP 0.6 s, 256³ 희소 골격 1 s, 꼭짓점 52 k 나무 0.9 s. 부피 한도 2.5×10⁸ 화소(분명한 오류).
- **가지치기**(Fiji 와 같음): Prune cycle method = none · shortest branch(줄기 수가 가장 적은 가지, 가운데 줄기에서 자름) · lowest intensity voxel · lowest intensity branch
  (세기 = 원본 이미지의 화소 주변 상자 평균 — Fiji 와 같음; 원본은 같은 크기의 회색 이미지, 8-bit 로 대조). Prune ends(끝점 가지 제거). 고리를 자르면 골격 전체를 다시 분석.
- 입력은 0/255 이진(0 아닌 화소 = 골격). '먼저 Skeletonize' 는 2D(ImageJ BinaryProcessor). 스택은 3D(z 간격 = z step, 보정 단위로). 결과: Skeleton 표(골격마다)·Branch information(나무마다 길이 내림차순)·Summary 한 줄,
  Tagged skeleton(끝 30·줄기 127·접점 70, Fire LUT), 선택하면 Longest shortest path 이미지(경로 96). intensity 열은 없다(Fiji 의 intensity 열은 골격 이미지의 값이라 늘 255).
- 표 상한: Skeleton 20 000 줄, Branch information 100 000 줄(넘는 것은 경고와 함께 넣지 않음). Slice 열 = c·t. z 간격 모름 → 1 화소 + 경고. hyperstack → 현재 c·t 의 z 스택만 + 경고(Fiji 는 모든 장을 한 3D 로 봄).
- .ijm: Skeletonize·가지치기는 Fiji 에서 입력을 바꾸므로 `Duplicate...` 사본에서 → `run("Analyze Skeleton (2D/3D)", "prune=[…] prune_0 calculate show original_image=[…]")` → 원본 다시 선택.
  `display`(labeled skeletons) 는 앱에 없어 넣지 않음. Fiji 플러그인을 넣은 ImageJ 로 실행해 같은 표 확인(가지치기·original_image 포함).

## 3g. 속도 개선 (v1.4) — 결과는 ImageJ 와 화소 하나까지 같음
- **Median 16/32-bit, 창 화소 n ≥ 25** (반지름 ≥ 2.5): 값을 순위로 바꿔(16-bit 정수는 값 그대로, 그 밖은 정렬한 서로 다른 값의 번호, ±∞ 포함) 세 단계 개수(순위·64·4096 칸) 이동 히스토그램.
  NaN 없는 원형 창은 개수가 늘 홀수 → 중앙값이 하나로 정해져 ImageJ getNaNAwareMedian 과 같다. 창에 NaN 이 있는 화소만 ImageJ 규칙(가운데 화소 추정, NaN 뺀 개수)으로 직접 계산.
  1024² r=15: 16-bit 13 s → 0.9 s, 32-bit 14 s → 2.1 s. 2048² 32-bit r=10 7.6 s (서로 다른 값 4 M 개, 최악).
- **Remove Outliers**: ImageJ RankFilters.filterLine(OUTLIERS) 의 순서 그대로(float 연산 v·s+thr < 창 최댓값 일 때만 중앙값, 추정값 = 그 줄에서 마지막으로 구한 중앙값, 줄 처음 = 첫 화소(NaN→0),
  getMedian 은 NaN 을 빼지 않음). 이전 구현은 NaN 이 있는 32-bit 에서 80 화소가 ImageJ 와 달랐다(NaN 창의 짝수 개 중앙값이 추정값에 따라 달라짐) → 같음.
- 시험(`test/test_ij_edge.js`, ImageJ 직접 실행): 16-bit(넓은 범위·같은 값 많음)·32-bit(NaN 덩어리·±∞·같은 값 많음 + NaN) × r 2.5·5·8.5 Median, × r 6·1, 정수·소수 문턱, Bright·Dark Remove Outliers — 모두 차이 0.
  빠른 경로 = 직접 경로 무작위 400 경우(NaN·±∞·같은 값) 차이 0.
- **Wand**: Worker 에서(형식 배열 스택). 4096² 한 덩어리 0.75 s, 그동안 화면 멈추지 않음. 연속 클릭은 앞 것이 끝날 때까지 알림. 결과 = 이전 주 스레드 구현(무작위 3000 경우 같음), ImageJ doWand(8-connected) 꼭짓점과 같음(Chromium 시험).

## 3h. Drop Shape — 접촉각 · 펜던트 드랍 (v1.6)
ImageJ 플러그인 DropSnake·LB-ADSA(Stalder 2006·2010)·Pendent_Drop 이 하던 일을 **Analyze** 패널 맨 아래 칸 **Drop Shape** 에(이미지 처리의 연장 — 따로 탭을 두지 않음, 패널 머리의 "Drop Shape ↓" 로 바로 감). 계산 `src/core/88_drop.js`(Worker), 화면 `src/ui/35_drop.js`.
- **윤곽**: 가우스(σ 1) → Otsu 문턱 → 방울 성분(영역 테두리 중앙값보다 어두운/밝은 쪽 자동, 구멍 메움 — 방울 안 반사점) → marching squares 부화소 등고선 → 법선 ±3 px 에서 경사 최대 자리로 다듬음(Canny 와 같은 정의, 비네팅으로 바탕 밝기가 고르지 않을 때 문턱 등고선이 밀리는 것을 없앰). 영역 테두리를 따라 생긴 가짜 윤곽은 맞춤에서 뺀다.
- **Young–Laplace**(무차원, apex 곡률 반지름 b 로 나눔, Rotenberg·Boyce·Neumann 1983 ADSA): dx/ds = cos φ, dz/ds = sin φ, dφ/ds = 2 ± Bo·z − sin φ/x (세실 +, 펜던트 −), Bo = Δρ g b²/γ. RK4(ds 0.002–0.003). 맞춤 = 윤곽 점과 이론 곡선의 최단 거리 제곱합, Levenberg–Marquardt.
- **펜던트**: 영역 위 테두리에 닿은 성분 → 위쪽 폭이 일정한 구간 = 바늘(두 가장자리 직선 맞춤 → 지름·기울기) → 그 아래 윤곽에 [apex x·y, b, Bo, 기울기 ω] 맞춤. 바늘이 붙는 높이 = 적도 위에서 이론 곡선 반지름 = 바늘 반지름인 곳. γ = Δρ g b²/Bo, 부피·겉넓이 = 회전체 적분(바늘 끝까지), Worthington Wo = Δρ g V/(π γ D_바늘)(Berry et al. 2015). 눈금 = 바늘 바깥 지름(mm) ÷ 바늘 화소 폭, 없으면 이미지 보정.
- **접촉각** — 기준선: 자동(윤곽 열마다 가장 위 점 → RANSAC 직선 = 기판 윗면) 또는 선 ROI(Use line ROI). 방법 4 가지, 모두 액체 쪽에서 잰 각:
  Young–Laplace(축대칭, 중력 포함 — 기준), 접선(접촉점 둘레 반지름 0.12·min(높이, 반폭) 안 윤곽을 호 길이 u 의 2 차식 x(u)·t(u) 로 맞춰 기준선까지 외삽 — θ > 90° 에서도 한 값; 기준선 쪽 3 px 은 흐림 구간이라 뺌; 문턱 등고선 사용), 타원(Fitzgibbon 1999 직접 맞춤, Halir–Flusser), 원(구면 모자 — Bo ≈ 0 에서만).
  바늘이 방울 안에 꽂혀 있으면(영역 위 테두리에 닿음) 바늘 끝보다 위 윤곽은 맞춤에서 뺀다.
- **경고**: 방법 사이 > 5°, 접선 왼·오른 차 > 5°, θ > 150° 접선, 잔차 rms > 1 px, Wo < 0.5, 기울기 > 5°, 눈금 없음.
- 불확도: 표의 'fit SD' 는 맞춤 공분산(잔차 분산 × (JᵀJ)⁻¹, b·Bo 상관 포함 gᵀCg)에서 온 **맞춤 통계만**이다. 펜던트는 바늘 가장자리 ±0.3 px 가 눈금에 주는 몫을 'needle-scale SD' 로 따로 보인다(보통 맞춤 SD 보다 10–20 배 큼). 바늘 지름·밀도·g·가장자리 정의의 계통 오차는 들어 있지 않다 — 같은 방울 여러 장(스택 Frames: All)의 흩어짐이 더 정직한 값.
- 결과: 아래 **Drop** 표(방법마다 한 줄, CSV·Copy), 스택은 시계열 그래프(Plots), 매크로 기록(op `drop` — 영역·기준선 좌표째).
- 검증(`node test/test_drop.js` 21, `python3 test/ui_drop.py` 38):
  ① 합성(정답 = `ref/gen_drop.py`, scipy DOP853 — 이 구현과 다른 적분기): 펜던트 16 장(Bo 0.1–0.45, γ 25–75 mN/m, 기울임 ±1.5°) γ 중앙 오차 0.06 %·최대 0.37 %, 부피 ≤ 0.53 %, 기울기 ≤ 0.024°. 세실 20 장(θ 20–160°, Bo 0.0005–0.6, 기울임 ±3°) YL ≤ 0.85°, 접선 θ ≤ 150° 에서 ≤ 1.40°·160° 에서 3.45°, 원(Bo ≈ 0) ≤ 0.15°. 바늘 꽂힌 세실 3 장 YL ≤ 0.08°.
  ② 실측(OpenDrop 예제, `sh ref/fetch_drop.sh` — GPL-3 이라 저장소·배포 제외): water_in_air001 은 OpenDrop 결과 화면과 γ −1.17 %, V −0.87 %, A −0.56 %, Wo +0.34 %, Bo +0.50 %. 같은 물 방울 5 장 평균 γ 73.07 mN/m(물 20–25 °C 의 IAPWS R1-76 값 72.0–72.8 mN/m 근처; 온도 미상). drop_on_surface 접선 118.8° / 123.8° ↔ OpenDrop 120.0° / 124.2°.
  Worthington 수가 작은 장(004: 0.22, 005: 0.10)은 γ 가 ±2 % 흔들린다 — Berry 2015 가 말한 대로이며 화면 경고로 알린다.

## 4. 화소 크기 (calibration) 규칙
- 출처가 확실한 값만 자동 적용: Nikon ND2 `dCalibration`(bCalibrated 일 때), Zeiss CZI `Scaling/Distance`, ImageJ ImageDescription(unit)+XResolution, OME-XML PhysicalSize, FEI/Thermo 34680·34682 `[Scan] PixelWidth`, Zeiss 34118 `Image Pixel Size`, TESCAN 50431, Hitachi 옆 `.txt` `PixelSize`(nm), JEOL 옆 `.txt` 축척 막대, DM·MRC·EMD·BCF(앱 기존 파서), Igor wave scaling, Nanonis SCAN_RANGE, Gwyddion xreal/xres, Nanoscope Scan Size(전역 크기·종횡비 규칙), ARDF IDEF 간격.
- TIFF/PNG 해상도 태그만 있으면 '후보'로 보이고 적용하지 않는다(72 dpi 같은 기본값이 흔함).
- FEI 데이터 막대(`[Image] ResolutionY` < 이미지 높이)는 경고와 한 번에 잘라내기.
- 크기 필터의 단위 변환: size(unit²)/(pw·ph) = 화소 수. Circularity 필터는 화소 단위 둘레·넓이(ImageJ 와 같음 — 보정과 무관).

## 5. 검증
- `node test/test_img.js` (378): 아래에 더해 파일 읽기 회귀 13(`test/io/audit`). ① ImageJ 1.54u 실행 결과 275 사례(8·16·32-bit 필터·문턱 16×2×3 histogram 모양 8 종·이진·maxima 위치·입자 33 열·ROI 8 종(이미지 밖으로 나간 ROI 포함)·계산기·투영) ② 정답을 아는 합성 이미지(줄무늬 방향 5, 부화소 이동 4, 고리 반지름, EDM 대 완전 탐색 거리, 원 ECD·원형도, 로그정규 회복, local Mean) ③ 파일 38(tifffile·pypng·numpy·igor2 로 만든 정답) + TIFF 쓰기→읽기(아주 작은 화소 크기, RowsPerStrip 기본값), 저장한 TIFF 를 ImageJ 가 보정·hyperstack 째 읽는 것 확인.
- 현미경 형식: 합성 nd2 12·czi 9 (`node test/test_img.js` 안, 정답 = nd2 패키지·czifile·libCZI), 실측 nd2 13·czi 28 → 135 대조 (`node test/test_real.js`, 정답 = nd2 패키지·libCZI, Bgr96Float 만 czifile).
- `python3 test/ui_img.py` (41): nd2·czi 열기(화소 = 정답, 고르기 창, 창이 떠 있을 때 또 열기), Register rigid(정답 회전·drift, 화소 = pystackreg), Worker 결과 = ImageJ, 되돌리기, Preview, Threshold 패널·Convert to Mask = ImageJ, Analyze Particles = ImageJ, 마우스 ROI·Ctrl+M = numpy, Find Maxima 개수, Macro 기록→다른 이미지, .ijm, 새로 고침 복원, 모바일, 콘솔 오류 0.
- AFM: `node test/test_afm.js` (164) — 합성 16 파일(정답 = Gwyddion·gwyfile·magic-afm·손 계산) + 실측 13 파일(`ref/real/afm`, 정답 = Gwyddion) + 손상 파일 오류·단위 해석. 정답 다시 만들기: `python3 -I ref/gen_afm.py test/io/afm ref/real/afm` (gwyddion·libgwyddion20-dev·gwyfile·magic-afm 필요).
- Color Threshold: `node test/test_color.js` (17, ImageJ 실행) — 2²⁴ 색 변환, 문턱 7 경우, Auto·Sample, RGB hyperstack TIFF. Chromium: 패널·Sample·Auto, 스택·hyperstack .ijm 을 ImageJ 로 실행해 대조.
- v1.4: LIF `node test/test_lif.js` (47: 실측 9 파일 = Bio-Formats·liffile, 합성 경계 11 = `ref/gen_lif_synth.py`), OIR `node test/test_oir.js` (54), Skeleton `node test/test_skeleton.js` (102, `FUZZ=3000` 5852),
  ImageJ 경계 `node test/test_ij_edge.js` (47: + 큰 창 Median·Remove Outliers 28), Chromium + Wand(ImageJ doWand)·Skeleton 가지치기·.ijm(Fiji 실행).
- ImageJ 기준 다시 만들기: ImageJ 소스(github.com/imagej/ImageJ) 를 javac 로 컴파일 → `python3 -I ref/gen_ref.py ij.jar ref/data`.

## 6. 독립 검토로 고친 것 (v1.0 안)
되돌리기가 보정을 지우고 엉뚱한 매크로 단계를 빼던 것(단계별 표지로) · Preview 계산 중 Apply 가 무시되던 것(전역 표시 → 매개변수) · 적용 뒤 Preview 가 한 번 더 겹쳐 보이던 것 ·
RGB Z Project 의 Sum/SD 가 깨진 이미지 · hyperstack Z Project 가 시간까지 합치던 것 · OME DimensionOrder 무시(채널·z 뒤섞임) · RowsPerStrip 기본값(2³²−1) 에서 멈춤 ·
TIFF 쓰기 XResolution 넘침 · 열 단위가 지금 고른 이미지를 따르던 것(측정 때 단위를 행에 저장) · 분포 통계가 0 이하 값을 조용히 빼던 것(선형 = 전부, 로그정규 = 양수, 알림) ·
열 때 표시 범위(ImageJ 와 같은 min–max) · 형 변환 범위를 매크로에 기록 · Duplicate 기록 · .ijm 에 문턱 명령 · FFT radial profile 중심 반 화소 · 이미지 밖으로 나간 ROI 의 둘레·Feret(ImageJ 와 같음) ·
batch 의 redirect 를 원본 상대 지정(__parent/__root)으로 · RGB 입자 분석 막음 · Space 고정 · Flip·Register 뒤 남던 겹쳐 그림 · 보정 단위 Gaussian 의 세로 σ.

## 6b. 독립 검토로 고친 것 (v1.1)
ND2 v2 XML 의 4096 자 넘는 태그(긴 Description)에서 열기 실패 → 따옴표를 지키며 태그 끝까지 읽기, `&#..;` 해독 · 16-bit RGB ND2 를 8-bit RGB 로 잘라 버림 → R·G·B 3 채널 16-bit ·
TIFF spacing 단위(nm 보정에서 z step 1000 배) · CZI 의 RGB 채널 + 회색 채널 섞임 → 따로 이미지 · CZI 고르기 창이 첫 scene 크기로 메모리 계산, 전체 한도 없음 → scene 별 크기·고른 것 합으로 한도 ·
'첫 n 장 평균' 기준에서 첫 장 미정렬 · Register 매크로 재생 시 기준 장 바뀜 · 위치 64 개 넘을 때 기본 간격 계산 · 고르기 창이 떠 있을 때 또 열면 앞 파일이 사라짐 → 차례로 ·
채널 수가 바뀐 뒤 남는 채널 이름 · 여러 scene 이미지마다 메타데이터 원문 복사(저장 공간) → 첫 이미지에만 · nd2 sub-loop 규칙 · CZI JPG 색 공간 변환 끔.

## 6c. 전수 검수 (v1.2) — 독립 검토 4 영역에서 확인된 결함 수정
- **계산 엔진(ImageJ 직접 실행으로 대조)**: NaN 이 있는 32-bit 에서 Find Maxima 봉우리 누락(정렬 키 (int)NaN), Max/Min·Open/Close·Top Hat 이 NaN 근처에서 틀림, Mean/Variance 의 NaN 이 줄 끝까지 번짐,
  Median 의 NaN 처리(ImageJ getNaNAwareMedian: 가운데 화소 추정·아래쪽 중앙값), 32-bit Image Calculator Divide 가 정수로 잘림, 측정·히스토그램·표시 범위에 ±∞ 포함(ImageJ 는 뺌),
  min=max 에서 8-bit 로 NaN, Z Project 의 NaN(평균·최대·최소는 뺌)·16-bit median 은 32-bit, Enhance Contrast 가 일정한 이미지를 0 으로, 임의 각도 회전 가장자리 1 px, 정수·축 평행 다각형 둘레(traced),
  얇은 이미지 Subtract Background → 분명한 오류, 8/16-bit ÷0 은 그대로·Reciprocal 1/0=NaN, Apply LUT 반올림, float 최빈값 칸 폭 1, 빈 histogram 문턱 범위, 일정한 16/32-bit 문턱 없음, 다각형 float32 좌표,
  NaN·∞ 가 있는 이미지의 정렬 → 분명한 오류, StackReg 최소 16 px, 일정한 이미지 Orientation = 방향 없음. 속도: Include holes 입자마다 전체 크기 배열(1024² 3.2 s → 상자 크기), redirect 없을 때 반복 복사.
- **파일**: TIFF·PNG·NPY 메모리 한도와 잘림·손상 → 분명한 오류(0 으로 채우지 않음), LZW 손상 검출, ImageJ 한 채널 표시 범위, multi-series OME → series 마다, OME PhysicalSizeZ, (H,W,3) uint16/float NPY,
  텍스트(';'·',' 소수점, NaN·inf, 줄 끝 구분자, 머리 줄), ImageJ signed 16-bit, 4-bit·1-bit 팔레트, ImageJ 이어진 RGB·부호 스택, int64·uint64, IBW 이름 위치(28)·음의 간격·4D·잘림,
  쓰기 해상도 유리수 정밀도·8/16-bit 반올림·자르기, StripByteCounts 없음, CMYK·YCbCr 알림, UTF-8 µ(태그·옆 파일), 경고 두 번, 파일이 밝힌 단위(micron → µm) 유지, CZI XML 경로 되돌림.
- **상태·매크로**: Scale 등 되돌리기에서 보정이 돌아오지 않던 것(이전 수정에서 주석이 코드를 삼킴 — 주입 전 node --check 추가), 계산 중 되돌리기·닫기·모두 지우기 막기, 계산 중 닫힌 이미지의 화소가 되살아남·고아 화소 삭제,
  중단이 기다리던 일·batch 를 멈추지 못함(중단 세대 번호), .ijm 이 ImageJ 에서 멈추거나 다르게 실행(Flip·Rotate·Scale·Duplicate·Bandpass·Enhance·Convert to Mask·Conversions·redirect·문자열·EDM·Measure ROI·Apply LUT·대화 상자 없는 단일 장 처리·이미지 파일만 열기),
  재생: Orientation 색 지도·Measure ROI·Apply LUT 범위·이진이 아닌 이미지의 Analyze·없는 redirect, batch 파일 실패 시 정리·옆 .txt, Convert to Mask 가 스택을 한 장으로, 저장 실패 후 옛 화소(판 번호 확인),
  복사 TSV·plot CSV 머리, 되돌리기의 composite·databar, Summary 시간축 Slice, n=1 SD, FFT 이미지 역공간 보정(1/unit), 되돌리기 메모리 계산(공유 배열), 매크로 덮어쓰기 확인, IndexedDB 연결 재사용.
- **표시·상호작용**: 큰 이미지에서 Adjust 패널이 다시 그릴 때마다 10 s 멈춤(통계 캐시·형식 배열 정렬), 2.4×10⁷ 화소 넘는 이미지 확대 시 원래 해상도(보이는 영역만)·PNG 내보내기 원래 해상도,
  ±∞ 이미지가 검게, 끝난 일꾼이 큰 배열을 붙잡음, 닫은 이미지의 표시 버퍼, 막대 끌기(한 프레임에 한 번)·16-bit 65536 칸 표, 터치로 다각형, min=max·max<min 표시(ImageJ 와 같은 이진),
  16-bit 표시 매핑 = ImageJ ShortProcessor, 재생 중 fps, 썸네일 갱신, 눈금 막대(보이는 이미지 안·단위 자동·글자 크기 비례), 오른쪽 단추로 ROI 그리지 않음·ROI 좌표 내림(ImageJ), Esc 로 창 닫기, HiLo composite.

## 6d. 독립 검토로 고친 것 (v1.3)
- AFM: `\*Force file list` 로 시작하는 Force Volume 을 힘 곡선으로 거부, 0 으로 채운 ARDF IDAT 줄에서 파일 전체 실패, `Data offset` 이 수가 아니면 머리를 이미지로 읽음,
  Nanoscope E(`Samps/line: 512 512`), Scan Size 두 수 사이 빈칸 여럿, ARDF TEXT 가 24 바이트보다 짧을 때, `º`(0xBA) 단위, 파일 고르기 창의 .100–.999, GWYO meta. — 각 경우 합성 파일과 Gwyddion 대조 시험 추가.
- Color Threshold: 변환 중 Auto·Sample 이 이전 색 공간·다른 이미지의 캐시로 계산(→ 준비될 때까지 끔·확인), 스택 .ijm 이 ImageJ 에서 멈춤·한 장만(→ 장마다 고리·Images to Stack·hyperstack 되살림, ImageJ 로 실행해 화소·차원 대조),
  큰 이미지에서 끌 때 전체 세기(→ 4 M 화소 넘으면 듬성 '≈', 놓으면 정확히), 변환 일이 쌓임(→ 한 번에 하나, 마지막 요청만 이어서, Uint8 로 보냄, Lab 선형화 표), 지금 장 다시 변환(→ 캐시 사용),
  실패 뒤 '변환 중' 고정, 재생 중 겹쳐 그림 사라짐, 닫은 이미지 캐시 남음, min > max 허용(→ ImageJ 처럼 min ≤ max).
- 함께 찾은 결함: Enhance Contrast(표시만) 가 이미지 ver 를 올려 새로 고침 때 그 이미지를 버림. RGB hyperstack TIFF 의 z·t 를 읽고 쓰지 못함(z 로 합침).

## 6e. 독립 검토로 고친 것 (v1.4)
- Analyze Skeleton(이식본 검토): GPL-3 절차 이식이라는 판정 → clean-room 으로 새로 씀(§3f). 그 과정의 반영: 다시 볼 접점 목록 shift·가지 이력 unshift 가 제곱 시간, LSP 한도·경고가 묻힘,
  .ijm 이 원본을 Skeletonize(Fiji 는 제자리), Slice 열에 채널 없음, z 간격 모를 때 조용히 1, LSP 이미지만 골라도 LSP 계산, 표 상한.
- LIF 13 건: 크기 다른 series 섞인 파일에서 큰 series 를 Z·T 일부로 못 엶, series 하나가 손상되면 파일 전체 실패, 잘린 파일, 색 카메라 RGB 가 B·G·R, XT·X-λ 의 세로 보정, 채널 이름·metaText 없음,
  같은 잎 이름, 잘린 머리의 RangeError, 잘못된 문자 참조, 차원 값 검사, LUT 하나 모르면 전부 버림, 진행 표시.
- OIR 21 건(구현 직후 검토): 한 면 블록을 UID 번호로 다시 정렬·잘라 냄, Z drift 파일의 z step, 4 바이트 깊이, fileinfo 만 있는 파일 거부, POIR 메모리 한도·창·취소, 넘치는 바이트·중복 UID, 축 설정 여럿,
  채널 이름·순서·GUID 대소문자, index 표지 없음, frame interval, RGB 순서, 선 주사 판정, 참조 이미지만, ZIP64·암호화·이름, 진행 표시 다시 그리기, 채널 색·Info.
- 함께 찾은 결함: Remove Outliers 가 NaN 있는 32-bit 에서 ImageJ 와 80 화소 다름(§3g), 이미지 고르기 칸을 건드리지 않고 실행하면 보이는 첫 이미지가 아니라 '없음'(Image Calculator 포함).

## 6f. v1.5 (2026-10-17)
- 화면 골격 영어(문서 `docs/UI-English-style.md`): ImageJ 메뉴·대화 상자 이름, 패널·표·버튼. 알림·설명·경고는 한국어. 연산 id·매개변수 키·상태·저장 형식·.ijm 은 그대로.
- IBW 값 단위(§3d), OIR 실측 검증(§3b-3).

## 6g. v1.6 (2026-10-17) — Drop Shape (§3h) · 독립 검토 9 건 + 낮음 7 건
- 높음: SD 가 맞춤 통계뿐인데 전체 불확도처럼 보임 → 'fit SD' 이름, b·Bo 공분산, 바늘 눈금 SD 따로. 세실 'γ est.' 가 θ 작거나 Bo ≈ 0 에서 수백 배 틀림 → Bo > 0.1·상대 SD < 10 %·θ ≥ 60°·한계 미접촉일 때만. 영역이 방울을 자르면 조용히 틀리거나 빠짐 → '잘렸습니다' 경고, 방법이 빠진 이유를 경고에.
- 중간: Δρ ≤ 0 이면 음수 γ → 계산 안 하고 경고. 덧그림이 결과를 만든 장이 아닌 장에 → 만든 장 기록. 세실 시계열 그래프가 안 생김·실패 장이 있으면 x 밀림 → 이름 비교·실제 장 번호. 카드에 다른 측정 줄이 섞임 → 실행 id. 매크로를 다른 크기 이미지에 → 경고. 방법 차이 경고가 YL 쪽으로만 → 비대칭이면 YL 도 틀릴 수 있음을 함께.
- 낮음: 화소가 바뀐 뒤 옛 덧그림 숨김, CSV 에 Δρ·g·단위, 단위 micron·pm·μm, 비등방 화소 경고, NaN 화소 = 영역 중앙값(0 이 아니라), 큰 배열 Math.max(...) 제거. 남김: 좁은 영역에서 자동 기준선 실패(오류로 알림 — 선 ROI 사용), curveDist 탐색 창(수렴 뒤 영향 작음).

## 6h. Analyze Particles 보호 장치 · ROI 안만 · 단위 (2026-10-19)
- 사용자 보고: TEM 나노입자 이미지를 기본값으로 Analyze → 입자 수가 이상하고, 설정을 잘못하면 오래 걸리며 페이지가 멈춤.
  원인 ① 기본 문턱이 Dark background 켜짐(ImageJ 기본은 꺼짐) → TEM 명시야(밝은 바탕·어두운 입자)에서 바탕 70 % 를 입자로 → 14,608 개.
  ② 흐림 없이 문턱 → 잡음 점 150,595 개: 계산(2.8 s)보다 15 만 개 행·윤곽을 Worker 에서 옮기고 표로 그리고 저장하는 데서 멈춤.
- 고침: 기본 Dark background 꺼짐(ImageJ 와 같음, 저장된 상태는 그대로). 'Max count'(기본 2000, Analyze 패널) — 엔진이 조건에 맞는 입자가 이 수를 넘으면 측정을 멈추고 개수만 세어 돌려줌(1–2 s) → 원인 안내(선택 화소 > 50 % 면 극성, 아니면 Gaussian·Size 최솟값). 선택 화소 > 50 % 면 Threshold 패널에 경고 줄. 표는 앞 3000 줄만 그림(CSV·Copy 는 전부).
- 영역 ROI(사각형·타원·다각형·자유) 가 있으면 그 안만 셈(ImageJ 와 같음) — 축척 막대·글자를 뺄 때. Exclude on edges 는 사각형 ROI 테두리에 닿은 것도. 매크로에 ROI 기록.
- 화소 단위 고르기(IO.niceUnit): 화소 값이 0.1 이상이 되는 가장 큰 단위(mm·µm·nm), pm 은 0.01 nm 미만에서만. 예전엔 0.55 nm/px 가 552.6 pm 으로 보여 Size 칸에 nm² 로 생각하고 넣은 값이 10⁶ 배 틀렸다.
- 시험 test/ui_pa_guard.py 9: 합성 TEM 명시야(원판 421, 붙은 쌍 21, 잡음 σ 25, 축척 막대·글자) — 흐림 없이 0.7 s 에 상한·안내, 극성 경고, Gaussian σ2 → Mask → Fill Holes → Watershed → ROI 안 = 421 개(정답), 매크로 ROI, Max count 올리면 17,945 개 2.9 s·표 3000 줄.

## 6i. 원클릭 Analyze particles (auto) · '영상' → '이미지' (2026-10-19)
- 사용자 요청: 권장값을 기본으로, TEM 이미지를 넣고 단추 하나로 합리적 옵션 → 크기 분포. paAuto(30_analyze.js): Duplicate(원본 보존) → Gaussian σ 2 → Threshold Default(Default 문턱으로 어두운 쪽이 > 50 % 면 Dark background) → Fill Holes → Watershed → 축척 막대 자동 제외(DIFF_OCR.findBars, 가장자리 20 % 안 막대, 막대 ± 0.35·길이 상자 안에 중심이 있는 덩어리만 지움; 눈금 없으면 readAuto conf ≥ 0.5 로 눈금) → Size ≥ ¼·성분 넓이 중앙값(9 px 이상 성분), Circ 0.5–1, edge 제외, redirect = 원본 → Size distribution 탭. 각 단계는 보통 연산이라 매크로에 기록.
- 기본값(새 상태 + 저장된 상태 한 번 이전 _d1019): Circularity 0.5–1, Exclude on edges, Dark background 꺼짐. 손으로 하는 Analyze 는 '고급'으로 아래에.
- 사용자 TEM 이미지: Auto 5.7 s → 432 개, ECD 23.7 ± 4.6 nm, D50 23.3 nm, 축척 막대 둘레 13 덩어리 제외. 합성 421 개: ROI 있음 416 · 없음 421.
- 화면·안내 문구의 '영상' 을 '이미지' 로(조사 맞춤: 이→가, 은→는, 을→를, 으로→로 …, 괄호 뒤 조사 포함) — image·diffraction·touch 모듈과 사용 가이드.

## 7. 한계 (다음 단계)
- vsi, ND2 v1(JPEG 2000), CZI JPEG XR, JPEG 압축 TIFF 는 아직 읽지 않는다. OIR 보조 파일(_00001…)을 함께 읽지 않는다. 선 주사·RGB OIR 은 실측 미검증. AFM 힘 곡선(Force·VOLM·curve map)은 읽지 않는다.
- 3D 연산(3D skeletonize·필터), bilinear·비강체(elastic) 정합은 없다.
- 2.4×10⁷ 화소를 넘는 이미지는 전체를 볼 때 줄여 보이고(확대하면 원래 해상도), 문턱 겹쳐 그림은 줄인 해상도.
- 터치 두 손가락 확대 없음.
- Drop: 거울 같은 기판(반사상)은 자동 기준선이 틀릴 수 있어 선 ROI 로 직접 정한다. 기울어진 기판의 비축대칭 방울(경사판 미끄럼각)·진동 방울·액막(captive bubble)은 아직 없다. 전진·후퇴각은 바늘 꽂힌 방울의 스택 시계열로 잰다(따로 자동 판정 없음).
