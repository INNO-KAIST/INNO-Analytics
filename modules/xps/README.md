# XPS 모듈 (INNO Analytics)

XPS 읽기·배경·피팅·불확도·정량·survey 판별 탭. **배포본(HTML)을 직접 고치지 말고 여기를 고친 뒤 `../inject_all.py` 로 주입한다.**
수식의 정의는 `SPEC.md` 하나 — 코드와 다르면 SPEC 을 고치고 코드를 맞춘다.

| 경로 | 내용 |
|---|---|
| `SPEC.md` | 수식 명세(선모양·배경·피팅·불확도·정량·대전 보정·잡음·survey 판별·라이브러리·자동 모형) |
| `src/xps_io.js` | 파일 읽기(VGD·VAMAS·Avantage xls/xlsx/avg·텍스트). DOM 없음 |
| `src/xps_core.js` | 계산 엔진. DOM 없음 — 주입 시 `<script id="xps-core-src">`, Worker 는 이 태그의 글을 Blob 으로 띄운다 |
| `src/xps_ui.js`, `src/xps.css` | 화면. v1.8: 위쪽 작업 공간 전환 **Spectra · Series · ARXPS · Imaging** + Peak Fit 옆 **Process** 탭 (SPEC.md §16) |
| `SPEC_multi.md`, `src/xps_multi.js` · `xps_arxps.js` · `xps_image.js` | series·다변량(PCA·MLPCA·TFA·MCR-ALS·기준 NNLS) · ARXPS · 분광 영상 엔진. DOM 없음 |
| `SPEC_proc.md`, `src/xps_proc.js` | 디컨볼루션 · 다시 표본화 · 평활/미분 · Fermi edge · 기준 봉우리 PSF 엔진. DOM 없음 |
| `src/xps_guide.js` | v1.9.3 피팅 안내(SPEC.md §17) — Fit check(피팅 뒤 진단 R1–R9 · 한 번 누르기 조치), 영역 레시피(Ce 3d Burroughs v/u · v1.9.6 Ni·Co·Fe·Mn·Cr·Cu 2p3/2 Biesinger 포락선 — 화학종 체크, PE 열, FWHM×, Fe₃O₄ 1:2, 화학종 분율 ± σ, Cu(II) shake-up 식 1–2; §17.7), 단계 피팅(stagedFit · 보통 피팅보다 나빠지지 않음). DOM 없음, core 다음에 주입(`xps-guide-src`, Worker 에도) |
| `src/xps_lib.js` | 원소 라이브러리 — **생성물**(`build/build_lib.py`). 직접 고치지 말 것 |
| `build/build_lib.py` | `data/` → `src/xps_lib.js`. 검사(Scofield 논문 대조·doublet 합·Auger 대조·순서 역전)가 어긋나면 중단 |
| `build/assemble_core.py`, `build/extract_kherve.py` | 원자료 정리(X-ray Data Booklet 전사 · KherveFitting 수치 추출) |
| `data/` | 원자료와 대조용 판독값(출처는 각 JSON 의 `source`) |
| `test/test_xps.js` | 회귀·검증 (§10 = 10-19 결함 회귀) |
| `test/all_elements.js` | 원소 독립성 — Li–U 143 원소·궤도 합성 시험(test_xps §7 이 부름, 단독 실행 가능) |
| `test/ref/ref_xps.py` → `ref.json` | 기준 구현 ① (SPEC 만 보고 작성, scipy) |
| `test/ref/xps_ref.py`, `make_fixtures.py`, `validate.py` → `fixtures/` | 기준 구현 ② (SPEC 만 보고 따로 작성, 닫힌 꼴·제3자 패키지 대조). 자세한 것은 `test/ref/README.md` |
| `test/survey_bench.js` | 합성 survey 판별 성능 측정 |
| `test/fitlog_hook.js`, `test/fitlog_compare.js` | LM 멈춤 규칙(SPEC §3.3b) 검증: `node -r ./test/fitlog_hook.js test/<시험>.js` 를 `XPS_FITLOG=<파일>`(·`XPS_NOSTALL=1` 옛 동작)로 돌려 모든 피팅을 기록 → 두 폴더 비교 |
| `test/stall_gen.js` | 멈춤 규칙 시험용 무작위 합성 스펙트럼 생성기(seed → 문제; 성분 1–6·모든 모양·ST/SB·배경 5 종·가중 3 종·약 25 % 는 모양 매개변수까지 자유) |
| `test/stall_neverworse.js` | 멈춤 규칙 '옛 LM 보다 나빠지지 않음' 회귀 시험: 같은 문제를 한 프로세스에서 규칙 끔·켬으로 번갈아 피팅 → 최종 χ² 새 ≤ 옛·(1+1e-7), 갈래 없는 피팅 비트 동일, 시간 합 (기본 seed 1–1100; 고정 표본은 `test_shapes2.js` §4c) |
| `test/stall_stress.js`, `test/stall_stress_cmp.js` | 같은 무작위 합성 스펙트럼을 한 엔진으로 피팅해 기록(옛 = `XPS_NOSTALL=1`) → 두 기록의 χ² 비·시간 비교 (SPEC §3.3b 표) |
| `test/test_xps_casa.js` | v1.7 CasaXPS/Avantage 대응 기능 시험: 측정 선모양·LCF·기준 잡음 MC(§12), spline 배경(§2), Simplex(§14), 여러 스펙트럼 일괄(§15), at% MC(§4), CasaXPS 2.3.27 실측 대조(§2b·§13), CasaXPS·Avantage 규약 기준 대조(§2b·§2c·§1.10·§4). 약 5 분, `QUICK=1` 이면 약 40 초 |
| `test/test_casa_v3.js` | CasaXPS 2.3.27 Round 3–6 실측 대조(SPEC_casa_v3·v4): 폭 숫자·asym 규칙, DS 넓이·위치, ST 모양·위치, 표기 왕복, casaModelToINNO. 약 2 초. `test/ref/check_casa_v4.js [core]` 는 같은 판독의 표 출력 |
| `test/test_guide.js` | v1.9.3 피팅 안내(SPEC §17): 레시피 자료·대전 보정 고정·정렬, 2차 차분 잡음(MA(1)), 규칙마다 양성 대조, 잘 놓인 피팅 거짓 경보(270 개), 조치, 단계 피팅 '나빠지지 않음', 합성 Ce 3d(Ce³⁺ 0·20·50 %) 회복 · 순진한 Auto suggest 와 한계 닿음 비교 · §8 포락선 레시피: 숫자 436 개를 `test/ref/biesinger/` 원문 전사와 대조, 제약 구조, 합성 2–3 화학종 혼합 18 개 회복, Cu shake-up, Fit check. 약 2–3 분, `QUICK=1` 약 30 초 |
| `test/test_real_ce.js` | 실측 사례 S3(CeO₂, 비공개 `test/data/real/ce`, 없으면 건너뜀): 순진한 피팅의 Fit check · 레시피 + 단계 피팅의 Ce³⁺·u‴ 분율·χ²ν·한계 닿음 · 대전 기준 어긋남 |
| `test/ui_xps_guide.py` | 브라우저 실측: 위성 영역 한 줄 안내 + Recipe 선택 · 순진한 피팅 → Fit check(배지·title·Apply recipe) → 레시피 → Staged fit(Worker) → Ce³⁺·u‴ 표시 · 빠진 선 → Add a line → 다시 피팅 · Ni 2p 포락선(화학종 체크 → Recipe → Staged fit → 분율) · Cu/Fe 함정 title · 390 px |
| `test/test_real_a2.js` | K-Alpha 실측 a2(비공개 `test/data/real/a2`, 없으면 건너뜀): Avantage 규약 v1.9.2 Smart·Shirley 배경 · Area(P) · at% · VGD prop 415 (TF 제거) — SPEC §2c·§4. src 에 avantage_a2 패치가 들어 있는지 먼저 확인 |
| `test/ref/ref_casa.py` → `ref_casa.json` | 위 기능의 독립 기준 구현(PchipInterpolator·quad·nnls·least_squares·CubicSpline·Nelder–Mead·numpy MC) |
| `test/ref/biesinger/*.txt` | Biesinger 2011 표 1·3–6 · 2017 표 4(·A1s/Bs·α′) 원문 전사 — v1.9.6 포락선 레시피의 숫자 출처, test_guide §8 이 대조 |
| `test/ref/casa_readings.json`, `test/data/casaval/` | CasaXPS 2.3.27 화면 실측값과 그 자료(합성 + BCT C1s) — CasaXPS 코드·자료 아님 |
| `test/ref/make_casaval2.py` | CasaXPS 추가 판독용 무잡음 합성 VAMAS(`test/data/casaval/V2_*.vms`: LA 폭 환산·f_G(n)·SGL 가설) |
| `test/ui_xps_casa.py` | 브라우저 실측: 기준 저장 → LCF → MC → 일괄 피팅·CSV → 측정 선모양 → spline 매듭(클릭·끌기·두 번 누르기)·민감도 → Simplex → at% MC → CasaXPS 규약 → 390 px |
| `test/ui_xps.py` | 브라우저(Chromium) 실측 |
| `test/ui_xps_multi.py` | v1.8 화면 실측: SDP·각도·MAP VAMAS·여러 열 텍스트(`test/data/synthetic/`) → Series(at % ±3 at.% · PCA · MCR) · ARXPS(d ≈ 1.8 nm · 다층 · RDP · MaxEnt) · Imaging(넓이 z · ROI → 피팅 → NNLS 상태 지도 · PCA · binning · 텍스트 stack) · Process(디컨볼루션 · Fermi edge) · 2 s 넘는 주 스레드 작업 없음 · 390 px · 다시 열기 |
| `test/ui_xps_fixes.py` | 10-19 감사 결함 회귀(템플릿 왕복·C 1s 기준 선택·피팅 중 수정·Auger 짝·정량 기본 포함·KE 축·쉼표 소수·배경 뺀 자료·390 px) — 항목마다 판정 |
| `test/data/` | 실측 회귀 자료(BCT xls·txt·VAMAS, K2312 VGD) — 반례 감시용, 모형을 맞추는 데 쓰지 않음 |

```bash
python3 -I build/build_lib.py            # 라이브러리 다시 만들기 (data/ 를 바꿨을 때)
python3 -I test/ref/ref_xps.py           # 기준 ① 다시 만들기 (SPEC 을 바꿨을 때, 약 1 분)
(cd test/ref && python3 make_fixtures.py all && python3 validate.py)   # 기준 ② (약 5 분)
python3 -I test/ref/ref_casa.py         # 기준 (v1.7 기능, 약 25 초)
node test/test_xps.js                    # 약 6 분
node test/test_xps_casa.js               # 약 4 분 (QUICK=1 약 40 초)
node test/test_guide.js                  # 약 2–3 분 (QUICK=1 약 30 초) — v1.9.3 피팅 안내
node test/test_real_ce.js                # 실측 S3 (자료 없으면 건너뜀)
python3 ../inject_all.py "../../INNO Analytics.html" "../../INNO Analytics.html"
python3 test/ui_xps.py "../../INNO Analytics.html" <그림 폴더>
python3 test/ui_xps_fixes.py "../../INNO Analytics.html" <그림 폴더>
python3 test/ui_xps_casa.py "../../INNO Analytics.html" <그림 폴더>
python3 test/ui_xps_multi.py "../../INNO Analytics.html" <그림 폴더>   # 약 3 분
python3 test/ui_xps_guide.py "../../INNO Analytics.html" <그림 폴더>
QUICK=1 node test/test_multi.js; QUICK=1 node test/test_arxps.js; node test/test_image.js; QUICK=1 node test/test_proc.js
```
