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
| `src/xps_lib.js` | 원소 라이브러리 — **생성물**(`build/build_lib.py`). 직접 고치지 말 것 |
| `build/build_lib.py` | `data/` → `src/xps_lib.js`. 검사(Scofield 논문 대조·doublet 합·Auger 대조·순서 역전)가 어긋나면 중단 |
| `build/assemble_core.py`, `build/extract_kherve.py` | 원자료 정리(X-ray Data Booklet 전사 · KherveFitting 수치 추출) |
| `data/` | 원자료와 대조용 판독값(출처는 각 JSON 의 `source`) |
| `test/test_xps.js` | 회귀·검증 (§10 = 10-19 결함 회귀) |
| `test/all_elements.js` | 원소 독립성 — Li–U 143 원소·궤도 합성 시험(test_xps §7 이 부름, 단독 실행 가능) |
| `test/ref/ref_xps.py` → `ref.json` | 기준 구현 ① (SPEC 만 보고 작성, scipy) |
| `test/ref/xps_ref.py`, `make_fixtures.py`, `validate.py` → `fixtures/` | 기준 구현 ② (SPEC 만 보고 따로 작성, 닫힌 꼴·제3자 패키지 대조). 자세한 것은 `test/ref/README.md` |
| `test/survey_bench.js` | 합성 survey 판별 성능 측정 |
| `test/test_xps_casa.js` | v1.7 CasaXPS/Avantage 대응 기능 시험: 측정 선모양·LCF·기준 잡음 MC(§12), spline 배경(§2), Simplex(§14), 여러 스펙트럼 일괄(§15), at% MC(§4), CasaXPS 2.3.27 실측 대조(§2b·§13), CasaXPS·Avantage 규약 기준 대조(§2b·§2c·§1.10·§4). 약 5 분, `QUICK=1` 이면 약 40 초 |
| `test/ref/ref_casa.py` → `ref_casa.json` | 위 기능의 독립 기준 구현(PchipInterpolator·quad·nnls·least_squares·CubicSpline·Nelder–Mead·numpy MC) |
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
python3 ../inject_all.py "../../INNO Analytics.html" "../../INNO Analytics.html"
python3 test/ui_xps.py "../../INNO Analytics.html" <그림 폴더>
python3 test/ui_xps_fixes.py "../../INNO Analytics.html" <그림 폴더>
python3 test/ui_xps_casa.py "../../INNO Analytics.html" <그림 폴더>
python3 test/ui_xps_multi.py "../../INNO Analytics.html" <그림 폴더>   # 약 3 분
QUICK=1 node test/test_multi.js; QUICK=1 node test/test_arxps.js; node test/test_image.js; QUICK=1 node test/test_proc.js
```
