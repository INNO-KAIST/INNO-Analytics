# XPS 모듈 (INNO Analytics)

XPS 읽기·배경·피팅·불확도·정량·survey 판별 탭. **배포본(HTML)을 직접 고치지 말고 여기를 고친 뒤 `../inject_all.py` 로 주입한다.**
수식의 정의는 `SPEC.md` 하나 — 코드와 다르면 SPEC 을 고치고 코드를 맞춘다.

| 경로 | 내용 |
|---|---|
| `SPEC.md` | 수식 명세(선모양·배경·피팅·불확도·정량·대전 보정·잡음·survey 판별·라이브러리·자동 모형) |
| `src/xps_io.js` | 파일 읽기(VGD·VAMAS·Avantage xls/xlsx/avg·텍스트). DOM 없음 |
| `src/xps_core.js` | 계산 엔진. DOM 없음 — 주입 시 `<script id="xps-core-src">`, Worker 는 이 태그의 글을 Blob 으로 띄운다 |
| `src/xps_ui.js`, `src/xps.css` | 화면 |
| `src/xps_lib.js` | 원소 라이브러리 — **생성물**(`build/build_lib.py`). 직접 고치지 말 것 |
| `build/build_lib.py` | `data/` → `src/xps_lib.js`. 검사(Scofield 논문 대조·doublet 합·Auger 대조·순서 역전)가 어긋나면 중단 |
| `build/assemble_core.py`, `build/extract_kherve.py` | 원자료 정리(X-ray Data Booklet 전사 · KherveFitting 수치 추출) |
| `data/` | 원자료와 대조용 판독값(출처는 각 JSON 의 `source`) |
| `test/test_xps.js` | 회귀·검증 145 항목 |
| `test/all_elements.js` | 원소 독립성 — Li–U 143 원소·궤도 합성 시험(test_xps §7 이 부름, 단독 실행 가능) |
| `test/ref/ref_xps.py` → `ref.json` | 기준 구현 ① (SPEC 만 보고 작성, scipy) |
| `test/ref/xps_ref.py`, `make_fixtures.py`, `validate.py` → `fixtures/` | 기준 구현 ② (SPEC 만 보고 따로 작성, 닫힌 꼴·제3자 패키지 대조). 자세한 것은 `test/ref/README.md` |
| `test/survey_bench.js` | 합성 survey 판별 성능 측정 |
| `test/ui_xps.py` | 브라우저(Chromium) 실측 |
| `test/data/` | 실측 회귀 자료(BCT xls·txt·VAMAS, K2312 VGD) — 반례 감시용, 모형을 맞추는 데 쓰지 않음 |

```bash
python3 -I build/build_lib.py            # 라이브러리 다시 만들기 (data/ 를 바꿨을 때)
python3 -I test/ref/ref_xps.py           # 기준 ① 다시 만들기 (SPEC 을 바꿨을 때, 약 1 분)
(cd test/ref && python3 make_fixtures.py all && python3 validate.py)   # 기준 ② (약 5 분)
node test/test_xps.js                    # 145 항목 (약 2 분)
python3 ../inject_all.py "../../INNO Analytics.html" "../../INNO Analytics.html"
python3 test/ui_xps.py "../../INNO Analytics.html" <그림 폴더>
```
