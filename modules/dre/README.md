# DRE 모듈 (INNO Analytics)

rate equation 시뮬레이션·피팅 탭. **배포본을 직접 고치지 말고 여기를 고친 뒤 `../inject_all.py` 로 주입한다.**

| 경로 | 내용 |
|---|---|
| `src/dre_core.js` | 엔진 (DOM 없음, 메인·Worker·node 공용). 주입 시 `<script id="dre-core-src">` 로 들어가고, Worker 는 이 태그의 글을 Blob 으로 띄운다 |
| `src/dre_ui.js`, `src/dre.css` | 화면 |
| `test/test_dre.js` | 회귀 시험 56 항목 (`node test/test_dre.js
node test/test_adequacy.js   # 적합성 A1–A6 (fixtures_v4: WJ-0128c 실측 포함)
node test/test_insight.js     # 상세 분석 재료(공분산·해 집단·유량 분해·χ² 단면)`) |
| `test/lee_ref.json`, `test/v3_ref.json` | scipy 독립 구현(`make_ref.py`, `v3_bi.py`, `pa3.py`)으로 만든 기준값 |
| `test/v2_*.js` | 합성 데이터 식별성 시험. 결과는 `v2_*.json` |

```bash
node test/test_dre.js
python ../inject_all.py "../../INNO Analytics.html" "../../INNO Analytics.html"
```
suite `build_all.sh` 로 Analytics 를 원본에서 다시 만들었다면 마지막 단계로 `inject_all.py` 를 돌린다.
