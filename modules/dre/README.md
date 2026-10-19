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

## 10-19 수정
- 화면 시뮬레이션: 방향마다 스윕 한 번으로 모든 채널을 읽는다(`DRE_CORE.observeMulti`). 계산 상한(점 60 ms, 스윕 1.5 s, 실패 8회면 중단 — `sweep(…, {maxFails, deadline, stOpt:{maxMs}})`)으로 오타 값에도 탭이 멈추지 않는다. 피팅은 상한 없이 예전과 같다. 수렴 실패는 S-curve 제목 옆 배지("Steady state did not converge at k/N points"), 탐색 범위에서 10³배 넘게 벗어난 값도 배지.
- λ 200–3000 nm 로 제한, 잘못된 λ·Range 입력은 알림 후 원래 값으로 되돌림. 템플릿 파장과 다르면 "Cross-sections are for X nm; only photon flux rescaled" 표시.
- ED 가져오기: 농도가 없는 종은 "c defaulted to 1 %" 표시.
- 그래프는 그릴 자리 폭으로 그리고(최소 320 px) 창 크기가 바뀌면 계산 없이 다시 그림. 매개변수 표: 탐색 범위를 값 아래 줄로, 오른쪽 열 minmax(450px,35%).
- 시험: `node test/test_simbudget.js`, UI `python3 test/ui_fix.py <빌드 html>` (ED·DRE 수정 항목 전체).
