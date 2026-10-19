# Energy Diagram 모듈 (INNO Analytics)

INNO Analytics 단일 HTML 에 'Energy Diagram' 탭을 주입하는 원본이다. **배포본(`INNO Analytics.html`)을 직접 고치지 말고 여기를 고친 뒤 다시 주입한다.**

## 구성

| 경로 | 내용 |
|---|---|
| `src/ed_core.js` | 물리·데이터 로직 (화학식 해석, 준위 생성, 다중포논 판정, ET/광흡수 판정, PA 루프 점검, 경로 제안). UI 가 없어 node 로 시험한다 |
| `src/ed_ui.js` | 화면 (SVG 도식, 편집 모드, 툴팁, 옆 패널, 내보내기) |
| `src/ed.css` | 스타일. 앱 토큰(`--accent` 등)을 따른다 |
| `src/ed_data.js` | **자동 생성** — `build/build_data.py` 가 `data/*.json` 을 압축해 만든다. 직접 수정 금지 |
| `data/` | 준위 원자료와 출처. `README_levels.md` 에 출처·방법·빈칸이 정리되어 있다 |
| `../inject_all.py` | 주입기 (Energy Diagram·DRE 공용, 멱등). `build/inject.py` 는 이것을 부른다 |
| `build/gen_ts.py`, `dn_ligand_field.py` | Tanabe–Sugano 격자 생성 (정확한 dⁿ 대각화) |
| `test/` | `test_core.js`(1,404 항목), `test_formula_warn.js`(농도 점검), `sweep.js`(27,266 점검), Playwright UI 시험 |

## 다시 만들기

```bash
python build/build_data.py                                   # data → src/ed_data.js
python ../inject_all.py "../../INNO Analytics.html" "../../INNO Analytics.html"   # 제자리 주입(멱등)
node test/test_core.js && node test/test_formula_warn.js && node test/sweep.js
```

UI 시험(`test/ui*.py`)은 파일 경로 `file:///home/claude/ed/out.html` 을 가리키므로, 쓰기 전에 그 줄을 실제 경로로 바꾼다.

**suite 파이프라인(`build_all.sh`: patch_* → … → ux_*)이 Analytics 를 원본에서 다시 만들면 이 탭이 빠진다.** 그 경우 마지막 단계로 `../inject_all.py` 를 한 번 더 돌린다.

## 본 앱 변경 (inject.py)

1. `VIEWS` 에 `["energy","Energy Diagram"]` 추가
2. 파일이 없어도 energy 탭은 환영 화면 대신 렌더
3. 탭 분기에 `window.ED_UI.render(root)` 추가
4. energy 탭에서 파일 출처 줄을 숨기고, `#app.ed-on` 으로 파일 열을 숨김
5. 휴대폰 하단 탭 격자 6칸 → 7칸 (아이콘 CSS 는 주입 블록에 있음)

## 10-19 수정
- 화학식 농도 점검: >100 %, ≤0 %, 같은 자리 도펀트 합 >100 % 경고. % 없는 숫자는 식 하나 안에서 한 가지로 읽는다(모두 ≤1 → 분율, 아니면 모두 %; 혼용이면 경고).
- `ED_UI._ensureModel()` — DOM 없이 저장 상태로 모델을 만든다(DRE 'Import from Energy Diagram' 이 탭을 열지 않아도 동작).
- Inter-ion 카드: 이온이 셋 이상이면 'Ion pair' 선택으로 한 쌍(양방향)만 보인다.
- 내보내기: 빈 도식은 "Nothing to export" 알림. CSV 에 `laser_cm1`·`transition_cm1`·`gap_cm1` 열 추가(donor_/acceptor_E 는 ET 행에만).
- 터치 기기: 준위 선택 영역을 이웃과 겹치지 않는 만큼(최대 24 px) 키움, 편집 안내의 'Esc = restart' 대신 Restart 버튼.
