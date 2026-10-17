# Diffraction · Phase DB — API (diff_db.js / diff_db_ui.js)

넣는 순서(같은 HTML 안 classic `<script>`): `diff_sg.js → diff_scat.js → diff_core.js → [diff_pack.js] → diff_db.js → diff_db_ui.js`, CSS `diff_db.css`.
`diff_pack.js` 는 없어도 된다. 나중에 `window.DIFF_PACK` 이 생기면 `list()/get()` 이 알아서 집어 든다(객체가 바뀌면 다시 읽음).

## 층 (layer)
| layer | 어디 | id | 쓰기 |
|---|---|---|---|
| `pack` (화면: Built-in) | `DIFF_PACK.phases` | `pack:…` | 읽기 전용 |
| `lab` | Firestore `inno-analytics-base` → IndexedDB 캐시 | `lab:…` (문서 id 와 같음) | 구성원, Firestore 규칙대로 |
| `local` | IndexedDB `inno-diffraction` (이 브라우저만) | `local:…` | 자유 |

IndexedDB `inno-diffraction` v1 (stores `local`, `lab`, `cal`, `meta`)은 diff_db.js 가 가진다. Diffraction 탭의 다른 저장은 **다른 DB 이름**을 쓸 것(같은 이름에 다른 store 를 만들려면 판 올림을 맞춰야 함).

## DIFF_DB (전역, Node 에선 `module.exports`)
반환 기록은 모두 사본(`layer` 붙음). 저장 결과는 `{ok:true, rec, validation}` 또는 `{ok:false, reason:"invalid"|"duplicate", validation, duplicate}`.

**시작·읽기**
- `init()` → Promise. IndexedDB 열기 + 캐시·메타 읽기 + pack. 다른 함수(비동기)는 스스로 `init()` 을 기다린다. `ready()` 같은 것.
- `list({layers?=["pack","lab","local"], includeHidden?=false, query?, elements?, elementMode?="subset"|"all"("contains" 도 같음)})` → 기록[] (이름순).
  - `elements`: `"Na Y F"` 또는 `["Na","Y","F"]`. `subset` = 상의 원소 ⊆ 고른 원소 (REDI 의 EDS 원소), `all` = 고른 원소를 모두 포함.
  - **조성비**(2026-10-17, REDI 의 '성분 + 조성비'): `elements` 에 수까지 — `"Fe3 O4"`·`"Fe3O4"`·`"Fe:3 O:4"`·`"Fe 3 O 4"`·`"Fe=42.1 O=57.9"`(EDS at.%). 모든 원소에 수가 있으면 넣은 원소끼리의 몫(합 1)과 상 화학식의 같은 원소 몫이 원소마다 `ratioTol`(기본 0.05 = 5 at.%p) 안이어야 통과. `useRatio:false` 면 비를 무시. 수가 일부만 있으면 비를 쓰지 않는다(`partial`).
  - `parseComposition(str|arr)` → `{els, ratio|null, partial, bad[]}`, `matchRatio(rec, ratio, tol)`. 화학식은 `DC.parseFormula` — 괄호·대괄호·배수(`Ca3(PO4)2`)와 쉼표 혼합 자리(`(Fe,Mn)2SiO4` → Fe 1·Mn 1)를 푼다.
  - `query`: 이름·id·COD id·DOI·공간군 기호·번호 부분 일치, 화학식은 공백 무시 부분 일치 + **원소 비 일치**(`TiO2` = `Ti2 O4`).
- `get(id)`, `layerOf(id)`, `counts()` → `{pack, lab, labHidden, local, cal}`, `packInfo()` → `{version, count}`, `elementsOf(rec)`, `parseElements(str|arr)`.
- `phase(id)` → `DC.phaseFromRecord` 결과 (판이 바뀔 때만 다시 만듦). 숨긴 Lab 항목도 된다(재현용).
- `validate(rec)` = `DC.validateRecord(contentOf(rec))`. `findDuplicate(rec, {exclude?:[id], layers?})` → 첫 중복 기록 또는 null (`DC.sameRecord`, 숨긴 것 제외).
- `contentOf(rec)` → PhaseRec 내용 필드만 정규화(`name, formula, cell, sg, symops?, atoms?, source, note?`). `contentHash(rec)` → 8 자리 FNV-1a.

**판정 재현 (결과에 반드시 저장)**
- `refOf(id)` → `{id, layer, version, updatedAt, hash, name, packVersion?(pack), projectId?(lab)}`.
  Diffraction 결과(상 판정·보정)에 쓴 상마다 이 객체를 그대로 저장한다. 예: `result.phaseRefs = ids.map(DIFF_DB.refOf)`.
- `getVersion(ref)` → Promise `{rec, exact, note?}` | null. 지금 판이 같으면 그대로, Local 은 기록 안의 history(최근 30 판), Lab 은 `phases/{id}/history/{version}` (구성원·온라인). `exact=false` 면 내용(hash)이 저장 때와 다르다 → 결과 화면에 경고.
- `verifyState(id)` (lab) → `{state:"unverified"|"verified"|"edited", verifiedVersion, version}` — 검증 뒤 주인이 고치면 규칙상 `verified` 가 그대로라 `verifiedVersion < version` 으로 "edited" 를 알린다.

**Local**
- `saveLocal(rec, {force?, checkDuplicate?})` — `rec.id` 가 기존 local 이면 고치기(version+1, 옛 내용을 `history` 로), 없으면 새로(중복이면 `reason:"duplicate"`, `force` 로 무시). 내용이 같으면 `{unchanged:true}`.
- `removeLocal(id)`, `copyToLocal(id)` (pack/lab → local 사본, `copiedFrom`).
- `importCIF(text, {fileName?})` → `[{block, file, rec?, validation?, duplicate?, error?}]` (저장은 안 함; 블록마다 따로, 깨진 블록은 `error`). 출처: COD id → `COD`, DOI → `DOI`, 학술지 정보 → `paper`(citation), 없으면 `manual`; `source.note = "CIF 파일 · data_블록"`.
- `exportJSON(ids?, {layers?=["lab","local"], calibrations?})` → 문자열 `{format:"inno-diffraction-phasedb", v:1, exportedAt, projectId, packVersion, phases[], calibrations?}` (Local 은 history 포함).
- `exportCSV(ids?, {layers?})` → 요약 CSV (`id,layer,name,formula,sg_number,…`).
- `importJSON(text|obj)` → `{added, updated, unchanged, skipped:[{id,name,reason}], errors:[…], calibrations}`. `local:` id 는 그대로 덮어쓰기, 다른 층 항목은 이미 있거나 중복이면 건너뛰고 아니면 Local 사본(`importedFrom`).

**Lab (Firestore)** — SDK 는 `connect()/signIn()/sync()` 때만 `import("https://www.gstatic.com/firebasejs/10.12.2/…")`.
- `connect()` → 상태. SDK 불러오기 + 앱(`inno-diffraction-<projectId>` 이름) + redirect 결과 + 이전 로그인 되살리기. 실패해도 throw 하지 않고 `status.error`.
- `firebaseStatus()` → `{projectId, expectedProjectId, configOk, configOverridden, sdkLoaded, connected, fileProtocol, online, user:{uid,email,emailVerified,name}|null, isMember(true|false|null), isAdmin, admins, lastSync, cacheProject, labCount, ignoredCache, wasSignedIn, storage:"indexeddb"|"memory", error:{code,message}|null}`.
  `error.code`: `file-protocol`, `sdk-load`, `not-member`, `email-unverified`, `unavailable`, `permission-denied` … (메시지는 한국어).
- `signIn()` — 팝업, `auth/popup-blocked`·`operation-not-supported-in-this-environment`·`internal-error`·`web-storage-unsupported` 이면 redirect (`{redirect:true}` 반환, 돌아오면 `connect()` 가 마무리). file:// 이면 `file-protocol` 오류. 구성원이면 바로 `sync()`.
- `signOut({clearCache?})` — 캐시는 기본으로 남김(오프라인 사용).
- `sync()` → `{phases, hidden, calibrations, lastSync}`. `getDocsFromServer` 만 씀(오프라인 캐시가 빈 목록을 돌려 캐시를 지우는 일 방지). 실패하면 캐시 그대로 + 오류.
- `addLab(rec, {force?, exclude?})` — 쓰기 권한 확인 → 검사 → 모든 층 중복 확인(중복이면 쓰지 않고 돌려줌) → `version 1, hidden false, verified false, createdBy uid` (관리자가 `rec.verified===true` 면 검증 상태로).
- `updateLab(id, rec, {checkDuplicate?})` — 트랜잭션: 옛 문서 → `history/{옛 version}` (create), 새 내용 `version+1`, `createdBy`·`verified`(비관리자) 유지. 주인·관리자만.
- `hideLab(id, hidden=true)` (지우기 대신), `setVerified(id, bool)` (관리자; `verifiedVersion` 기록), `promoteLocal(id, {force?})` (local → lab, local 에 `promotedTo`), `labHistory(id)` → `[{version, updatedAt, updatedByEmail, archivedAt, archivedByEmail, archivedReason, rec}]` (최신 먼저).
- 보정: `listCalibrations({layers?})`, `saveCalibration(cal, {lab?:bool})` — 필드 `name, microscope, kV, cameraLength(표시값 mm), mode, scale(1/nm/px), ellip, phi, sigmaRel, standard(상 id), standardName, sourceFile, note` (그 밖의 필드·넘긴 `id` 는 무시, 새 id `local:…`/`lab:…`); `standard` 가 있으면 `standardRef = refOf(standard)` 를 함께 저장. Lab 은 만든 사람·관리자만 고침.
- 구성원(관리자): `listMembers()`, `addMember(email, name)` (소문자로), `setMemberActive(email, bool)`, `listAdmins()`, `addAdmin(email)`, `removeAdmin(email)` (마지막 한 명은 `last-admin` 오류).
- `config()`, `setConfig(cfg|null)` — 다른 projectId 면 `configOk:false`, 화면에 위험 경고, **모든 쓰기 거부**(`wrong-project`), 다른 프로젝트에서 받은 Lab 캐시는 숨김.

**이벤트**: `onChange(cb)` → 해제 함수. `cb({type})`, type = `init | pack | local | lab | cal | status | members` (+ `ids`, `sync`, `import`).

**시험용**: `_setIDB(factory)`, `_setSDK({app,auth,firestore})`, `_setEnv({protocol, online})`, `_reset()`.

## DIFF_DB_UI
- `render(container, {onUse?(ids, refs), autoConnect?=true})` — Phase DB 화면 전체(상태 줄·도구 막대·목록·상세·보정·구성원). `onUse` 를 주면 목록 아래 "Use selected (N)" 단추가 생기고 고른 id 와 `refOf` 를 넘긴다(판정 후보로).
  `autoConnect`: 전에 로그인했던 브라우저면 화면을 열 때 SDK 를 불러 로그인을 되살리고, 마지막 동기화가 10 분 넘었으면 조용히 sync.
- `openEditor(rec?)` — 추가(인자 없음)·편집 대화상자. `rerender()`, `_state()`.
- 화면 선호(켜 둔 층·원소)는 `localStorage["inno-diff-db-ui"]` (실패해도 무시).

## 시험
- `node test/test_db.js` (fake-indexeddb, firebase@10.12.2 는 `npm i --no-save fake-indexeddb firebase@10.12.2` — node_modules 는 .gitignore)
- `python3 test/ui_db.py <그림 폴더> [앱 빌드 html]` — Playwright, 네트워크 막음, Firebase 는 `test/fake_firebase.js`.
- `test/fake_firebase.js` 는 게시된 Firestore 규칙을 그대로 집행하는 메모리 가짜(오류 = 거부). 실제 Firestore 에뮬레이터 시험은 아직 없음.
