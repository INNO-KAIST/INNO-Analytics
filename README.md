# INNO Analytics

STEM 이미지·EDS 스펙트럼·SI(원소 맵)·격자 FFT·EDS 간섭 검사·Energy Diagram·DRE(속도식 모델)·XPS·Image Analysis(ImageJ 방식) 분석을 브라우저 안에서 하는 단일 HTML 앱. 서버·로그인 없음.

## 폴더 지도
```
INNO Analytics/
├─ INNO Analytics.html      ← 앱 본체(이것 하나로 동작). GitHub Pages 배포 대상
├─ README.md                ← 이 폴더 안내
├─ data/
│   └─ kfactors_NaLnF4.json ← EDS 정량용 k-인자 세트(앱의 ‘k세트 불러오기’로 읽음, 흡수 보정 없이 잰 값)
├─ docs/
│   ├─ 사용 가이드.md
│   └─ HISTORY.md           ← 작업 이력·교훈·문서 색인
├─ modules/                 ← 탭 모듈 원본 (HTML 에 주입됨)
│   ├─ inject_all.py        ← python modules/inject_all.py "INNO Analytics.html" "INNO Analytics.html"
│   ├─ energy_diagram/  dre/  interf/  xps/  image/   (src·data·build·test)
│   ├─ dre_study/           ← DRE 연구 기록(스크립트·결과 tgz·그림)
│   └─ terms_map.py
└─ _drive_backup/           ← Drive 보관용 (RESTORE.md · 백업 만들기.bat)
```

## 쓰는 법 · 고치는 법 · 배포
- 본체(`INNO Analytics.html`)는 **HTML 자체가 원본**이다. 다만 Energy Diagram·DRE·간섭 검사 추천·XPS·Image Analysis 탭은
  `modules/*/src` 가 원본이고 `modules/inject_all.py` 가 HTML 에 표지(`<!-- INNO-ED:BEGIN -->` 등) 사이로 주입한다(멱등).
  모듈을 고쳤다면: `python modules/inject_all.py "INNO Analytics.html" "INNO Analytics.html"`.
- 배포: `python "../INNO Deploy/inno_deploy.py" push analytics` → https://inno-kaist.github.io/INNO-Analytics/

- 사용법: `docs/` 의 가이드 · 작업 이력과 교훈: `docs/HISTORY.md`
- 백업(Google Drive): `_drive_backup/백업 만들기.bat` 실행 → `_drive_backup` 폴더를 Drive 에 올림. 복구는 `_drive_backup/RESTORE.md`.
