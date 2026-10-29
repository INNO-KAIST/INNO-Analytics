"""INNO Analytics 단일 HTML 에 Energy Diagram·DRE·XPS·Image Analysis 탭을 주입한다 (멱등).
사용: python inject_all.py <입력 html> <출력 html>   (입력=출력 가능)
- 이전에 주입한 블록(INNO-ED / INNO-DRE / INNO-IF / INNO-XPS / INNO-IMG / INNO-TOUCH 표지)을 지우고 다시 넣는다.
- INNO-TOUCH = 공용 터치 제스처 window.INNO_TOUCH (touch/src) — head 에 CSS·JS 로.
- XPS 원소 라이브러리 src/xps_lib.js 는 modules/xps/build/build_lib.py 가 만든다(없으면 중단).
- 본 앱 코드 변경은 정규식 치환 + 개수 확인. 원본(주입 전)과 이미 주입된 판본 모두에서 같은 결과가 나온다."""
import sys, re, os
HERE=os.path.dirname(os.path.abspath(__file__))
ED=os.path.join(HERE,'energy_diagram','src'); DF=os.path.join(HERE,'diffraction','src'); TO=os.path.join(HERE,'touch','src'); DR=os.path.join(HERE,'dre','src'); IF=os.path.join(HERE,'interf','src'); XP=os.path.join(HERE,'xps','src'); IM=os.path.join(HERE,'image','src')
inp,outp=sys.argv[1],sys.argv[2]
s=open(inp,encoding='utf-8').read()
s=re.sub(r'<!-- INNO-(ED|DRE|IF|XPS|IMG|TOUCH|DIFF|NAV|TIPS):BEGIN[^>]*-->.*?<!-- INNO-\1:END -->\n?','',s,flags=re.S)
def R(desc,pat,rep,flags=0):
    global s
    n=len(re.findall(pat,s,flags))
    assert n==1, f'[{desc}] 일치 {n}개 (1개여야 함)'
    s=re.sub(pat,rep,s,flags=flags)
# XPS 는 맨 왼쪽(사용자 요청 10-15), Energy Diagram·DRE 는 메타데이터 뒤
# 탭 이름은 영어(docs/UI-English-style.md). 옛 한글 이름·옛 영어 이름·이미 주입된 판본 모두 같은 결과가 나오도록 라벨은 [^"]* 로 받는다.
# 키(xps,image,imgx,maps,spectrum,lattice,interf,meta,energy,dre)는 바꾸지 않는다.
# EM Metadata 는 탭에서 빼고 EM Image 안의 'Metadata' 보기로 (10-19 사용자 요청 — 탭 수 줄이기). 키 meta 는 보기 이름으로 남는다.
TAB_LIST='const VIEWS=[["xps","XPS"],["image","EM Image"],["imgx","Image Analysis"],["maps","EDS Maps"],["spectrum","EDS Spectrum"],["lattice","Diffraction"],["interf","EDS Interference"],["energy","Energy Diagram"],["dre","DRE Kinetics"]];'
_L=r'"[^"]*"'
R('탭 목록',r'const VIEWS=\[(?:\["xps",'+_L+r'\],)?\["image",'+_L+r'\],(?:\["imgx",'+_L+r'\],)?\["maps",'+_L+r'\],\["spectrum",'+_L+r'\],\["lattice",'+_L+r'\],\["interf",'+_L+r'\],(?:\["meta",'+_L+r'\],?)?(?:,?\["energy",'+_L+r'\])?(?:,\["dre",'+_L+r'\])?(?:,\["xps",'+_L+r'\])?\];',
  lambda m:TAB_LIST)
# 휴대폰 하단 탭(열 10개, 폭 ≈39px): 짧은 이름. 라벨 글자는 DOM 에 남기고(접근성·테스트) CSS 가 data-short 를 보여 준다.
NAV_SHORT_JS='const NAV_SHORT={xps:"XPS",image:"EM Img",imgx:"Analysis",maps:"Maps",spectrum:"Spectra",lattice:"Diffr.",interf:"Interf.",energy:"Energy",dre:"DRE"};'
if 'const NAV_SHORT=' not in s:
    s=s.replace(TAB_LIST, TAB_LIST+'\n/* 휴대폰 하단 탭용 짧은 이름 — 키가 없으면 전체 이름 */\n'+NAV_SHORT_JS,1)
else:
    s,_n=re.subn(r'const NAV_SHORT=\{[^}]*\};',lambda m:NAV_SHORT_JS,s); assert _n==1, 'NAV_SHORT 일치 %d개'%_n
R('탭 버튼 라벨(짧은 이름)',r'return `<button data-v="\$\{k\}" class="\$\{(?:state\.view===k|\(state\.view===k\|\|\(k==="image"&&state\.view==="meta"\)\))\?"active":""\}"(?: aria-label="\$\{esc\(t\)\}")? \$\{dis\}>(?:<span class="nl" data-short="\$\{esc\(NAV_SHORT\[k\]\|\|t\)\}">\$\{t\}</span>|\$\{t\})\$\{cnt\}</button>`;',
  lambda m:'return `<button data-v="${k}" class="${(state.view===k||(k==="image"&&state.view==="meta"))?"active":""}" aria-label="${esc(t)}" ${dis}><span class="nl" data-short="${esc(NAV_SHORT[k]||t)}">${t}</span>${cnt}</button>`;')
R('파일 없이 열리는 탭',r'if\(!rec && state\.view!=="interf"(?: && state\.view!=="energy")?(?: && state\.view!=="dre")?(?: && state\.view!=="xps")?(?: && state\.view!=="imgx")?(?: && state\.view!=="lattice")?\)\{ renderWelcome\(root\); return; \}',
  'if(!rec && state.view!=="interf" && state.view!=="energy" && state.view!=="dre" && state.view!=="xps" && state.view!=="imgx" && state.view!=="lattice"){ renderWelcome(root); return; }')
# Diffraction 탭(키 lattice): 파일 없이도 열림(Phase DB·SAED 파일 열기), 분기 = DIFF_UI (FFT 모드 안에서 기존 renderLattice)
s,_n=re.subn(r'    if\(k==="lattice"&&!\(rec&&rec\.images&&rec\.images\.length\)\) dis="disabled";\n','',s); assert _n in (0,1)
assert 'k==="lattice"&&!(rec' not in s, '격자 탭 비활성 줄이 남음'
R('Diffraction 분기',r'    else if\(state\.view==="lattice"\) (?:renderLattice\(root\)|\(window\.DIFF_UI\?window\.DIFF_UI\.render\(root\):renderLattice\(root\)\));\n',
  '    else if(state.view==="lattice") (window.DIFF_UI?window.DIFF_UI.render(root):renderLattice(root));\n')
R('탭 분기',r'    else if\(state\.view==="interf"\) renderInterf\(root\);\n(?:    else if\(state\.view==="energy"\) window\.ED_UI\.render\(root\);\n)?(?:    else if\(state\.view==="dre"\) window\.DRE_UI\.render\(root\);\n)?(?:    else if\(state\.view==="xps"\) window\.XPS_UI\.render\(root\);\n)?(?:    else if\(state\.view==="imgx"\) window\.IMG_UI\.render\(root\);\n)?',
  '    else if(state.view==="interf") renderInterf(root);\n    else if(state.view==="energy") window.ED_UI.render(root);\n    else if(state.view==="dre") window.DRE_UI.render(root);\n    else if(state.view==="xps") window.XPS_UI.render(root);\n    else if(state.view==="imgx") window.IMG_UI.render(root);\n')
R('출처 줄·파일 열',r'  renderSidebar\(\);\n(?:  const _noFile=state\.view==="energy"\|\|state\.view==="dre"(?:\|\|state\.view==="xps")?(?:\|\|state\.view==="imgx")?;\n)?  renderProv\((?:rec|state\.view==="energy"\?null:rec|_noFile\?null:rec)\);\n(?:  \{ const a=document\.getElementById\("app"\); if\(a\) a\.classList\.toggle\("ed-on",(?:state\.view==="energy"|_noFile)\); \}\n)?',
  '  renderSidebar();\n  const _noFile=state.view==="energy"||state.view==="dre"||state.view==="xps"||state.view==="imgx";\n  renderProv(_noFile?null:rec);\n  { const a=document.getElementById("app"); if(a) a.classList.toggle("ed-on",_noFile); }\n')
R('간섭 검사 라인 추천 훅',r'(?:  const _run=\(\)=>\{ run\(\); if\(window\.IF_SUG\) IF_SUG\.after\(\); \};\n)?  \$\("#ifRun"\)\.onclick=(?:run|_run);\n(?:  if\(window\.IF_SUG\) IF_SUG\.setup\(root\);\n)?  if\(P\.els\.length>=(?:1|2)\) (?:run|_run)\(\);\n',
  '  const _run=()=>{ run(); if(window.IF_SUG) IF_SUG.after(); };\n  $("#ifRun").onclick=_run;\n  if(window.IF_SUG) IF_SUG.setup(root);\n  if(P.els.length>=1) _run();\n')
R('휴대폰 하단 탭',r'grid-template-columns:repeat\((?:6|7|8|9|10),1fr\);gap:0!important;','grid-template-columns:repeat(9,1fr);gap:0!important;')
# EM Image 탭: 영상이 없어도(스펙트럼만 있는 파일) 파일이 있으면 열림 → 그 안의 Metadata 보기로
R('EM Image 탭 활성 조건',r'if\(k==="image"\)\{ const n=rec\?rec\.images\.length:0; cnt=`<span class="cnt">\$\{n\}</span>`; if\(!(?:n|rec)\) dis="disabled"; \}',
  'if(k==="image"){ const n=rec?rec.images.length:0; cnt=`<span class="cnt">${n}</span>`; if(!rec) dis="disabled"; }')
# EM Image 안 보기 전환 막대 (Image · Metadata) + 탭 막대 맞춤(NAV_FIT: 한 줄 · 짧은 이름 · More)
R('EM Image 보기 전환',r'  try\{\n(?:    if\(rec&&\(state\.view==="image"\|\|state\.view==="meta"\)\)[^\n]*\n)?    if\(state\.view==="image"\) renderImage\(root\);\n',
  '  try{\n    if(rec&&(state.view==="image"||state.view==="meta")){ if(state.view==="image"&&!rec.images.length) state.view="meta"; const sw=h(`<div class="seg emsub" role="tablist" aria-label="EM Image view"><button data-sv="image" class="${state.view==="image"?"active":""}" ${rec.images.length?"":"disabled title=\\"No image in this file\\""}>Image</button><button data-sv="meta" class="${state.view==="meta"?"active":""}">Metadata</button></div>`); sw.querySelectorAll("button").forEach(b=>{ if(!b.disabled) b.onclick=()=>{ state._viewUser=b.dataset.sv; go(b.dataset.sv); }; }); root.appendChild(sw); }\n    if(state.view==="image") renderImage(root);\n')
R('출처 줄: 화소 크기 설정 링크 이름',r'<a href="#" id="provToLat">Set in (?:Lattice FFT|Diffraction)</a>','<a href="#" id="provToLat">Set in Diffraction</a>')
R('탭 막대 맞춤 호출',r'  nav\.querySelectorAll\("button"\)\.forEach\(b=>\{ if\(!b\.disabled\) b\.onclick=\(\)=>\{ state\._viewUser=b\.dataset\.v; go\(b\.dataset\.v\); if\(window\.matchMedia&&matchMedia\("\(max-width:720px\)"\)\.matches\) window\.scrollTo\(0,0\); \}; \}\);\n(?:  if\(window\.NAV_FIT\) NAV_FIT\(\);\n)?',
  lambda m:m.group(0).split('\n')[0]+'\n  if(window.NAV_FIT) NAV_FIT();\n')
R('전역 끌어놓기 → XPS·Image 탭이면 그 탭으로',r'if\(ev==="drop"\) ?(?:loadFiles\(e\.dataTransfer\.files\);|\{ if\(state\.view==="xps"&&window\.XPS_UI\) XPS_UI\._load\(e\.dataTransfer\.files\);(?: else if\(state\.view==="imgx"&&window\.IMG_UI\) IMG_UI\._load\(e\.dataTransfer\.files\);)? else loadFiles\(e\.dataTransfer\.files\); \})',
  'if(ev==="drop"){ if(state.view==="xps"&&window.XPS_UI) XPS_UI._load(e.dataTransfer.files); else if(state.view==="imgx"&&window.IMG_UI) IMG_UI._load(e.dataTransfer.files); else loadFiles(e.dataTransfer.files); }')
rd=lambda p:open(p,encoding='utf-8').read()
css=rd(f'{ED}/ed.css')+'\n'+rd(f'{DR}/dre.css')+'\n'+rd(f'{XP}/xps.css')+'\n'+rd(f'{IM}/img.css')+'\n'+rd(f'{DF}/diff.css')+('\n'+rd(f'{DF}/diff_db.css') if os.path.exists(f'{DF}/diff_db.css') else '')+'''
/* 파일 열을 쓰지 않는 탭(Energy Diagram·DRE·XPS — XPS 는 자체 파일 목록)에서는 숨기고 폭을 넘긴다 */
#app.ed-on .side,#app.ed-on .sidetog{display:none!important}
#app.ed-on .main{padding-left:0}
#nav button[data-v="energy"]{--ico:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M4 20h6M4 14h6M4 5h6M14 20h6M14 9h6'/%3E%3Cpath d='M7 18V8'/%3E%3Cpath d='M17 18v-6' stroke-dasharray='2 2'/%3E%3C/svg%3E")}
#nav button[data-v="dre"]{--ico:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M3 20h18'/%3E%3Cpath d='M4 19c5 0 7-1 9-5s3-9 7-10'/%3E%3Ccircle cx='8' cy='18' r='1'/%3E%3Ccircle cx='13' cy='14' r='1'/%3E%3Ccircle cx='17' cy='6' r='1'/%3E%3C/svg%3E")}
#nav button[data-v="xps"]{--ico:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M3 20h18'/%3E%3Cpath d='M3 18c3 0 4-1 5-6s2-8 3-8 2 6 3 9 2 5 7 5'/%3E%3Cpath d='M14 13c1-3 2-5 3-5'/%3E%3C/svg%3E")}
#nav button[data-v="imgx"]{--ico:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Crect x='3' y='3' width='18' height='18' rx='2'/%3E%3Ccircle cx='9' cy='9' r='2.2'/%3E%3Ccircle cx='15.5' cy='14' r='3'/%3E%3Ccircle cx='8' cy='16.5' r='1.3'/%3E%3C/svg%3E")}
@media (max-width:720px){ #nav button{font-size:9.5px!important}
  #nav button .nl{font-size:0;line-height:0;display:block;max-width:100%}
  #nav button .nl::after{content:attr(data-short);font-size:9.5px;line-height:1.2;white-space:nowrap;letter-spacing:-.2px} }
/* 탭 이름이 영어로 길어졌다 — 데스크톱에서 두 줄로 꺾이지 않도록 좌우 여백을 줄인다 */
@media (min-width:721px){ #nav button{padding-left:10px;padding-right:10px} }
@media (min-width:721px) and (max-width:1500px){ #nav button{font-size:13px;padding-left:8px;padding-right:8px} }
@media (max-width:380px){ #nav button .nl::after{font-size:8.5px;letter-spacing:-.3px} }
'''
head='<!-- INNO-ED:BEGIN css -->\n<style id="ed-css">\n'+css+'</style>\n<!-- INNO-ED:END -->\n'
# 공용 터치 제스처(window.INNO_TOUCH) — 모든 탭·나중의 Diffraction 이 같이 쓴다. head 에 두어 어느 탭 코드보다 먼저 정의된다.
tjs=rd(f'{TO}/inno_touch.js'); assert '</script' not in tjs.lower(), 'inno_touch: script 종료 태그'
head+='<!-- INNO-TOUCH:BEGIN css+js -->\n<style id="inno-touch-css">\n'+rd(f'{TO}/touch.css')+'</style>\n<script id="inno-touch-js">\n'+tjs+'</script>\n<!-- INNO-TOUCH:END -->\n'
NV=os.path.join(HERE,'nav'); njs=rd(f'{NV}/nav_fit.js'); assert '</script' not in njs.lower()
head+='<!-- INNO-NAV:BEGIN css+js -->\n<style id="inno-nav-css">\n'+rd(f'{NV}/nav_fit.css')+'</style>\n<script id="inno-nav-js">\n'+njs+'</script>\n<!-- INNO-NAV:END -->\n'
tjs=rd(f'{NV}/tips.js'); assert '</script' not in tjs.lower()
head+='<!-- INNO-TIPS:BEGIN css+js -->\n<style id="inno-tips-css">\n'+rd(f'{NV}/tips.css')+'</style>\n<script id="inno-tips-js">\n'+tjs+'</script>\n<!-- INNO-TIPS:END -->\n'
i=s.index('</head>'); s=s[:i]+head+s[i:]
edjs=''.join(rd(f'{ED}/{f}')+'\n' for f in ('ed_data.js','ed_core.js','ed_ui.js'))
core=rd(f'{DR}/dre_core.js'); ui=rd(f'{DR}/dre_ui.js'); ifjs=rd(f'{IF}/interf_suggest.js')
assert os.path.exists(f'{XP}/xps_lib.js'), 'xps_lib.js 없음: python3 -I modules/xps/build/build_lib.py 먼저'
xlib=rd(f'{XP}/xps_lib.js'); xio=rd(f'{XP}/xps_io.js'); xcore=rd(f'{XP}/xps_core.js'); xui=rd(f'{XP}/xps_ui.js')
# XPS 엔진(v1.8; v1.9.3 guide): core → guide → multi → arxps → image → proc 순서. 페이지 전역 + XPS 탭 Worker(xps_ui.js workerURL 가 이 태그들의 글을 이어 붙임)
xguide=rd(f'{XP}/xps_guide.js'); xmulti=rd(f'{XP}/xps_multi.js'); xarx=rd(f'{XP}/xps_arxps.js'); ximg=rd(f'{XP}/xps_image.js'); xproc=rd(f'{XP}/xps_proc.js')
icore=''.join(rd(os.path.join(IM,'core',f))+'\n' for f in sorted(os.listdir(os.path.join(IM,'core'))) if f.endswith('.js'))
iio=rd(f'{IM}/zstd_vendor.js')+'\n'+rd(f'{IM}/img_io.js')+'\n'+rd(f'{IM}/img_micro.js')+'\n'+rd(f'{IM}/img_afm.js')+'\n'+rd(f'{IM}/img_lif.js')+'\n'+rd(f'{IM}/img_oir.js'); iui='(function(){\n"use strict";\n'+''.join(rd(os.path.join(IM,'ui',f))+'\n' for f in sorted(os.listdir(os.path.join(IM,'ui'))) if f.endswith('.js'))+'})();\n'
DFS=[('diff-src-sg','diff_sg.js'),('diff-src-scat','diff_scat.js'),('diff-src-core','diff_core.js'),('diff-src-dyn','diff_dyn.js'),('diff-src-ocr','diff_ocr.js'),('diff-src-pattern','diff_pattern.js'),('diff-src-index','diff_index.js'),('diff-pack','diff_pack.js'),('diff-db','diff_db.js'),('diff-db-ui','diff_db_ui.js'),('diff-work','diff_work.js'),('diff-ui','diff_ui.js')]
dfjs=[(i,rd(f'{DF}/{fn}')) for i,fn in DFS if os.path.exists(f'{DF}/{fn}')]
for name,src in tuple(dfjs)+(('touch',tjs),('img_core',icore),('img_io',iio),('img_ui',iui),('ed',edjs),('dre_core',core),('dre_ui',ui),('interf',ifjs),('xps_lib',xlib),('xps_io',xio),('xps_core',xcore),('xps_guide',xguide),('xps_multi',xmulti),('xps_arxps',xarx),('xps_image',ximg),('xps_proc',xproc),('xps_ui',xui)):
    assert '</script' not in src.lower(), name+': script 종료 태그'
    # 문법 확인 (node 가 있으면): 주석이 코드를 삼키는 등의 편집 실수를 주입 전에 막음
    import shutil as _sh, subprocess as _sp, tempfile as _tf
    if _sh.which('node'):
        with _tf.NamedTemporaryFile('w',suffix='.js',delete=False,encoding='utf-8') as _f: _f.write(src); _tn=_f.name
        _r=_sp.run(['node','--check',_tn],capture_output=True,text=True); os.unlink(_tn)
        assert _r.returncode==0, name+': 문법 오류\n'+_r.stderr[:800]
body=('<!-- INNO-ED:BEGIN js -->\n<script id="ed-js">\n'+edjs+'</script>\n<!-- INNO-ED:END -->\n'
      '<!-- INNO-DRE:BEGIN js -->\n<script id="dre-core-src">\n'+core+'</script>\n<script id="dre-ui">\n'+ui+'</script>\n<!-- INNO-DRE:END -->\n'
      '<!-- INNO-IF:BEGIN js -->\n<script id="interf-suggest">\n'+ifjs+'</script>\n<!-- INNO-IF:END -->\n'
      '<!-- INNO-XPS:BEGIN js -->\n<script id="xps-lib">\n'+xlib+'</script>\n<script id="xps-io">\n'+xio+'</script>\n<script id="xps-core-src">\n'+xcore+'</script>\n<script id="xps-guide-src">\n'+xguide+'</script>\n<script id="xps-multi-src">\n'+xmulti+'</script>\n<script id="xps-arxps-src">\n'+xarx+'</script>\n<script id="xps-image-src">\n'+ximg+'</script>\n<script id="xps-proc-src">\n'+xproc+'</script>\n<script id="xps-ui">\n'+xui+'</script>\n<!-- INNO-XPS:END -->\n'
      '<!-- INNO-IMG:BEGIN js -->\n<script id="img-core-src">\n'+icore+'</script>\n<script id="img-io">\n'+iio+'</script>\n<script id="img-ui">\n'+iui+'</script>\n<!-- INNO-IMG:END -->\n'
      '<!-- INNO-DIFF:BEGIN js -->\n'+''.join('<script id="'+i+'">\n'+t+'</script>\n' for i,t in dfjs)+'<!-- INNO-DIFF:END -->\n')
k=s.index('<script>window.INNO_UX') if '<script>window.INNO_UX' in s else s.rindex('</body>')
s=s[:k]+body+s[k:]
open(outp,'w',encoding='utf-8').write(s)
print('ok',len(s))
