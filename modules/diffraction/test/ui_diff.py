"""Diffraction 탭 UI 시험 (Playwright · Chromium). 사용: python3 ui_diff.py <주입된 html> <그림 폴더> [실측 폴더]
실측 폴더(기본 ../ref/real): GOLD-NP-DIFF.dm3 (Au 링), Si001_…dm3 (Si [001] 스폿), FeCr_ODS_SAD_….tif (8-bit, 스케일바) — 없으면 그 항목은 건너뜀.
상 목록: DIFF_PACK 이 없으면 test/cif 의 시험 구조로 임시 묶음을 만든다. 데스크톱 1400×900 + 휴대폰 390×844 (가로 스크롤 없음, 콘솔 오류 0).
선택자는 id / data-* 만."""
import sys, os, glob, json
from playwright.sync_api import sync_playwright
html=sys.argv[1]; out=sys.argv[2]; os.makedirs(out,exist_ok=True)
HERE=os.path.dirname(os.path.abspath(__file__)); REAL=sys.argv[3] if len(sys.argv)>3 else os.path.join(HERE,'..','ref','real')
CIFS={os.path.basename(p)[:-4]:open(p).read() for p in glob.glob(os.path.join(HERE,'cif','*.cif'))}
res=[]; errs=[]
def chk(c,name,info=''): res.append((bool(c),name)); print(('PASS' if c else 'FAIL'),name,info)
PACK='''(cifs)=>{ if(window.DIFF_PACK&&DIFF_PACK.phases&&DIFF_PACK.phases.length) return DIFF_PACK.phases.length; const DC=DIFF_CORE; const L=[];
  for(const [k,t] of Object.entries(cifs)){ const p=DC.parseCIF(t)[0]; L.push(Object.assign(DC.recordFromPhase(p,{id:"pack:"+k,name:k}),{layer:"pack"})); }
  window.DIFF_PACK={version:"test",phases:L}; return L.length; }'''
def wait_res(pg,ms=60000):
    pg.wait_for_function('()=>{ const S=DIFF_UI._S(); return S&&S.res&&!S.busy&&!S.idBusy; }',timeout=ms); pg.wait_for_timeout(400)
def hscroll(pg): return pg.evaluate('document.documentElement.scrollWidth>document.documentElement.clientWidth+1')
def run(pg,tag):
    pg.goto('file://'+os.path.abspath(html)); pg.wait_for_timeout(1200)
    n=pg.evaluate(PACK,CIFS); chk(n>=5,f'{tag}: phase list available',n)
    # 파일 없이 Diffraction 탭 → SAED 열기 안내
    pg.click('#nav button[data-v="lattice"]'); pg.wait_for_timeout(500)
    chk(pg.query_selector('#dfOpen') is not None and pg.evaluate('state.view')=='lattice',f'{tag}: Diffraction tab opens without a file (Open file…)')
    pg.screenshot(path=f'{out}/{tag}_empty.png')
    au=os.path.join(REAL,'GOLD-NP-DIFF.dm3')
    if os.path.exists(au):
        pg.set_input_files('#fin',au); pg.wait_for_timeout(800); pg.click('#nav button[data-v="lattice"]'); wait_res(pg)
        S=pg.evaluate('(()=>{ const S=DIFF_UI._S(); return {type:S.type,c:S.res.rings&&S.res.rings.center,ratio:S.res.rings&&S.res.rings.ratio,top:S.rcands&&S.rcands[0]&&{n:S.rcands[0].phaseName,k:S.rcands[0].kappa,s:S.rcands[0].scale},mode:document.querySelector("#dfMode .active").dataset.m}; })()')
        print(S)
        chk(S['mode']=='saed','Au: SAED mode chosen from the 1/nm calibration',S['mode'])
        chk(S['type']=='rings','Au: ring pattern detected',S['type'])
        top2=pg.evaluate('DIFF_UI._S().rcands.slice(0,2).map(c=>c.phaseName)'); tie=pg.evaluate('!!document.querySelector("#dfResult .note.warn")')
        FCC={'Gold','Aluminium','Copper','Nickel','Platinum','Silver'}
        chk(all(n in FCC for n in top2),'Au: without elements the top 2 are fcc metals (same ring ratios — honest ambiguity)',top2)
        chk(tie or top2[0] in ('Gold',) or top2[0].startswith('Au'),'Au: tie note shown when an fcc look-alike is next',(top2,tie))
        pg.fill('#dfEls','Au'); pg.dispatch_event('#dfEls','change'); pg.wait_for_timeout(300); wait_res(pg)
        S=pg.evaluate('(()=>{ const S=DIFF_UI._S(); return {top:S.rcands&&S.rcands[0]&&{n:S.rcands[0].phaseName,k:S.rcands[0].kappa,s:S.rcands[0].scale}}; })()')
        chk(S['top'] and (S['top']['n'].startswith('Au') or S['top']['n']=='Gold'),'Au: with Elements = Au the top ring candidate is Au',S['top'])
        chk(S['top'] and abs(S['top']['s']/0.011116446-1.089)<0.01,'Au: implied scale ÷ DM calibration ≈ 1.089 (DM tag off by ~9 %), ranked without trusting the tag',S['top'])
        pg.screenshot(path=f'{out}/{tag}_au.png',full_page=True); chk(not hscroll(pg),f'{tag} Au: no page horizontal scroll')
        # 보정 저장 (prompt/confirm 자동 응답)
        pg.once('dialog',lambda d: d.accept('Test 375mm'))
        if pg.query_selector('#dfCalib'):
            pg.click('#dfCalib'); pg.wait_for_timeout(800)
            cal=pg.evaluate('(()=>{ try{ return JSON.parse(localStorage.getItem("inno.diff.cal")||"[]"); }catch(e){ return []; } })()')
            got=pg.evaluate('window.DIFF_DB&&DIFF_DB.listCalibrations?1:0')
            chk(got or (cal and abs(cal[-1]['scale']-0.01211)/0.01211<0.01),'Au: Use as calibration stores scale ≈ 0.0121 1/nm/px',cal[-1] if cal else got)
            # 저장한 보정으로 다시 읽기 → κ ≈ 1 (단축 기준 반지름과 평균 축척을 섞으면 √e 만큼 어긋남)
            pg.click('#dfCal button[data-c="saved"]'); pg.wait_for_timeout(800)
            opts=pg.eval_on_selector_all('#dfSaved option','els=>els.map(e=>[e.value,e.textContent])'); pick=[o for o in opts if 'Test 375mm' in o[1]]
            if pick:
                pg.select_option('#dfSaved',pick[-1][0]); pg.wait_for_timeout(300); pg.select_option('#dfRankCal','on'); pg.wait_for_timeout(300); wait_res(pg)
                k=pg.evaluate('(()=>{ const S=DIFF_UI._S(); const c=S.rcands[0]; const cs=S.savedCal; return {n:c.phaseName,s:c.scale,mean:cs.scale,e:cs.ellip,ratio:S.res.rings.ratio}; })()')
                kap=k['s']/(k['mean']*k['ratio']**0.5); chk(abs(kap-1)<0.002 and (k['n']=='Gold' or k['n'].startswith('Au')),'Au: re-read with the saved calibration → κ = 1.000 and Au first (calibration in ranking)',(kap,k))
            else: chk(False,'Au: saved calibration listed',opts)
            pg.select_option('#dfRankCal','auto'); pg.click('#dfCal button[data-c="auto"]'); pg.wait_for_timeout(300)
    si=glob.glob(os.path.join(REAL,'Si001*.dm3'))
    if si:
        pg.evaluate('DIFF_UI._P.elements=""')
        pg.set_input_files('#fin',si[0]); pg.wait_for_timeout(800); pg.click('#nav button[data-v="lattice"]'); pg.wait_for_timeout(300)
        mode=pg.evaluate('document.querySelector("#dfMode .active").dataset.m'); chk(mode=='saed','Si: SAED chosen by file name although the tag unit is nm',mode)
        wait_res(pg)
        S=pg.evaluate('(()=>{ const S=DIFF_UI._S(); const L=S.res.spots&&S.res.spots.lattices[0]; return {type:S.type,len:L&&L.len,ang:L&&L.ang12,n:L&&L.nIndexed,warn:document.querySelector("#dfBody .note.warn")&&document.querySelector("#dfBody .note.warn").textContent}; })()')
        print(S); chk(S['type']=='spots' and S['len'] and abs(S['len'][0]-336.9)<1 and abs(S['len'][1]-361.4)<1,'Si: spot lattice 337/361 px',S)
        chk(S['warn'] and 'nm' in S['warn'],'Si: warns that the file calibration is real-space (nm)',S['warn'])
        # 보정 없이 상 찾기, 찌그러짐 허용 8 % → Si [001] 이 후보 안
        pg.fill('#dfDist','8'); pg.dispatch_event('#dfDist','change'); pg.wait_for_timeout(300); wait_res(pg)
        C=pg.evaluate('DIFF_UI._S().cands.slice(0,8).map(c=>[c.phaseName,c.zone.join(","),+c.score.toFixed(2)])'); print(C)
        chk(any(c[0].startswith('Si') and sorted(abs(int(v)) for v in c[1].split(','))==[0,0,1] for c in C),'Si: Si [001] among the top 8 (distortion 8 %)',C)
        pg.click('#dfTool button[data-t="center"]'); pg.wait_for_timeout(150)
        bx=pg.eval_on_selector('#dfCv','e=>{ const r=e.getBoundingClientRect(); return [r.left,r.top,r.width,r.height]; }')
        c0=pg.evaluate('DIFF_UI._S().res.spots.lattices[0].center'); v=pg.evaluate('DIFF_UI._S().view')
        pg.mouse.click(bx[0]+(c0['x']+12-v['x'])*v['zoom'],bx[1]+(c0['y']-9-v['y'])*v['zoom']); pg.wait_for_timeout(500); wait_res(pg)
        L2=pg.evaluate('(()=>{ const L=DIFF_UI._S().res.spots.lattices[0]; return L&&{c:L.center,len:L.len}; })()')
        chk(L2 and abs(L2['c']['x']-c0['x'])<1.5 and abs(L2['c']['y']-c0['y'])<1.5 and abs(L2['len'][0]-336.9)<1,'Si: Set center 15 px off snaps to the Friedel centre and refits',L2)
        pg.screenshot(path=f'{out}/{tag}_si.png',full_page=True); chk(not hscroll(pg),f'{tag} Si: no page horizontal scroll')
    fe=glob.glob(os.path.join(REAL,'FeCr*.tif'))
    if fe:
        pg.set_input_files('#fin',fe[0]); pg.wait_for_timeout(800); pg.click('#nav button[data-v="lattice"]'); pg.wait_for_timeout(300)
        mode=pg.evaluate('document.querySelector("#dfMode .active").dataset.m'); chk(mode=='saed','FeCr tif: SAED chosen',mode)
        if mode=='saed':
            wait_res(pg); bar=pg.evaluate('DIFF_UI._S().res.scaleBar'); chk(bar and abs(bar['lenPx']-177)<=1,'FeCr: scale bar 177 px found',bar)
            lab=pg.evaluate('(()=>{ const S=DIFF_UI._S(); return {lab:S.res.scaleLabel,val:S.barVal,ocr:!!S.barFromOCR}; })()')
            chk(lab['lab'] and lab['lab']['text']=='5 1/nm' and lab['val']==5 and lab['ocr'],'FeCr: scale label "5 1/nm" read and prefilled (not "SAD 1")',lab)
            pg.click('#dfCal button[data-c="bar"]'); pg.wait_for_timeout(200); pg.fill('#dfBarVal','5'); pg.dispatch_event('#dfBarVal','change'); pg.wait_for_timeout(300)
            sc=pg.evaluate('DIFF_UI._S().barVal'); chk(sc==5,'FeCr: bar value accepted',sc)
            pg.select_option('#dfTopN','16'); pg.wait_for_timeout(300); wait_res(pg)
            L=pg.evaluate('(()=>{ const L=DIFF_UI._S().res.spots.lattices[0]; return L&&[L.len[0]*5/177,L.len[1]*5/177,L.ang12]; })()'); print(L)
            chk(L and abs(L[0]-8.38)<0.1 and abs(L[1]-8.45)<0.1,'FeCr: strongest-16 lattice = bcc 211-type (8.4 1/nm)',L)
            pg.screenshot(path=f'{out}/{tag}_fecr.png',full_page=True); chk(not hscroll(pg),f'{tag} FeCr: no page horizontal scroll')
    # Phase DB 모드
    pg.click('#dfMode button[data-m="db"]'); pg.wait_for_timeout(600); pg.screenshot(path=f'{out}/{tag}_db.png',full_page=True)
    chk(pg.evaluate('document.querySelector("#dfMode .active").dataset.m')=='db',f'{tag}: Phase DB mode renders')
with sync_playwright() as p:
    b=p.chromium.launch()
    for tag,kw in [('desk',dict(viewport={'width':1400,'height':900})),('phone',dict(viewport={'width':390,'height':844},has_touch=True,is_mobile=True,device_scale_factor=2))]:
        pg=b.new_page(**kw)
        pg.on('pageerror',lambda e: errs.append('pageerror: '+str(e)))
        pg.on('console',lambda m: errs.append('console: '+m.text) if m.type=='error' and 'ERR_' not in m.text and 'gstatic' not in m.text else None)
        run(pg,tag); pg.close()
    chk(not errs,'no console/page errors',errs[:5])
    b.close()
print(sum(1 for r in res if r[0]),'/',len(res),'pass')
sys.exit(0 if all(r[0] for r in res) else 1)
