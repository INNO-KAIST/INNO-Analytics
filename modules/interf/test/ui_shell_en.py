"""셸·EM 탭 영어 골격 + 간섭 검사 추천 라인 UI 시험. 사용: python3 ui_shell_en.py <주입된 html> <그림 폴더>
(합성 스펙트럼·맵 레코드를 끼워 넣어 Spectrum/Maps 탭까지 그린다. 선택자는 data-v / id 만 쓴다)"""
import sys, os, json, math
from playwright.sync_api import sync_playwright
html=sys.argv[1]; out=sys.argv[2]; os.makedirs(out,exist_ok=True)
TIF=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','image','test','io','t16.tif')
errs=[]; res=[]
def chk(c,name,info=''): res.append((bool(c),name,str(info))); print(('PASS' if c else 'FAIL'),name,info)
SYN='''(()=>{
  const rec=curRec(); const orig=rec.h5.get.bind(rec.h5);
  const n=2048, disp=10, off=0; const y=new Float64Array(n);
  const peaks=[[0.677,900],[1.04,300],[1.74,200],[8.04,700],[1.5,150],[9.713,600],[7.18,250],[6.0,150]];
  for(let i=0;i<n;i++){ const E=(off+disp*i)/1000; let v=40*Math.exp(-E/3)+2; for(const [pe,a] of peaks){ const sg=(0.055+0.002*pe)/2.355*1.0; v+=a*Math.exp(-0.5*((E-pe)/sg)**2); } y[i]=Math.round(v*(0.9+0.2*((i*2654435761>>>0)%1000)/1000)); }
  const store={"/syn/spectrum/0/Data":{shape:[n,1],dtype:"<Q",value:y}};
  rec.h5={get(p){ const d=store[p]; if(d) return {shape:d.shape,dtype:d.dtype,value:d.value,slice(){return d.value;}}; return orig(p); }, close(){}};
  const md={Detectors:{"Detector-0":{DetectorType:"AnalyticalDetector",DetectorName:"Synthetic EDS",Dispersion:String(disp),OffsetEnergy:String(off),BeginEnergy:"150",FwhmMnKa:"128"}}};
  rec.spectra.push({kind:"spectrum",uid:"/syn/spectrum/0",base:"/syn/spectrum/0",shape:[n,1],dtype:"<Q",md,px:null,det:"Synthetic EDS",role:"Summed"});
  // fake Velox element maps so the Maps tab is capable
  const im0=rec.images[0];
  for(const t of ["Lu","Yb","F"]) rec.images.push(Object.assign({},im0,{uid:"m_"+t,isMap:true,title:t}));
  rec.bcfRegions=[{name:"Lu",z:71,line:"L",E:7.655},{name:"Yb",z:70,line:"L",E:7.415}];
  state.curSpec=rec.spectra.length-1; state.curFile=state.files.indexOf(rec); renderView();
})()'''
with sync_playwright() as p:
    b=p.chromium.launch()
    pg=b.new_page(viewport={'width':1400,'height':900})
    pg.on('pageerror',lambda e: errs.append('pageerror: '+str(e)))
    pg.on('console',lambda m: errs.append('console: '+m.text) if m.type=='error' and 'ERR_' not in m.text else None)
    pg.goto('file://'+os.path.abspath(html)); pg.wait_for_timeout(1500)
    labels=pg.eval_on_selector_all('#nav button','els=>els.map(e=>[e.dataset.v,e.textContent.trim(),e.getAttribute("aria-label"),e.querySelector(".nl")&&e.querySelector(".nl").dataset.short])')
    print(labels)
    exp=['XPS','EM Image','Image Analysis','EDS Maps','EDS Spectrum','Diffraction','EDS Interference','EM Metadata','Energy Diagram','DRE Kinetics']
    chk([l[1].rstrip('0123456789') for l in labels]==exp,'nav shows 10 English names',[l[1] for l in labels])
    chk([l[0] for l in labels]==['xps','image','imgx','maps','spectrum','lattice','interf','meta','energy','dre'],'view keys unchanged')
    pg.screenshot(path=f'{out}/d_welcome.png')
    # open an EM file (TIFF)
    pg.set_input_files('#fin',TIF); pg.wait_for_timeout(1500)
    chk(pg.evaluate('state.files.length')==1,'EM file (tif) opened', pg.evaluate('state.files[0].name'))
    for v in ['image','lattice','interf','meta']:
        pg.click(f'#nav button[data-v="{v}"]'); pg.wait_for_timeout(700)
        pg.screenshot(path=f'{out}/d_{v}.png',full_page=True)
        chk(pg.evaluate('state.view')==v,f'tab {v} renders')
    # interf with elements -> run + suggestion
    pg.click('#nav button[data-v="interf"]'); pg.wait_for_timeout(300)
    pg.fill('#ifEls','Na, F, Yb, Lu, Tm, Cu, C, O'); pg.click('#ifRun'); pg.wait_for_timeout(1500)
    pg.screenshot(path=f'{out}/d_interf_run.png',full_page=True)
    chk(pg.query_selector('#ifSug') is not None,'interf suggestion block rendered')
    chk('All overlaps' in pg.inner_text('#ifOut'),'interf all-overlaps summary (heading regex matches English)', pg.evaluate("(document.querySelector('details.if-all summary')||{}).textContent"))
    # synthetic spectrum + maps
    pg.evaluate(SYN); pg.wait_for_timeout(500)
    for v in ['spectrum','maps']:
        pg.click(f'#nav button[data-v="{v}"]'); pg.wait_for_timeout(1000)
        chk(pg.evaluate('state.view')==v,f'tab {v} renders (synthetic)')
    pg.click('#nav button[data-v="spectrum"]'); pg.wait_for_timeout(500)
    pg.fill('#inEls','Na, F, Yb, Lu, Tm'); pg.click('#btnAutoId'); pg.wait_for_timeout(2500)
    pg.screenshot(path=f'{out}/d_spectrum.png',full_page=True)
    pg.click('#nav button[data-v="maps"]'); pg.wait_for_timeout(800)
    pg.screenshot(path=f'{out}/d_maps.png',full_page=True)
    # k-set dialog
    pg.click('#nav button[data-v="spectrum"]'); pg.wait_for_timeout(500)
    if pg.query_selector('#qLib'):
        pg.click('#qLib'); pg.wait_for_timeout(400); pg.screenshot(path=f'{out}/d_kset.png'); 
        ok=pg.query_selector('#kfOk'); 
        if ok: ok.click()
        pg.wait_for_timeout(200)
        if pg.query_selector('#qStd'):
            pg.click('#qStd'); pg.wait_for_timeout(400); pg.screenshot(path=f'{out}/d_kstd.png'); pg.click('#stdCancel')
    for v in ['xps','imgx','energy','dre']:
        pg.click(f'#nav button[data-v="{v}"]'); pg.wait_for_timeout(800)
        pg.screenshot(path=f'{out}/d_{v}.png'); chk(pg.evaluate('state.view')==v,f'tab {v} renders')
    # lattice with spots on the real image: sidebar present
    pg.click('#nav button[data-v="lattice"]'); pg.wait_for_timeout(800); pg.screenshot(path=f'{out}/d_lattice2.png',full_page=True)
    # ---- mobile
    pg.set_viewport_size({'width':390,'height':844}); pg.wait_for_timeout(600)
    pg.click('#nav button[data-v="image"]'); pg.wait_for_timeout(500)
    pg.screenshot(path=f'{out}/m_image.png')
    info=pg.evaluate('''()=>[...document.querySelectorAll('#nav button')].map(b=>{const r=b.getBoundingClientRect(); return {v:b.dataset.v,w:Math.round(r.width*10)/10,sw:b.scrollWidth,cw:b.clientWidth,txt:getComputedStyle(b.querySelector('.nl'),'::after').content,h:Math.round(r.height)}})''')
    print(json.dumps(info,ensure_ascii=False))
    chk(all(i['sw']<=i['cw']+0 for i in info),'mobile nav: no label overflow at 390px',[ (i['v'],i['sw'],i['cw']) for i in info if i['sw']>i['cw']])
    chk(pg.evaluate('document.documentElement.scrollWidth<=window.innerWidth'),'mobile: no horizontal page scroll')
    pg.screenshot(path=f'{out}/m_nav.png',clip={'x':0,'y':844-70,'width':390,'height':70})
    for v in ['lattice','interf','meta','spectrum','maps']:
        pg.click(f'#nav button[data-v="{v}"]'); pg.wait_for_timeout(600); pg.screenshot(path=f'{out}/m_{v}.png',full_page=False)
    pg.set_viewport_size({'width':360,'height':780}); pg.wait_for_timeout(400)
    info2=pg.evaluate('''()=>[...document.querySelectorAll('#nav button')].map(b=>[b.dataset.v,b.scrollWidth,b.clientWidth])''')
    print('360px',info2); pg.screenshot(path=f'{out}/m360_nav.png',clip={'x':0,'y':780-70,'width':360,'height':70})
    b.close()
print('ERRORS:',errs)
chk(not errs,'no console/page errors',errs[:5])
print(sum(1 for r in res if r[0]),'/',len(res),'pass')
