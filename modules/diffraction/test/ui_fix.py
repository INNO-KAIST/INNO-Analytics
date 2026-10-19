"""감사(2026-10-19) 결함의 UI 회귀 시험 (Playwright · Chromium). 사용: python3 ui_fix.py <주입된 html> <출력 폴더>
시험 영상은 여기서 합성한다 (실측·내부 표본 없음): 8-bit 로그 표시 Au 링 + "5 1/nm" 이름표 · 16-bit 저계수 Si 스폿 · HRTEM 비슷한 8-bit 무늬.
확인: 모드 자동(SAED/FFT) · 링 판정 · 이름표 단위 확인(10 배 틀린 값 자동 사용 안 함) · 파일마다 보정 출처·패턴 Auto · 필터/κ 칩 ·
카메라 길이 보정 · 얼룩 격자 거부(탭 안 죽음) · 수동 hkl 경고 · 링 중심 유지 · Cancel · 8-bit 안내 · Phase ID 머리글 · 가로 스크롤·콘솔 오류 없음"""
import sys, os, json, math
import numpy as np, tifffile
from PIL import Image, ImageDraw, ImageFont
from playwright.sync_api import sync_playwright
html=sys.argv[1]; out=sys.argv[2]; os.makedirs(out,exist_ok=True)
res=[]; errs=[]
def chk(c,name,info=''): res.append((bool(c),name)); print(('PASS' if c else 'FAIL'),name,info if info!='' else '')

# ---------- 합성 영상 ----------
rg=np.random.default_rng(3)
def au_rings(w=1024,cx=505,cy=520,scale=0.012,peak=300,bg=10):
    a=0.40782; y,x=np.mgrid[0:w,0:w]; r=np.hypot(x-cx,y-cy); I=bg*(1+30*np.exp(-r/40))
    for (h,k,l),m,F in [((1,1,1),8,1.0),((2,0,0),6,0.8),((2,2,0),12,0.55),((3,1,1),24,0.45),((2,2,2),8,0.4),((4,0,0),6,0.3),((3,3,1),24,0.28),((4,2,0),24,0.26)]:
        g=math.sqrt(h*h+k*k+l*l)/a; rp=g/scale; A=peak*F*m/(g*g)*4*math.exp(-g*g/144); I=I+A*np.exp(-(r-rp)**2/8)
    return rg.poisson(I).astype(np.float64)
def log8(a):
    b=np.log1p(a); lo,hi=np.percentile(b,0.1),np.percentile(b,99.9); return np.clip((b-lo)/(hi-lo)*255,0,255).astype(np.uint8)
F=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',28)
def with_bar(a8,val,unit,scale):
    im=Image.fromarray(a8); d=ImageDraw.Draw(im); L=int(round(val/scale)); x0,y0=30,a8.shape[0]-40
    d.rectangle([x0,y0,x0+L-1,y0+7],fill=255); d.text((x0,y0-36),f"{val} {unit}",fill=255,font=F); return np.array(im)
def si_spots(w=512,peak=200):  # Si [011] 비슷한 스폿 (16-bit, 저계수: 최댓값 ≤ 255)
    y,x=np.mgrid[0:w,0:w]; cx,cy=250.3,261.7; I=2+6*np.exp(-((x-cx)**2+(y-cy)**2)/(2*60**2)); a1=np.array([96.0,40.0]); a2=np.array([-62.0,110.0])
    for i in range(-4,5):
        for j in range(-4,5):
            if (i+j)%2: continue
            p=np.array([cx,cy])+i*a1/2+j*a2/2
            if not(6<p[0]<w-6 and 6<p[1]<w-6): continue
            A=peak*math.exp(-0.06*(i*i+j*j)); I=I+A*np.exp(-((x-p[0])**2+(y-p[1])**2)/(2*1.6**2))
    return np.clip(rg.poisson(I),0,255).astype(np.uint16)
P=os.path.join(out,'in'); os.makedirs(P,exist_ok=True)
AU8=os.path.join(P,'au_log8_bar.tif'); tifffile.imwrite(AU8,with_bar(log8(au_rings()),5,'1/nm',0.012))
AUOFF=os.path.join(P,'au_off_log8.tif'); tifffile.imwrite(AUOFF,log8(au_rings(cx=230,cy=300)))
SI16=os.path.join(P,'si_lowcount16.tif'); tifffile.imwrite(SI16,si_spots())
yy,xx=np.mgrid[0:512,0:512]; HR=os.path.join(P,'hrtem_like.tif'); tifffile.imwrite(HR,np.clip(128+40*np.sin(xx*0.6)*np.sin(yy*0.5)+rg.normal(0,10,(512,512)),0,255).astype(np.uint8))

def wait_res(pg,ms=120000):
    pg.wait_for_function('()=>{ const S=DIFF_UI._S(); return S&&(S.res||S.err)&&!S.busy&&!S.idBusy; }',timeout=ms); pg.wait_for_timeout(300)
def openf(pg,f):
    pg.set_input_files('#fin',f); pg.wait_for_timeout(900); pg.click('#nav button[data-v="lattice"]'); pg.wait_for_timeout(400)
def mode(pg): return pg.evaluate('(document.querySelector("#dfMode .active")||{dataset:{}}).dataset.m')
def hscroll(pg): return pg.evaluate('document.documentElement.scrollWidth>document.documentElement.clientWidth+1')
def S(pg,expr): return pg.evaluate('(()=>{ const S=DIFF_UI._S(), P=DIFF_UI._P; return '+expr+'; })()')

def run(pg,tag):
    pg.goto('file://'+os.path.abspath(html)); pg.wait_for_timeout(1500)
    # 지난 세션 설정이 남아 있다고 가정: Scale bar · Rings · κ 98 %
    pg.evaluate('()=>{ const P=DIFF_UI._P; P.calSrc="manual"; P.pattern="spots"; P.elements=""; P.kappaPct=100; }')
    # 1. 8-bit 로그 SAED 사본 (교정 없음) → SAED, 링, 이름표 단위 확인
    openf(pg,AU8); m=mode(pg); chk(m=='saed',f'{tag}: 8-bit log ring export opens in SAED (not FFT)',m)
    if m!='saed': pg.click('#dfMode button[data-m="saed"]'); pg.wait_for_timeout(300)
    wait_res(pg)
    st=S(pg,'{cal:P.calSrc,pat:P.pattern,type:S.type,why:S.res.typeWhy,lab:S.res.scaleLabel,bar:S.barVal,sug:S.barSug,eight:S.meta.eightBit}'); print(json.dumps(st,ensure_ascii=False,default=str)[:600])
    chk(st['cal']=='auto' and st['pat']=='auto',f'{tag}: new file resets calibration source and pattern to Auto',(st['cal'],st['pat']))
    chk(st['type']=='rings',f'{tag}: 8-bit log ring export classified as rings',st['why'])
    lab=st['lab'] or {}
    chk(st['bar'] is None or abs(st['bar']-5)<1e-6,f'{tag}: bar value never pre-filled 10x wrong (label {lab.get("text")} conf {lab.get("conf")})',st['bar'])
    if st['sug']:
        btn=pg.query_selector('[data-baruse="5"]'); chk(btn is not None,f'{tag}: unsure unit → both readings offered as buttons',st['sug'])
        if btn: btn.click(); pg.wait_for_timeout(300); wait_res(pg); chk(abs(S(pg,'S.barVal')-5)<1e-9,f'{tag}: choosing "5 1/nm" sets the bar value')
    # Phase ID 머리글: 보정을 순위에 쓰지 않으면 '보정값 +' 라고 하지 않음
    hd=pg.evaluate('[...document.querySelectorAll("#dfResult h3")].map(h=>h.textContent).join(" | ")'); cal=S(pg,'DIFF_UI._calib()')
    chk(cal['rank'] or '보정값 +' not in hd,f'{tag}: Phase ID header honest about calibration use',(hd[:120],cal['rank'],cal['prior']))
    # 칩: Elements · κ
    pg.fill('#dfEls','Au'); pg.dispatch_event('#dfEls','change'); wait_res(pg)
    chk(pg.query_selector('#dfElChip') is not None,f'{tag}: Elements filter chip visible')
    pg.fill('#dfKappa','98'); pg.dispatch_event('#dfKappa','change'); wait_res(pg)
    chk(pg.query_selector('#dfKChip') is not None,f'{tag}: κ ≠ 100 % chip visible')
    pg.click('#dfKChip button'); wait_res(pg); pg.click('#dfElChip button'); wait_res(pg)
    chk(S(pg,'P.kappaPct')==100 and S(pg,'P.elements')=='' and pg.query_selector('#dfKChip') is None,f'{tag}: chips clear κ and Elements')
    # 카메라 길이 보정
    pg.click('#dfCal button[data-c="camlen"]'); pg.wait_for_timeout(300)
    for k,v in [('#dfCamL','300'),('#dfCamKV','200'),('#dfCamPx','14'),('#dfCamBin','1')]: pg.fill(k,v); pg.dispatch_event(k,'change'); pg.wait_for_timeout(150)
    wait_res(pg); cal=S(pg,'DIFF_UI._calib()'); lam=1.226426/math.sqrt(200e3*(1+0.97848e-6*200e3)); want=14e-3/(lam*300)
    chk(cal['src']=='camlen' and cal['scale'] and abs(cal['scale']/want-1)<1e-6 and cal['prior'],f'{tag}: camera-length calibration = pixel/(λL) = {want:.5f} 1/nm/px, used as ±15 % prior',cal)
    pg.select_option('#dfUnit','A'); pg.wait_for_timeout(300)
    chk('1/Å' in pg.evaluate('document.querySelector("#dfBody").innerText'),f'{tag}: Å⁻¹ display units selectable')
    pg.select_option('#dfUnit','nm'); pg.wait_for_timeout(200); pg.click('#dfCal button[data-c="auto"]'); pg.wait_for_timeout(200); wait_res(pg)
    # 링 중심 직접 지정 → Analyze 해도 유지 → Auto center 로 지움
    pg.click('#dfTool button[data-t="center"]'); pg.wait_for_timeout(150)
    bx=pg.eval_on_selector('#dfCv','e=>{ const r=e.getBoundingClientRect(); return [r.left,r.top]; }'); v=S(pg,'S.view'); c=S(pg,'S.res.rings.center')
    pg.mouse.click(bx[0]+(c['x']+3-v['x'])*v['zoom'],bx[1]+(c['y']-2-v['y'])*v['zoom']); pg.wait_for_timeout(400); wait_res(pg)
    rc=S(pg,'S.ringCenter'); pg.click('#dfRun'); pg.wait_for_timeout(300); wait_res(pg)
    chk(rc and S(pg,'S.ringCenter') is not None and pg.query_selector('#dfAutoC') is not None,f'{tag}: manual ring centre kept on Analyze',rc)
    pg.click('#dfAutoC'); pg.wait_for_timeout(300); wait_res(pg); chk(S(pg,'S.ringCenter') is None,f'{tag}: Auto center clears the manual centre')
    # 링을 스폿으로 강제 → 얼룩 격자 거부 (탭이 멈추지 않음)
    pg.click('#dfPat button[data-p="spots"]'); pg.wait_for_timeout(300); wait_res(pg,60000)
    st=S(pg,'{type:S.type,info:S.idInfo,n:(S.cands||[]).length,err:S.err}'); print(st)
    chk(st['type']=='spots' and not st['err'] and ((st['info'] or {}).get('refused')=='speckle' or (st['info'] or {}).get('noLattice')),f'{tag}: rings forced to Spots → speckle lattice refused, no crash',st)
    pg.screenshot(path=f'{out}/{tag}_au8.png',full_page=True); chk(not hscroll(pg),f'{tag} au8: no page horizontal scroll')
    # 2. 치우친 8-bit 링 → rings, 새 파일이라 패턴 Auto
    openf(pg,AUOFF)
    if mode(pg)!='saed': pg.click('#dfMode button[data-m="saed"]'); pg.wait_for_timeout(300)
    wait_res(pg); st=S(pg,'{pat:P.pattern,type:S.type,why:S.res.typeWhy}'); chk(st['pat']=='auto' and st['type']=='rings',f'{tag}: off-centre 8-bit rings → Auto + rings',st)
    # 3. 16-bit 저계수 스폿: 8-bit 안내 없음 · 수동 hkl 경고 · Cancel
    openf(pg,SI16)
    if mode(pg)!='saed': pg.click('#dfMode button[data-m="saed"]'); pg.wait_for_timeout(300)
    wait_res(pg); txt=pg.evaluate('document.querySelector("#dfBody").innerText')
    chk('8비트' not in txt and S(pg,'S.meta.eightBit')==False,f'{tag}: no "8-bit display copy" note on low-count 16-bit data',S(pg,'{bits:S.meta.bits,max:S.meta.rawMax}'))
    st=S(pg,'{type:S.type,L:S.res.spots&&S.res.spots.lattices[0]&&S.res.spots.lattices[0].len}'); print(st)
    if pg.query_selector('#dfMGo'):
        pg.wait_for_timeout(500); pg.fill('#dfMh1','1 0 0'); pg.fill('#dfMh2','0 1 0'); pg.click('#dfMGo'); pg.wait_for_timeout(400)
        mres=S(pg,'S.mres'); chk(mres and (mres.get('mismatch') or mres.get('err')) and pg.query_selector('#dfMNote,.note.warn') is not None,f'{tag}: manual hkl check warns on large angle/ratio mismatch',(mres or {}).get('note'))
    pg.evaluate('()=>{ DIFF_UI._identify(); }'); pg.wait_for_timeout(30)
    busy=S(pg,'S.idBusy'); pg.evaluate('()=>DIFF_UI._cancel()'); pg.wait_for_timeout(200)
    st=S(pg,'{busy:S.idBusy,info:S.idInfo}')
    chk(not st['busy'] and (not busy or (st['info'] or {}).get('cancelled')),f'{tag}: Cancel stops the phase search (was busy {busy})',st)
    pg.click('#dfId'); wait_res(pg); chk(S(pg,'(S.cands||[]).length')>0 or S(pg,'!!S.idInfo'),f'{tag}: search works again after Cancel (worker recreated)',S(pg,'(S.cands||[]).length'))
    pg.screenshot(path=f'{out}/{tag}_si16.png',full_page=True); chk(not hscroll(pg),f'{tag} si16: no page horizontal scroll')
    # 4. HRTEM 비슷한 8-bit 무늬 → FFT
    openf(pg,HR); chk(mode(pg)=='fft',f'{tag}: HRTEM-like 8-bit image opens in FFT',mode(pg))
with sync_playwright() as p:
    b=p.chromium.launch()
    for tag,kw in [('desk',dict(viewport={'width':1400,'height':900})),('phone',dict(viewport={'width':390,'height':844},has_touch=True,is_mobile=True,device_scale_factor=2))]:
        pg=b.new_page(**kw)
        pg.on('pageerror',lambda e: errs.append('pageerror: '+str(e)))
        pg.on('console',lambda m: errs.append('console: '+m.text) if m.type=='error' and 'ERR_' not in m.text and 'gstatic' not in m.text else None)
        pg.on('dialog',lambda d: d.accept())
        run(pg,tag); pg.close()
    chk(not errs,'no console/page errors',errs[:5])
    b.close()
print(sum(1 for r in res if r[0]),'/',len(res),'pass')
sys.exit(0 if all(r[0] for r in res) else 1)
