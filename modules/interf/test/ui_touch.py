"""휴대폰 터치 시험 (Playwright · Chromium) — 공용 제스처 window.INNO_TOUCH 를 쓰는 모든 캔버스.
사용: python3 ui_touch.py <주입된 html> <그림 폴더> [--shots-only]
  모바일 문맥(390×844, hasTouch, isMobile)에서 CDP Input.dispatchTouchEvent 로 손가락 여러 개를 합성한다.
  캔버스마다(EM Image · EDS Maps · EDS Spectrum · Lattice FFT 실공간·FFT · Image Analysis · XPS)
    - 두 손가락 핀치(중심도 조금 옮김) → 확대 비가 손가락 거리 비와 같고, 시작 중심 아래의 내용이 끝 중심 아래에 있다(3 px 안)
    - 한 손가락 끌기 → 시작점 아래의 내용이 끝점 아래로 (3 px 안)
    - 두 번 탭 → 처음 보기(맞춤)
  그리고 Image Analysis: 손가락 사각형 ROI 좌표 = 기대값, 그리는 중 둘째 손가락 → 그리던 ROI 취소, 다각형 탭·Finish,
         Point 도구 탭 = 점 · 길게 누르기 = 지우기, Hand 탭 → 상태 줄에 화소 값
  XPS: 'Add component' 뒤 탭 → 그 BE 에 성분, 탭 → 값 풍선 / EM: 탭 → 화소 값, 손가락 거리 측정
  Lattice: 확대한 FFT 를 탭 → 그 자리 스폿 측정
  쪽: 가로 스크롤 없음(390·360 px), 캔버스 밖 손가락 끌기 = 쪽 스크롤, 주요 단추·탭 높이 ≥ 40 px, 콘솔 오류 0.
  --shots-only: 시험 없이 탭별 휴대폰 화면만 저장 (고치기 전·후 비교용)
선택자는 id / data-* 만 쓴다."""
import sys, os, json, glob
from playwright.sync_api import sync_playwright
html, out = sys.argv[1], sys.argv[2]; SHOTS = '--shots-only' in sys.argv; os.makedirs(out, exist_ok=True)
H = os.path.dirname(os.path.abspath(__file__)); M = os.path.join(H, '..', '..')
TIF = os.path.join(M, 'image', 'test', 'io', 't16.tif'); IAF = os.path.join(M, 'image', 'test', 'io', 'tfei.tif')
XPSF = sorted(glob.glob(os.path.join(M, 'xps', 'test', 'data', 'bct600', '*.xls')))
res = []; errs = []
def chk(c, name, info=''):
    res.append((bool(c), name, str(info))); print(('PASS' if c else 'FAIL'), name, info, flush=True)
SYN = '''(()=>{
  const rec=curRec(); const orig=rec.h5.get.bind(rec.h5);
  const n=2048, disp=10, off=0; const y=new Float64Array(n);
  const peaks=[[0.677,900],[1.04,300],[1.74,200],[8.04,700],[1.5,150],[9.713,600],[7.18,250],[6.0,150]];
  for(let i=0;i<n;i++){ const E=(off+disp*i)/1000; let v=40*Math.exp(-E/3)+2; for(const [pe,a] of peaks){ const sg=(0.055+0.002*pe)/2.355; v+=a*Math.exp(-0.5*((E-pe)/sg)**2); } y[i]=Math.round(v*(0.9+0.2*((i*2654435761>>>0)%1000)/1000)); }
  const store={"/syn/spectrum/0/Data":{shape:[n,1],dtype:"<Q",value:y}};
  rec.h5={get(p){ const d=store[p]; if(d) return {shape:d.shape,dtype:d.dtype,value:d.value,slice(){return d.value;}}; return orig(p); }, close(){}};
  const md={Detectors:{"Detector-0":{DetectorType:"AnalyticalDetector",DetectorName:"Synthetic EDS",Dispersion:String(disp),OffsetEnergy:String(off),BeginEnergy:"150",FwhmMnKa:"128"}}};
  rec.spectra.push({kind:"spectrum",uid:"/syn/spectrum/0",base:"/syn/spectrum/0",shape:[n,1],dtype:"<Q",md,px:null,det:"Synthetic EDS",role:"Summed"});
  const im0=rec.images[0];
  for(const t of ["Lu","Yb","F"]) rec.images.push(Object.assign({},im0,{uid:"m_"+t,isMap:true,title:t}));
  state.curSpec=rec.spectra.length-1; state.curFile=state.files.indexOf(rec); renderView();
})()'''
# 보기 ↔ 화면 좌표 (시험 쪽에서 앱과 따로 계산 — 앱의 변환을 그대로 믿지 않는다)
PROBES = r'''(()=>{
const R=s=>document.querySelector(s).getBoundingClientRect();
const css=(sel,get)=>({fwd(x,y){ const c=document.querySelector(sel), r=c.getBoundingClientRect(); return [(x-r.left)/r.width*c.width,(y-r.top)/r.height*c.height]; },
  inv(u,w){ const c=document.querySelector(sel), r=c.getBoundingClientRect(); return [r.left+u/c.width*r.width, r.top+w/c.height*r.height]; }, zoom:get, el:sel});
window.TP={
 em:{el:'#cv', vr(){ const S=state.img, c=loadImageData(); const z=S.zoom, sw=c.w/z, sh=c.h/z; let sx=S.cx*c.w-sw/2, sy=S.cy*c.h-sh/2; sx=Math.max(0,Math.min(c.w-sw,sx)); sy=Math.max(0,Math.min(c.h-sh,sy)); return {sx,sy,sw,sh}; },
   fwd(x,y){ const r=R('#cv'), v=this.vr(); return [v.sx+(x-r.left)/r.width*v.sw, v.sy+(y-r.top)/r.height*v.sh]; },
   inv(u,w){ const r=R('#cv'), v=this.vr(); return [r.left+(u-v.sx)/v.sw*r.width, r.top+(w-v.sy)/v.sh*r.height]; }, zoom(){ return state.img.zoom; }},
 maps:Object.assign(css('#mcv',()=>MAPFIT.view.s),{el:'#mstage'}),
 latR:css('#latReal',()=>LATZ.real.s), latF:css('#latFFT',()=>LATZ.fft.s),
 ia:{el:'#imOverlay', v(){ return IMG_UI._state().ui.views[IMG_UI._cur().id]; },
   fwd(x,y){ const r=R('#imOverlay'), v=this.v(); return [v.x+(x-r.left)/v.zoom, v.y+(y-r.top)/v.zoom]; },
   inv(u,w){ const r=R('#imOverlay'), v=this.v(); return [r.left+(u-v.x)*v.zoom, r.top+(w-v.y)*v.zoom]; }, zoom(){ return this.v().zoom; }},
 xps:{el:'#xpsSvg', m(){ return document.querySelector('#xpsSvg').getScreenCTM(); },
   fwd(x,y){ const G=XPS_UI._geo(), m=this.m(); const q=new DOMPoint(x,y).matrixTransform(m.inverse()); return [G.inv(q.x),y]; },
   inv(be,y){ const G=XPS_UI._geo(), m=this.m(); return [new DOMPoint(G.X(be),0).matrixTransform(m).x,y]; }, zoom(){ const v=XPS_UI._geo().v; return 1/(v.x1-v.x0); }},
 sp:{el:'#spPlot', rng(){ if(SPVIEW.x0!=null) return [Math.max(window.__sp[0],SPVIEW.x0),Math.min(window.__sp[1],SPVIEW.x1)]; return window.__sp; },
   fwd(x,y){ const r=R('#spPlot svg'), [a,b]=this.rng(); return [a+((x-r.left)*860/r.width-52)/796*(b-a),y]; },
   inv(E,y){ const r=R('#spPlot svg'), [a,b]=this.rng(); return [r.left+(52+(E-a)/(b-a)*796)*r.width/860,y]; }, zoom(){ const [a,b]=this.rng(); return 1/(b-a); }}
};})()'''
with sync_playwright() as p:
    b = p.chromium.launch()
    ctx = b.new_context(viewport={'width': 390, 'height': 844}, has_touch=True, is_mobile=True, device_scale_factor=2)
    pg = ctx.new_page()
    pg.on('pageerror', lambda e: errs.append('pageerror: ' + str(e)))
    pg.on('console', lambda m: errs.append('console: ' + m.text) if m.type == 'error' and 'ERR_' not in m.text else None)
    pg.on('dialog', lambda d: d.accept())
    cdp = ctx.new_cdp_session(pg)
    def T(typ, pts):
        cdp.send('Input.dispatchTouchEvent', {'type': typ, 'touchPoints': [{'x': float(x), 'y': float(y), 'id': i, 'radiusX': 1, 'radiusY': 1, 'force': 1} for i, x, y in pts]})
    def pinch(cx, cy, d0, d1, shift=(0, 0), steps=10):
        P = lambda t: (cx + shift[0] * t, cy + shift[1] * t, d0 + (d1 - d0) * t)
        x, y, d = P(0); T('touchStart', [(0, x - d / 2, y)]); pg.wait_for_timeout(20)
        T('touchStart', [(0, x - d / 2, y), (1, x + d / 2, y)]); pg.wait_for_timeout(20)
        for k in range(1, steps + 1):
            x, y, d = P(k / steps); T('touchMove', [(0, x - d / 2, y), (1, x + d / 2, y)]); pg.wait_for_timeout(17)
        T('touchEnd', [(0, x - d / 2, y)]); pg.wait_for_timeout(20); T('touchEnd', []); pg.wait_for_timeout(120)
        return x, y
    def drag(x, y, dx, dy, steps=10, hold=0):
        T('touchStart', [(0, x, y)]); pg.wait_for_timeout(20 + hold)
        for k in range(1, steps + 1): T('touchMove', [(0, x + dx * k / steps, y + dy * k / steps)]); pg.wait_for_timeout(17)
        T('touchEnd', []); pg.wait_for_timeout(120)
    def tap(x, y, wait=60):
        T('touchStart', [(0, x, y)]); pg.wait_for_timeout(30); T('touchEnd', []); pg.wait_for_timeout(wait)
    def dtap(x, y):
        tap(x, y, 70); tap(x, y, 500)
    def longpress(x, y, ms=700):
        T('touchStart', [(0, x, y)]); pg.wait_for_timeout(ms); T('touchEnd', []); pg.wait_for_timeout(150)
    def box(sel):
        pg.evaluate("s=>document.querySelector(s).scrollIntoView({block:'center'})", sel); pg.wait_for_timeout(150)
        return pg.evaluate("s=>{ const r=document.querySelector(s).getBoundingClientRect(); return [r.left,r.top,r.width,r.height]; }", sel)
    def hscroll(): return pg.evaluate('document.documentElement.scrollWidth>window.innerWidth+1')
    def go(v):
        if v == 'meta':  # EM Metadata 는 탭이 아니라 EM Image 안의 보기 (10-19)
            pg.click('#nav button[data-v="image"]'); pg.wait_for_timeout(500); pg.click('.emsub button[data-sv="meta"]'); pg.wait_for_timeout(700); return
        pg.click(f'#nav button[data-v="{v}"]'); pg.wait_for_timeout(700)
    def gestures(key, label, cx, cy, d0=70, d1=180, shift=(14, 8), pan=(-50, -36), only_x=False, tol=3.0):
        """핀치·끌기·두 번 탭 — key = window.TP 의 항목"""
        z0 = pg.evaluate(f'TP.{key}.zoom()'); c0 = pg.evaluate(f'TP.{key}.fwd({cx},{cy})')
        ex, ey = pinch(cx, cy, d0, d1, shift)
        z1 = pg.evaluate(f'TP.{key}.zoom()'); s = pg.evaluate(f'TP.{key}.inv({c0[0]},{c0[1]})')
        err = abs(s[0] - ex) if only_x else max(abs(s[0] - ex), abs(s[1] - ey))
        chk(abs(z1 / z0 - d1 / d0) < 0.05 * d1 / d0, f'{label}: pinch zoom ×{d1/d0:.2f}', f'{z0:.4g} → {z1:.4g} (×{z1/z0:.3f})')
        chk(err <= tol, f'{label}: pinch keeps the content under the pinch centre', f'{err:.2f} px')
        c1 = pg.evaluate(f'TP.{key}.fwd({ex},{ey})'); drag(ex, ey, pan[0], pan[1])
        s2 = pg.evaluate(f'TP.{key}.inv({c1[0]},{c1[1]})'); err2 = abs(s2[0] - ex - pan[0]) if only_x else max(abs(s2[0] - ex - pan[0]), abs(s2[1] - ey - pan[1]))
        chk(err2 <= tol and abs(pg.evaluate(f'TP.{key}.zoom()') / z1 - 1) < 1e-6, f'{label}: one-finger pan moves the view with the finger', f'{err2:.2f} px')
        return z0, (ex + pan[0], ey + pan[1])
    def dtap_reset(key, label, x, y, z0):
        dtap(x, y); z = pg.evaluate(f'TP.{key}.zoom()'); chk(abs(z / z0 - 1) < 1e-3, f'{label}: double-tap resets the view', f'{z:.4g} vs {z0:.4g}')
    def sizes(tab):
        """휴대폰에서 주요 단추·도구·탭 높이 (40 px 미만이면 목록)"""
        small = pg.evaluate('''()=>{ const sel=['#nav button','.btn','.seg button','.imx .im-toolbar .tb','.xps .xps-tabs button','.imx .im-tabs button'];
          const out=[]; document.querySelectorAll(sel.join(',')).forEach(e=>{ const r=e.getBoundingClientRect(); if(!r.width||!r.height) return; const cs=getComputedStyle(e); if(cs.visibility==='hidden'||cs.display==='none') return;
            if(e.closest('[hidden]')||e.closest('.side.collapsed')) return; if(r.height<39.5) out.push((e.id||e.dataset.v||e.dataset.t||e.className)+':'+Math.round(r.height)+' '+(e.textContent||'').trim().slice(0,16)); }); return out; }''')
        chk(not small, f'{tab}: main buttons, tool palette and tabs are ≥ 40 px tall on the phone', small[:8])

    pg.goto('file://' + os.path.abspath(html)); pg.wait_for_timeout(1500)
    chk(pg.evaluate("matchMedia('(pointer: coarse)').matches"), 'mobile context reports pointer: coarse')
    if not SHOTS: chk(pg.evaluate('!!(window.INNO_TOUCH&&INNO_TOUCH.attach&&INNO_TOUCH.cssZoom&&INNO_TOUCH.affine)'), 'window.INNO_TOUCH API present')
    # ---------- EM 파일 ----------
    pg.set_input_files('#fin', TIF); pg.wait_for_timeout(1500)
    go('image'); pg.screenshot(path=f'{out}/p_image.png')
    if not SHOTS:
        pg.evaluate(PROBES)
        x, y, w, h = box('#cv'); cx, cy = x + w * 0.45, y + h * 0.5
        sizes('EM Image')
        z0, end = gestures('em', 'EM Image', cx, cy)
        pg.screenshot(path=f'{out}/t_image_zoomed.png')
        tap(cx, cy, 500); ro = pg.inner_text('#hoverRO'); chk('=' in ro, 'EM Image: tap shows the pixel value (hover readout)', ro)
        dtap_reset('em', 'EM Image', cx, cy, z0)
        # 거리 도구: 한 손가락으로 선을 그림 → 측정 카드
        pg.click('#toolSeg button[data-t="line"]'); pg.wait_for_timeout(200); x, y, w, h = box('#cv')
        a = pg.evaluate(f'TP.em.fwd({x+w*0.3},{y+h*0.4})'); bb = pg.evaluate(f'TP.em.fwd({x+w*0.7},{y+h*0.6})')
        drag(x + w * 0.3, y + h * 0.4, w * 0.4, h * 0.2); pg.wait_for_timeout(300)
        p01 = pg.evaluate('[state.img.p0,state.img.p1]')
        chk(p01[0] and abs(p01[0]['x'] - a[0]) < 1.5 and abs(p01[1]['x'] - bb[0]) < 1.5 and abs(p01[1]['y'] - bb[1]) < 1.5 and pg.query_selector('#profPlot') is not None, 'EM Image: finger draws the distance/profile line', p01)
        # 그리는 중 둘째 손가락 → 핀치로 넘어가고 그리던 선은 되돌림
        T('touchStart', [(0, x + w * 0.2, y + h * 0.2)]); pg.wait_for_timeout(30); T('touchMove', [(0, x + w * 0.25, y + h * 0.25)]); pg.wait_for_timeout(30)
        T('touchStart', [(0, x + w * 0.25, y + h * 0.25), (1, x + w * 0.6, y + h * 0.6)]); pg.wait_for_timeout(30)
        T('touchMove', [(0, x + w * 0.2, y + h * 0.2), (1, x + w * 0.7, y + h * 0.7)]); pg.wait_for_timeout(30); T('touchEnd', [(0, x + w * 0.2, y + h * 0.2)]); T('touchEnd', []); pg.wait_for_timeout(150)
        p02 = pg.evaluate('[state.img.p0,state.img.p1,state.img.zoom]')
        chk(p02[0] == p01[0] and p02[1] == p01[1] and p02[2] > 1.2, 'EM Image: second finger cancels the line being drawn and pinches', p02)
        pg.click('#btnReset'); pg.click('#toolSeg button[data-t="none"]'); pg.wait_for_timeout(200)
        # 쪽 스크롤: 캔버스 밖(타일)에서 끌면 쪽이 움직이고, 캔버스 위에서는 안 움직인다
        pg.evaluate('window.scrollTo(0,0)'); pg.wait_for_timeout(200)
        tx, ty, tw, th = pg.evaluate("(()=>{ const r=document.querySelector('.tiles').getBoundingClientRect(); return [r.left,r.top,r.width,r.height]; })()")
        T('touchStart', [(0, tx + tw / 2, ty + th / 2 + 40)]);
        for k in range(1, 11): T('touchMove', [(0, tx + tw / 2, ty + th / 2 + 40 - 25 * k)]); pg.wait_for_timeout(16)
        T('touchEnd', []); pg.wait_for_timeout(600); sy = pg.evaluate('scrollY')
        chk(sy > 50, 'page still scrolls with a finger outside the canvases', f'scrollY {sy}')
        x, y, w, h = box('#cv'); sy0 = pg.evaluate('scrollY'); drag(x + w / 2, y + h / 2, 0, -120); sy1 = pg.evaluate('scrollY')
        chk(abs(sy1 - sy0) < 2, 'finger on the image canvas does not scroll the page', f'{sy0} → {sy1}')
    chk(not hscroll(), 'EM Image: no page horizontal scroll')
    # ---------- Lattice ----------
    go('lattice')
    if pg.query_selector('#latPxIn'):
        pg.fill('#latPxIn', '50'); pg.click('#latPxGo'); pg.wait_for_timeout(1200)
    pg.screenshot(path=f'{out}/p_lattice.png')
    if not SHOTS:
        pg.evaluate(PROBES); sizes('Lattice FFT')
        x, y, w, h = box('#latReal'); z0, end = gestures('latR', 'Lattice real space', x + w * 0.5, y + h * 0.5, d0=60, d1=150, pan=(30, 20))
        dtap_reset('latR', 'Lattice real space', x + w * 0.5, y + h * 0.5, z0)
        x, y, w, h = box('#latFFT'); z0, end = gestures('latF', 'Lattice FFT', x + w * 0.5, y + h * 0.5, d0=60, d1=150, pan=(30, 20))
        pg.screenshot(path=f'{out}/t_lattice_zoomed.png')
        # 스폿 고르기: 자동 탐색으로 스폿 자리를 알아낸 뒤 지우고, 확대한 FFT 위의 그 자리를 손가락으로 탭
        pg.click('#latAuto'); pg.wait_for_timeout(800)
        spot = pg.evaluate('(()=>{ const L=latState(); return L.spots.length?[L.spots[0].m.a.x,L.spots[0].m.a.y,L.spots[0].m.d]:null; })()')
        if spot:
            pg.click('#latClr'); pg.wait_for_timeout(500); x, y, w, h = box('#latFFT')
            dtap(x + w / 2, y + h / 2); x, y, w, h = box('#latFFT')       # 맞춤 → 스폿 자리를 중심으로 핀치 → 그 자리를 탭
            spxy = lambda: pg.evaluate(f"(()=>{{ const c=document.querySelector('#latFFT'); return TP.latF.inv({spot[0]}*c._sc,{spot[1]}*c._sc); }})()")
            q = spxy(); pinch(q[0], q[1], 50, 125); zf = pg.evaluate('LATZ.fft.s'); q = spxy()
            tap(q[0] + 1.5, q[1] - 1.0, 900)
            got = pg.evaluate('(()=>{ const L=latState(); return L.spots.map(s=>[s.m.a.x,s.m.a.y,s.m.d,s.pick]); })()')
            chk(len(got) == 1 and abs(got[0][2] - spot[2]) / spot[2] < 0.01 and got[0][3] == 'click', 'Lattice FFT: finger tap on the (zoomed) FFT measures that spot', {'auto': spot, 'tap': got, 'zoom': zf})
            chk(abs(pg.evaluate('LATZ.fft.s') - zf) < 1e-9 or zf == 1, 'Lattice FFT: zoom kept after the spot re-render', pg.evaluate('LATZ.fft.s'))
        else:
            chk(False, 'Lattice FFT: auto-find gave a spot to tap', 'no spot')
        x, y, w, h = box('#latFFT'); dtap_reset('latF', 'Lattice FFT', x + w * 0.5, y + h * 0.5, 1.0)
    chk(not hscroll(), 'Lattice FFT: no page horizontal scroll')
    # ---------- EDS Spectrum · Maps (합성 자료) ----------
    pg.evaluate(SYN); pg.wait_for_timeout(500)
    go('spectrum'); pg.screenshot(path=f'{out}/p_spectrum.png')
    if not SHOTS:
        pg.evaluate(PROBES); sizes('EDS Spectrum')
        pg.evaluate('(()=>{ const S=loadSpectrum(), em=state.spec.eMax; let a=null,b=null; for(let i=S.i0|0;i<S.nch;i++){ if(S.E[i]>em) break; if(a==null) a=S.E[i]; b=S.E[i]; } window.__sp=[a,b]; })()')
        x, y, w, h = box('#spPlot svg'); cx = x + w * 0.4; cy = y + h * 0.5
        z0, end = gestures('sp', 'EDS Spectrum', cx, cy, shift=(10, 0), pan=(40, 0), only_x=True, tol=2.0)
        pg.screenshot(path=f'{out}/t_spectrum_zoomed.png')
        dtap_reset('sp', 'EDS Spectrum', cx, cy, z0)
    chk(not hscroll(), 'EDS Spectrum: no page horizontal scroll')
    go('maps'); pg.wait_for_timeout(500); pg.screenshot(path=f'{out}/p_maps.png')
    if not SHOTS:
        pg.evaluate(PROBES); sizes('EDS Maps')
        x, y, w, h = box('#mcv'); z0, end = gestures('maps', 'EDS Maps', x + w * 0.5, y + h * 0.5, pan=(-30, -20))
        sb = pg.evaluate("document.querySelector('#msbar').style.display"); pg.screenshot(path=f'{out}/t_maps_zoomed.png')
        dtap_reset('maps', 'EDS Maps', x + w * 0.5, y + h * 0.5, z0)
    chk(not hscroll(), 'EDS Maps: no page horizontal scroll')
    # ---------- Image Analysis ----------
    go('imgx'); pg.wait_for_selector('.imx'); pg.set_input_files('#imFin', [IAF]); pg.wait_for_function('IMG_UI._state().images.length>=1', timeout=20000); pg.wait_for_timeout(600)
    pg.screenshot(path=f'{out}/p_imgx.png')
    if not SHOTS:
        pg.evaluate(PROBES); sizes('Image Analysis')
        pg.click('#imTools [data-t="hand"]'); pg.wait_for_timeout(150)
        x, y, w, h = box('#imOverlay'); cx, cy = x + w * 0.5, y + h * 0.5
        z0, end = gestures('ia', 'Image Analysis', cx, cy)
        tap(cx, cy, 500); st = pg.inner_text('#imStatus'); chk('value' in st, 'Image Analysis: tap shows the pixel value in the status bar', st)
        longpress(cx + 20, cy + 10); st2 = pg.inner_text('#imStatus'); chk('value' in st2 and st2 != st, 'Image Analysis: long-press reads the value under the finger', st2)
        pg.screenshot(path=f'{out}/t_imgx_zoomed.png')
        dtap_reset('ia', 'Image Analysis', cx, cy, z0)
        # 사각형 ROI (기본 도구) — 한 손가락 끌기
        pg.click('#imTools [data-t="rect"]'); pg.wait_for_timeout(150); x, y, w, h = box('#imOverlay')
        A = (x + w * 0.3, y + h * 0.35); B = (x + w * 0.62, y + h * 0.7)
        ia = pg.evaluate(f'TP.ia.fwd({A[0]},{A[1]})'); ib = pg.evaluate(f'TP.ia.fwd({B[0]},{B[1]})')
        drag(A[0], A[1], B[0] - A[0], B[1] - A[1]); pg.wait_for_timeout(200)
        roi = pg.evaluate('(()=>{ const S=IMG_UI._state(), r=S.roiSel[IMG_UI._cur().id]; return r?{type:r.type,pts:r.pts}:null; })()')
        exp = [[int(ia[0] // 1), int(ia[1] // 1)], [int(ib[0] // 1), int(ib[1] // 1)]]
        chk(roi and roi['type'] == 'rect' and roi['pts'] == exp, 'Image Analysis: rectangle ROI drawn by finger has the expected coordinates', {'got': roi, 'exp': exp})
        # 두 손가락은 언제나 핀치: 그리던 새 사각형은 취소되고 앞의 ROI 가 남는다
        zb = pg.evaluate('TP.ia.zoom()'); C = (x + w * 0.2, y + h * 0.2)
        T('touchStart', [(0, C[0], C[1])]); pg.wait_for_timeout(30); T('touchMove', [(0, C[0] + 20, C[1] + 25)]); pg.wait_for_timeout(30)
        T('touchStart', [(0, C[0] + 20, C[1] + 25), (1, C[0] + 120, C[1] + 125)]); pg.wait_for_timeout(30)
        for k in range(1, 6): T('touchMove', [(0, C[0] + 20 - 8 * k, C[1] + 25 - 8 * k), (1, C[0] + 120 + 8 * k, C[1] + 125 + 8 * k)]); pg.wait_for_timeout(17)
        T('touchEnd', [(0, C[0] - 20, C[1] - 15)]); T('touchEnd', []); pg.wait_for_timeout(150)
        roi2 = pg.evaluate('(()=>{ const S=IMG_UI._state(), r=S.roiSel[IMG_UI._cur().id]; return r?{type:r.type,pts:r.pts}:null; })()')
        chk(roi2 == roi and pg.evaluate('TP.ia.zoom()') > zb * 1.3, 'Image Analysis: two fingers pinch even with a drawing tool (the half-drawn ROI is cancelled)', {'roi': roi2, 'zoom': [zb, pg.evaluate('TP.ia.zoom()')]})
        # 손잡이: 손가락이면 크게(허용 14 px) — 모서리에서 10 px 떨어진 곳을 끌어도 크기 조절
        pg.evaluate('(()=>{ const im=IMG_UI._cur(); const S=IMG_UI._state(); delete S.ui.views[im.id]; IMG_UI.rerender(); })()'); pg.wait_for_timeout(300); x, y, w, h = box('#imOverlay')
        roi = pg.evaluate('(()=>{ const r=IMG_UI._state().roiSel[IMG_UI._cur().id]; return r.pts; })()'); br = pg.evaluate(f'TP.ia.inv({roi[1][0]},{roi[1][1]})')
        drag(br[0] + 9, br[1] + 9, 30, 20); roi3 = pg.evaluate('(()=>{ const r=IMG_UI._state().roiSel[IMG_UI._cur().id]; return r.pts; })()')
        chk(roi3[0] == roi[0] and roi3[1][0] > roi[1][0] and roi3[1][1] > roi[1][1], 'Image Analysis: enlarged ROI handle can be dragged with a finger', {'before': roi, 'after': roi3})
        # 다각형: 탭 세 번 → Finish
        pg.click('#imTools [data-t="poly"]'); pg.wait_for_timeout(150); x, y, w, h = box('#imOverlay')
        pts = [(x + w * 0.3, y + h * 0.3), (x + w * 0.7, y + h * 0.35), (x + w * 0.5, y + h * 0.7)]
        for q in pts: tap(q[0], q[1], 400)
        vis = pg.evaluate("!document.querySelector('#imPolyBar').hidden"); pg.screenshot(path=f'{out}/t_imgx_poly.png')
        pg.click('#imPolyDone'); pg.wait_for_timeout(200)
        poly = pg.evaluate('(()=>{ const r=IMG_UI._state().roiSel[IMG_UI._cur().id]; return r?{type:r.type,n:r.pts.length}:null; })()')
        chk(vis and poly == {'type': 'poly', 'n': 3}, 'Image Analysis: polygon by taps, finished with the Finish button', {'bar': vis, 'roi': poly})
        # Point 도구: 탭 = 점, 길게 누르기 = 가까운 점 지우기
        pg.click('#imTools [data-t="point"]'); pg.wait_for_timeout(150); x, y, w, h = box('#imOverlay')
        tap(x + w * 0.4, y + h * 0.4, 120); tap(x + w * 0.6, y + h * 0.5, 120)
        n1 = pg.evaluate('(()=>{ const r=IMG_UI._state().roiSel[IMG_UI._cur().id]; return r&&r.type==="point"?r.pts.length:0; })()')
        longpress(x + w * 0.6 + 4, y + h * 0.5 + 3)
        n2 = pg.evaluate('(()=>{ const r=IMG_UI._state().roiSel[IMG_UI._cur().id]; return r&&r.type==="point"?r.pts.length:0; })()')
        chk(n1 == 2 and n2 == 1, 'Image Analysis: Point tool — tap adds, long-press removes (the touch form of Alt+click)', [n1, n2])
        pg.click('#imTools [data-t="rect"]')
    chk(not hscroll(), 'Image Analysis: no page horizontal scroll')
    # ---------- XPS ----------
    go('xps'); pg.wait_for_selector('.xps'); pg.set_input_files('#xpsFin', XPSF)
    pg.wait_for_function('XPS_UI._state().specs.length>=6', timeout=30000); pg.wait_for_timeout(500)
    pg.evaluate("(()=>{ const S=XPS_UI._state(); const s=S.specs.find(x=>/C\\s*1s/i.test(x.region))||S.specs[0]; S.sel=s.uid; XPS_UI.rerender(); })()"); pg.wait_for_timeout(500)
    pg.screenshot(path=f'{out}/p_xps.png')
    if not SHOTS:
        pg.evaluate(PROBES); sizes('XPS')
        x, y, w, h = box('#xpsSvg'); cx, cy = x + w * 0.5, y + h * 0.35
        z0, end = gestures('xps', 'XPS', cx, cy, shift=(10, 0), pan=(40, 0), only_x=True, tol=2.0)
        tap(cx, cy, 500); tipv = pg.evaluate("(()=>{ const t=document.querySelector('#xpsTip'); return t&&!t.hidden?t.textContent:null; })()")
        chk(tipv and 'BE' in tipv, 'XPS: tap shows the value tip (hover info)', tipv)
        pg.screenshot(path=f'{out}/t_xps_zoomed.png')
        dtap_reset('xps', 'XPS', cx, cy, z0)
        # 성분 놓기: Add component → 그래프를 탭한 자리의 BE
        n0 = pg.evaluate('XPS_UI._pr(XPS_UI._cur()).comps.length'); pg.click('#xpsAdd'); pg.wait_for_timeout(200)
        x, y, w, h = box('#xpsSvg'); tx = x + w * 0.55; ty = y + h * 0.3
        exp = pg.evaluate(f'TP.xps.fwd({tx},{ty})[0]-XPS_UI._geo().sh'); tap(tx, ty, 400)
        cs = pg.evaluate('XPS_UI._pr(XPS_UI._cur()).comps.map(c=>c.be)')
        span = pg.evaluate('(()=>{ const v=XPS_UI._geo().v; return v.x1-v.x0; })()')
        chk(len(cs) == n0 + 1 and abs(cs[-1] - exp) < span / 300, 'XPS: tap places a component at the tapped BE', {'be': cs[-1] if cs else None, 'expected': exp})
        chk(not pg.evaluate('XPS_UI._state().ui.addMode'), 'XPS: place mode ends after one tap (as with the mouse)')
        # 성분 표지를 한 손가락으로 끌기 (앱의 끌기 처리 그대로) → BE 가 손가락만큼
        pg.wait_for_timeout(300); x, y, w, h = box('#xpsSvg')
        mk = pg.evaluate("(()=>{ const p=XPS_UI._pr(XPS_UI._cur()), c=p.comps[p.comps.length-1]; const e=document.querySelector('#xpsSvg [data-mk=\"'+c.id+'\"]'); const r=e.getBoundingClientRect(); return [r.left+r.width/2,r.top+r.height/2,c.be]; })()")
        bx = pg.evaluate(f'TP.xps.fwd({mk[0]+30},{mk[1]})[0]-XPS_UI._geo().sh') - pg.evaluate(f'TP.xps.fwd({mk[0]},{mk[1]})[0]-XPS_UI._geo().sh')
        drag(mk[0], mk[1], 30, 0); be2 = pg.evaluate('(()=>{ const p=XPS_UI._pr(XPS_UI._cur()); return p.comps[p.comps.length-1].be; })()')
        chk(abs((be2 - mk[2]) - bx) < span / 200, 'XPS: component marker dragged with one finger', {'before': mk[2], 'after': be2, 'expected shift': bx})
    chk(not hscroll(), 'XPS: no page horizontal scroll')
    # ---------- 나머지 탭: 단추·탭 크기, 가로 넘침 ----------
    for v, nm in [('interf', 'EDS Interference'), ('meta', 'EM Metadata'), ('energy', 'Energy Diagram'), ('dre', 'DRE Kinetics')]:
        go(v); pg.screenshot(path=f'{out}/p_{v}.png')
        if not SHOTS: sizes(nm)
        chk(not hscroll(), f'{nm}: no page horizontal scroll')
        if v == 'energy' and not SHOTS:
            pg.fill('#edF', 'NaYF4:20%Yb,2%Er'); pg.click('#edApply'); pg.wait_for_timeout(600)
            # 준위를 손가락으로 탭 → 풍선(hover 대신), 도식 밖 탭 → 닫힘
            hit = pg.evaluate("(()=>{ const e=[...document.querySelectorAll('.lvhit.tipable')].find(x=>{ const r=x.getBoundingClientRect(); return r.width>4&&r.height>2; }); if(!e) return null; e.scrollIntoView({block:'center'}); return 1; })()")
            if hit:
                pg.wait_for_timeout(200); hit = pg.evaluate("(()=>{ const e=[...document.querySelectorAll('.lvhit.tipable')].find(x=>{ const r=x.getBoundingClientRect(); return r.width>4&&r.height>2&&r.top>60&&r.bottom<window.innerHeight-80; }); const r=e.getBoundingClientRect(); return [r.left+Math.min(r.width/2,20),r.top+r.height/2]; })()")
                tap(hit[0], hit[1], 300); on = pg.evaluate("!!(document.querySelector('#edTip')&&document.querySelector('#edTip').classList.contains('on'))")
                txt = pg.evaluate("(document.querySelector('#edTip')||{}).textContent||''")
                chk(on and len(txt) > 3, 'Energy Diagram: finger tap on a level opens its tooltip', txt[:60])
                tap(12, 12, 200); off = pg.evaluate("!(document.querySelector('#edTip')&&document.querySelector('#edTip').classList.contains('on'))")
                chk(off, 'Energy Diagram: tap outside closes the tooltip')
            else: chk(False, 'Energy Diagram: a level is visible to tap')
    # ---------- 360 px 폭 ----------
    pg.set_viewport_size({'width': 360, 'height': 780}); pg.wait_for_timeout(300)
    over = []
    for v in ['xps', 'image', 'imgx', 'maps', 'spectrum', 'lattice', 'interf', 'meta', 'energy', 'dre']:
        go(v)
        if hscroll(): over.append(v)
    pg.screenshot(path=f'{out}/p_360_last.png')
    chk(not over, 'no page horizontal scroll on any tab at 360 px', over)
    pg.wait_for_timeout(300)
    chk(not errs, 'no console/page errors', errs[:5])
    b.close()
n = sum(1 for r in res if r[0]); print(f'{n} / {len(res)} pass'); json.dump(res, open(os.path.join(out, 'result.json'), 'w'), ensure_ascii=False, indent=1)
sys.exit(0 if n == len(res) else 1)
