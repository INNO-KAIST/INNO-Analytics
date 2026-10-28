"""v1.8 Series · ARXPS · Imaging · Processing — 브라우저(Chromium) 실측, 항목마다 판정
사용: python3 test/ui_xps_multi.py <주입된 html> <그림 폴더>
자료: test/data/synthetic/ (make_vamas.py 합성 — 측정 자료 아님) + 이 시험이 만드는 Fermi edge·영상 텍스트
1 Series: SDP VAMAS → 식각 시간 좌표 → Shirley 넓이 → at % (truth ±3 at.%) · PCA·요인 수 · MCR-ALS · 기준 NNLS · 정보 깊이 보정 · CSV
2 여러 열 텍스트 → series 12 단계 · PCA · MCR
3 ARXPS: Angle 변수 → 0° 피팅 성분을 기준으로 NNLS → 덮개층 d ≈ 1.8 nm · Hill · 다층 · RDP · MaxEnt(Worker)
4 Imaging: MAP 8×6 → 넓이 지도(truth) · ROI → 스펙트럼 → 피팅 → NNLS 상태 지도(truth) · PCA·재구성 · binning · 텍스트 stack(nx 바꾸기) · PNG/CSV
5 Process: 디컨볼루션(해상도 이득·잡음 띠) · Fermi edge(E_F·분해능) → PSF · 다시 표본화 · 미분
6 콘솔 오류 없음 · 2 s 넘는 주 스레드 작업 없음 · 390 px 가로 넘침 없음"""
import sys, os, json, math, random, tempfile
from playwright.sync_api import sync_playwright

html, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'synthetic')
T_dp = json.load(open(os.path.join(D, 'depth_profile_truth.json')))
T_ar = json.load(open(os.path.join(D, 'arxps_truth.json')))
T_map = json.load(open(os.path.join(D, 'map_truth.json')))
errs, res, bad = [], {}, []


def check(name, cond, info=''):
    res[name] = {'ok': bool(cond), 'info': info}
    if not cond:
        bad.append(name)
    print(('  ok   ' if cond else '  FAIL ') + name + ('  (' + str(info) + ')' if info != '' else ''), flush=True)


tmp = tempfile.mkdtemp()
random.seed(5)
# Fermi edge (Au 가전자대 흉내): 선형 DOS × Fermi–Dirac(300 K) ⊗ Gauss(FWHM 0.40 eV), E_F = 0 (BE), Poisson
kT = 8.617333e-5 * 300; FW = 0.40; sg = FW / (2 * math.sqrt(2 * math.log(2)))
xs = [-1.5 + 0.02 * i for i in range(176)]
def fd(e): return 1 / (1 + math.exp(-e / kT))       # BE 축: 점유(BE > 0)
def model(e):
    s = 0; w = 0
    for k in range(-60, 61):
        t = k * sg / 15; g = math.exp(-0.5 * (t / sg) ** 2); s += g * (2000 * (1 + 0.2 * (e - t)) * fd(e - t)); w += g
    return 300 + s / w
with open(os.path.join(tmp, '[Au foil] VB.txt'), 'w') as fo:
    fo.write('Binding Energy (eV)\tCounts\n')
    for x in xs:
        lam = model(x); fo.write(f'{x:.3f}\t{max(0, round(lam + math.sqrt(lam) * random.gauss(0, 1)))}\n')

with sync_playwright() as p:
    b = p.chromium.launch()
    ctx = b.new_context(viewport={'width': 1500, 'height': 1000}, accept_downloads=True)
    ctx.add_init_script("""window.__lt=[]; try{ new PerformanceObserver(l=>{ l.getEntries().forEach(e=>window.__lt.push(e.duration)); }).observe({entryTypes:['longtask']}); }catch(e){}""")
    pg = ctx.new_page()
    pg.on('console', lambda m: errs.append('console.' + m.type + ': ' + m.text) if m.type == 'error' and 'ERR_TUNNEL' not in m.text and 'ERR_NAME' not in m.text else None)
    pg.on('pageerror', lambda e: errs.append('pageerror: ' + str(e)))
    pg.goto('file://' + os.path.abspath(html))
    pg.wait_for_timeout(800)
    pg.click('#nav button[data-v="xps"]')
    pg.wait_for_selector('.xps', timeout=10000)
    J = lambda js, *a: pg.evaluate(js, *a)

    def idle(t=180000):
        pg.wait_for_timeout(150)
        pg.wait_for_function("document.querySelector('#xpsProg').hidden", timeout=t)
        pg.wait_for_timeout(150)

    def setv(sel, v):   # 숫자 칸·선택 칸에 값 → change
        pg.eval_on_selector(sel, "(e,v)=>{ e.value=v; e.dispatchEvent(new Event('change',{bubbles:true})); }", str(v))
        pg.wait_for_timeout(120)

    def chk(sel, on=True):
        pg.eval_on_selector(sel, "(e,v)=>{ if(e.checked!==v){ e.checked=v; e.dispatchEvent(new Event('change',{bubbles:true})); } }", on)
        pg.wait_for_timeout(120)

    def dl(click_sel):
        with pg.expect_download(timeout=20000) as d:
            pg.click(click_sel)
        f = d.value.path(); return open(f, 'rb').read()

    # ================= 1 Series (SDP VAMAS) =================
    print('== 1 Series: SDP depth profile', flush=True)
    pg.set_input_files('#xpsFin', [os.path.join(D, 'depth_profile_sdp.vms')]); pg.wait_for_timeout(800)
    check('SDP: 24 spectra (12 levels × Si 2p, O 1s)', J("XPS_UI._state().specs.length") == 24)
    pg.click('#xpsWs button[data-ws="series"]'); pg.wait_for_timeout(400)
    B = J("(()=>{ const D=XPS_UI._ser.build(); return {label:D.label,unit:D.unit,regions:D.regions.map(r=>({region:r.region,coord:r.coord}))}; })()")
    check('series by VAMAS variable "Etch Time" (s), coordinates = generator', B['label'] == 'Etch Time' and all(r['coord'] == T_dp['levels'] for r in B['regions']), f"{B['label']} [{B['unit']}] {B['regions'][0]['coord']}")
    # 생성기의 at% 는 n = A/S (투과·IMFP 보정 없음) — 같은 정의로 비교
    J("(()=>{ const S=XPS_UI._state(); S.quant.imfp={model:'none'}; S.quant.tf={model:'none'}; })()")
    for reg, inf in T_dp['regionInfo'].items():   # 생성기와 같은 창·S
        k = f'ser|reg|{reg}|'
        setv(f'[data-st="{k}lo"]', inf['lo']); setv(f'[data-st="{k}hi"]', inf['hi']); setv(f'[data-st="{k}S"]', inf['S'])
    pg.click('#xsrRun'); pg.wait_for_timeout(400)
    R = J("(()=>{ const R=XPS_UI._res().ser; return {rows:R.rows,at:R.comp.map(l=>l.map(q=>q.at)),s:R.comp.map(l=>l.map(q=>q.sAt))}; })()")
    e = 0
    for j, lv in enumerate(R['at']):
        for k, nm in enumerate(R['rows']):
            t = T_dp['at'][j][T_dp['regions'].index(nm)]
            if lv[k] is not None: e = max(e, abs(lv[k] - t))
    check('composition vs etch time within 3 at.% of the generator (Shirley areas)', e < 3, f'max |error| {e:.2f} at.%')
    check('σ(at %) reported per level', all(x is not None and x > 0 for x in R['s'][0]), f"σ level 0: {R['s'][0]}")
    csv = dl('#xsrFigC [data-fx="csv"]').decode('utf-8-sig')
    check('composition CSV', csv.count('\n') >= 13 and 'Si 2p_at_pct' in csv, csv.splitlines()[0][:120])
    svg = dl('#xsrFigC [data-fx="svg"]').decode()
    check('composition SVG (variables frozen)', svg.startswith('<svg') and 'var(--' not in svg, f'{len(svg)} bytes')
    # 정보 깊이 보정
    setv('[data-st="ser|rate"]', 0.02); chk('[data-st="ser|idc"]'); pg.click('#xsrRun'); pg.wait_for_timeout(300)
    idc = J("(()=>{ const R=XPS_UI._res().ser; return R.idc?{n:R.idc.at.length,lam:R.idc.lam,x1:R.x[1]}:null; })()")
    check('information-depth correction on a depth scale', idc and idc['n'] == 2 and abs(idc['x1'] - 0.6) < 1e-9 and all(l > 0 for l in idc['lam']), idc)
    chk('[data-st="ser|idc"]', False); setv('[data-st="ser|rate"]', '')
    # PCA (Si 2p) → 요인 수
    pg.click('#xsrPca'); idle()
    P = J("(()=>{ const M=XPS_UI._res().mv; return {nn:M.o.NF.nNoise,nIND:M.o.mal.nIND,nF:M.o.mal.nF,n:XPS_UI._state().ser.mv.n,eig:Array.from(M.o.eig).slice(0,4)}; })()")
    check('PCA: noise-limited factor count = 2 (Si⁰, Si⁴⁺; level-independent background)', P['nn'] == 2 and P['n'] == 2, P)
    pg.click('#xsrRes'); idle()
    M = J("(()=>{ const P=XPS_UI._ser.profiles(); const r=XPS_UI._res().mv.res.o; return {lof:r.lof,f0:Array.from(P.frac[0]),f1:Array.from(P.frac[1])}; })()")
    # 표면(0 s)은 산화물, 깊은 곳은 Si⁰: 두 성분의 분율이 반대로 바뀌어야 한다
    flip = (M['f0'][0] - 50) * (M['f0'][-1] - 50) < 0 and abs(M['f0'][0] - M['f0'][-1]) > 80
    check('MCR-ALS (2 comp.) separates oxide / substrate along the profile', flip and M['lof'] < 2, f"LOF {M['lof']:.3f} %, comp1 {M['f0'][0]:.1f} → {M['f0'][-1]:.1f} %")
    pg.screenshot(path=os.path.join(out, 'series.png'), full_page=True)

    # ================= 2 여러 열 텍스트 =================
    print('== 2 Multi-column text', flush=True)
    pg.set_input_files('#xpsFin', [os.path.join(D, 'depth_profile_multicol.txt')]); pg.wait_for_timeout(800)
    sets = J("[...document.querySelectorAll('[data-st=\"ser|set\"] option')].map(o=>o.value)")
    mk = [s for s in sets if 'multicol' in s]
    check('multi-column text → one data set (12 spectra)', bool(mk) and J("XPS_UI._state().specs.filter(s=>s.meta&&s.meta.mcol).length") == 12, sets)
    setv('[data-st="ser|set"]', mk[0]); pg.wait_for_timeout(300)
    B2 = J("(()=>{ const D=XPS_UI._ser.build(); return {label:D.label,unit:D.unit,reg:D.regions.map(r=>r.region),coord:D.regions[0].coord}; })()")
    check('column headers → "Etch Time" coordinate 0…480 s, region Si 2p', B2['label'] == 'Etch Time' and B2['unit'] == 's' and B2['coord'][-1] == 480 and B2['reg'] == ['Si 2p'], B2)
    pg.click('#xsrPca'); idle()
    nn2 = J("XPS_UI._res().mv.o.NF.nNoise")
    pg.click('#xsrRes'); idle()
    f2 = J("(()=>{ const P=XPS_UI._ser.profiles(); return P.frac.map(a=>[a[0],a[a.length-1]]); })()")
    check('text series: PCA + MCR-ALS run (factor count ≥ 2, fractions swap)', nn2 is not None and nn2 >= 2 and abs(f2[0][0] - f2[0][1]) > 50, f'n_noise {nn2}, {f2}')

    # ================= 3 ARXPS =================
    print('== 3 ARXPS', flush=True)
    pg.click('#xpsWs button[data-ws="spectra"]'); pg.wait_for_timeout(200)
    pg.set_input_files('#xpsFin', [os.path.join(D, 'arxps_angles.vms')]); pg.wait_for_timeout(800)
    # 0° 스펙트럼을 두 성분(Si⁰, Si⁴⁺)으로 피팅 → 그 성분이 기준 모양
    J("""(()=>{ const S=XPS_UI._state(); const s=S.specs.find(x=>x.sample.indexOf('1.8 nm')>=0&&x.meta.vamas.takeOffPolar===0); S.sel=s.uid; S.ui.tab='fit'; const p=XPS_UI._pr(s); const C=XPS_CORE;
      p.lo=95; p.hi=108.5; p.bg={type:'shirley',n:5}; const a=C.newComp('GLP',99.3,{F:0.9,m:0,name:'Si0',area:4000}), b=C.newComp('GLP',103.5,{F:1.5,m:0,name:'Si4',area:8000}); a.slot=0; b.slot=1;
      p.comps=[a,b]; p.cons={}; p.cons[a.id]={m:{mode:'fixed'}}; p.cons[b.id]={m:{mode:'fixed'}}; XPS_UI.rerender(); })()""")
    J("XPS_UI._fit('one')"); idle()
    check('0° spectrum fitted (2 components)', J("!!(XPS_UI._pr().result&&!XPS_UI._pr().result.stale)"))
    pg.click('#xpsWs button[data-ws="arxps"]'); pg.wait_for_timeout(300)
    sets = J("[...document.querySelectorAll('[data-st=\"arx|set\"] option')].map(o=>o.value)")
    setv('[data-st="arx|set"]', [s for s in sets if '1.8 nm' in s][0])
    th = J("Array.from(document.querySelectorAll('.xws-side .muted.xs')).map(e=>e.textContent).join(' | ')")
    setv('[data-st="arx|lo"]', 95); setv('[data-st="arx|hi"]', 108.5)
    for cb in J("[...document.querySelectorAll('[data-st^=\"arx|refs|c:\"]')].map(e=>e.dataset.st)"):
        chk(f'[data-st="{cb}"]')
    pg.click('#xarGet'); pg.wait_for_timeout(300)
    T = J("(()=>{ const T=XPS_UI._res().arx.tab; return {names:T.names,th:T.th,n:T.I.length}; })()")
    check('angle series → 8 angles from the "Angle" variable, 2 species', T['th'] == T_ar['angles'] and T['n'] == 2, T)
    so, ss = [i for i, n in enumerate(T['names']) if n.startswith('Si4')][0], [i for i, n in enumerate(T['names']) if n.startswith('Si0')][0]
    setv('[data-st="arx|o"]', so); setv('[data-st="arx|s"]', ss)
    setv('[data-st="arx|Lo"]', T_ar['Lo']); setv('[data-st="arx|Ls"]', T_ar['Ls']); setv('[data-st="arx|Rinf"]', T_ar['Rinf'])
    pg.click('#xarOv'); pg.wait_for_timeout(300)
    O = J("(()=>{ const r=XPS_UI._res().arx.ov; return {d:r.r.d,sd:r.r.sdScaled,phi:r.r.phi,hill:r.h.d}; })()")
    check('overlayer thickness d ≈ 1.8 nm (truth 1.8)', abs(O['d'] - T_ar['d']) < max(3 * O['sd'], 0.05), O)
    check('Hill plot drawn (biased, L_o ≠ L_s)', J("!!document.querySelector('#xarFigH svg')") and O['hill'] > 1.5, f"Hill d {O['hill']:.3f}")
    est = J("[XPS_UI._arx.al(0),XPS_UI._arx.al(1)]")
    check('attenuation-length calculator (Seah–Dench inorganic, SiO₂ defaults)', all(v and 1 < v < 6 for v in est), est)
    rdp = J("(()=>{ const R=XPS_UI._res().arx.rdp; return R.order; })()")
    check('relative depth plot: oxide nearer the surface than Si⁰', rdp[0] == so, rdp)
    pg.click('details:has(#xarMl) > summary'); pg.wait_for_timeout(200)
    setv(f'[data-st="arx|ml|sp|{so}|layer"]', 0); setv(f'[data-st="arx|ml|sp|{ss}|layer"]', 1)
    setv(f'[data-st="arx|ml|sp|{so}|L"]', T_ar['Lo']); setv(f'[data-st="arx|ml|sp|{ss}|L"]', T_ar['Ls'])
    setv(f'[data-st="arx|ml|sp|{so}|Iinf"]', T_ar['Rinf']); setv(f'[data-st="arx|ml|sp|{ss}|Iinf"]', 1)
    pg.click('#xarMl'); pg.wait_for_timeout(300)
    ml = J("(()=>{ const r=XPS_UI._res().arx.ml.r; return {d:r.d[0],sd:r.sdScaled[0],phi:r.phi}; })()")
    check('multilayer (1 layer, angle factor profiled) agrees with the overlayer fit', abs(ml['d'] - T_ar['d']) < max(3 * ml['sd'], 0.1), ml)
    pg.click('#xarMe'); idle(300000)
    me = J("(()=>{ const r=XPS_UI._res().arx.me.r; const c=r.c[%d]; return {alpha:r.alpha,chi2:r.chi2,N:r.N,nr:r.notReached,top:c[0],deep:c[c.length-1]}; })()" % so)
    check('MaxEnt profile (Worker): oxide-rich surface, substrate at depth, α by discrepancy', me['top'] > 0.6 and me['deep'] < 0.3 and me['alpha'] > 0 and not me['nr'], me)
    pg.screenshot(path=os.path.join(out, 'arxps.png'), full_page=True)

    # ================= 4 Imaging =================
    print('== 4 Imaging', flush=True)
    pg.click('#xpsWs button[data-ws="spectra"]'); pg.wait_for_timeout(200)
    pg.set_input_files('#xpsFin', [os.path.join(D, 'map_c1s.vms')]); pg.wait_for_timeout(1000)
    st = J("(()=>{ const s=XPS_UI._ims().stacks[0]; return {nx:s.nx,ny:s.ny,nE:s.nE,ws:XPS_UI._state().ui.ws,inList:XPS_UI._state().specs.filter(x=>x.meta.vamas&&x.meta.vamas.mode==='MAP').length}; })()")
    check('VAMAS MAP → 8×6 stack in Imaging (not 48 list entries)', st['nx'] == 8 and st['ny'] == 6 and st['ws'] == 'img' and st['inList'] == 0, st)
    setv('[data-st="imu|lo"]', 281); setv('[data-st="imu|hi"]', 293); setv('[data-st="imu|map"]', 'area'); setv('[data-st="imu|bg"]', 'shirley')
    pg.click('#ximRun'); idle()
    A = J("Array.from(XPS_UI._res().img.maps[0].v)"); sA = J("Array.from(XPS_UI._res().img.maps[0].s)")
    tot = [a + bb for a, bb in zip(T_map['areaA'], T_map['areaB'])]
    z = [(x - t) / s_ for x, t, s_ in zip(A, tot, sA)]; zm = sum(z) / len(z); zs = math.sqrt(sum((q - zm) ** 2 for q in z) / (len(z) - 1))
    check('Shirley area map vs generator (A + B): z = (area − truth)/σ', max(abs(q) for q in z) < 4 and abs(zm) < 0.6 and 0.5 < zs < 1.6, f'max |z| {max(abs(q) for q in z):.2f}, mean {zm:.2f}, sd {zs:.2f}, max rel. {100 * max(abs(x - t) / t for x, t in zip(A, tot)):.1f} % (low-signal px)')
    chk('[data-st="imu|sig"]'); pg.wait_for_timeout(200)
    check('σ map shown', J("XPS_UI._res().img.maps[0].s && document.querySelector('#ximFigM .t').textContent.indexOf('σ')>=0"))
    chk('[data-st="imu|sig"]', False)
    png = dl('#ximFigM [data-fx="png"]'); check('map PNG', png[:4] == b'\x89PNG', f'{len(png)} bytes')
    csv = dl('#ximFigM [data-fx="csv"]').decode('utf-8-sig'); check('map CSV (48 px)', len(csv.strip().splitlines()) == 50, csv.splitlines()[1])
    # ROI 끌기 (지도 전체)
    box = pg.query_selector('#ximFigM canvas.xim').bounding_box()
    pg.mouse.move(box['x'] + 2, box['y'] + 2); pg.mouse.down(); pg.mouse.move(box['x'] + box['width'] / 2, box['y'] + box['height'] / 2, steps=4); pg.mouse.move(box['x'] + box['width'] - 2, box['y'] + box['height'] - 2, steps=4); pg.mouse.up(); pg.wait_for_timeout(400)
    roi = J("(()=>{ const r=XPS_UI._res().roi; return r?{n:r.n,rect:r.rect}:null; })()")
    check('ROI drawn by dragging → mean spectrum', roi and roi['n'] == 48, roi)
    pg.click('#ximRoiS'); pg.wait_for_timeout(300)
    J("""(()=>{ const S=XPS_UI._state(); const s=S.specs[S.specs.length-1]; S.sel=s.uid; const p=XPS_UI._pr(s); const C=XPS_CORE; p.lo=281; p.hi=293; p.bg={type:'shirley',n:3};
      const a=C.newComp('GLP',284.8,{F:1.1,m:0,name:'A',area:2000}), b=C.newComp('GLP',288.6,{F:1.3,m:0,name:'B',area:1000}); a.slot=0; b.slot=1; p.comps=[a,b]; p.cons={}; p.cons[a.id]={m:{mode:'fixed'}}; p.cons[b.id]={m:{mode:'fixed'}}; })()""")
    J("XPS_UI._fit('one')"); idle()
    check('ROI spectrum added to the list and fitted', J("!!(XPS_UI._pr().result&&!XPS_UI._pr().result.stale)"))
    setv('[data-st="imu|map"]', 'nnls'); pg.wait_for_timeout(200)
    for cb in J("[...document.querySelectorAll('[data-st^=\"imu|refs|c:\"]')].map(e=>e.dataset.st)"):
        chk(f'[data-st="{cb}"]')
    pg.click('#ximRun'); idle()
    N = J("(()=>{ const m=XPS_UI._res().img.maps; return {a:Array.from(m.find(q=>q.key==='c0').v),b:Array.from(m.find(q=>q.key==='c1').v),sa:Array.from(m.find(q=>q.key==='c0').s)}; })()")
    eA = max(abs(x - t) / t for x, t in zip(N['a'], T_map['areaA'])); eB = max(abs(x - t) / t for x, t in zip(N['b'], T_map['areaB']))
    check('NNLS chemical-state maps (fitted ROI components as references) vs generator', eA < 0.12 and eB < 0.2, f'max rel. error A {100 * eA:.1f} %, B {100 * eB:.1f} %')
    check('state-map σ finite', all(math.isfinite(v) and v > 0 for v in N['sa']))
    setv('[data-st="imu|map"]', 'pca'); pg.click('#ximRun'); idle()
    pc = J("(()=>{ const P=XPS_UI._res().pca; return {nn:P.NF&&P.NF.nNoise,k:P.P.U.length}; })()")
    check('image PCA + MLPCA factor count (2 states + varying background → 3)', pc['nn'] is not None and 2 <= pc['nn'] <= 4, pc)
    n0 = J("XPS_UI._ims().stacks.length"); pg.click('#ximRec'); idle()
    check('PCA reconstruction → new stack', J("XPS_UI._ims().stacks.length") == n0 + 1)
    pg.click('#ximBin'); pg.wait_for_timeout(300)
    bs = J("(()=>{ const s=XPS_UI._ims().stacks[XPS_UI._ims().stacks.length-1]; return [s.nx,s.ny]; })()")
    check('2×2 binning → 4×3 stack', bs == [4, 3], bs)
    # 텍스트 영상 stack: 에너지 + 48 열 (행 우선) → 처음엔 48×1, nx = 8 로 바꾸기
    tx = J("""(()=>{ const s=XPS_UI._ims().stacks[0]; const L=['Binding Energy (eV)\\t'+[...Array(48)].map((_,p)=>'px'+p).join('\\t')]; for(let e=0;e<s.nE;e++){ const r=[s.E[e].toFixed(3)]; for(let p=0;p<48;p++) r.push(s.data[p*s.nE+e]); L.push(r.join('\\t')); } return L.join('\\n'); })()""")
    fp = os.path.join(tmp, 'c1s_stack.txt'); open(fp, 'w').write(tx)
    pg.set_input_files('#ximFin', [fp]); pg.wait_for_timeout(800)
    sh = J("(()=>{ const s=XPS_UI._ims().stacks[XPS_UI._ims().stacks.length-1]; return [s.nx,s.ny,s.nE]; })()")
    setv('#ximNx', 8); pg.wait_for_timeout(200)
    sh2 = J("(()=>{ const s=XPS_UI._ims().stacks[XPS_UI._ims().stacks.length-1]; return [s.nx,s.ny]; })()")
    check('multi-column text image stack, reshaped to nx = 8', sh[0] * sh[1] == 48 and sh2 == [8, 6], f'{sh} → {sh2}')
    pg.screenshot(path=os.path.join(out, 'imaging.png'), full_page=True)

    # ================= 5 Process =================
    print('== 5 Processing', flush=True)
    pg.click('#xpsWs button[data-ws="spectra"]'); pg.wait_for_timeout(200)
    pg.set_input_files('#xpsFin', [os.path.join(tmp, '[Au foil] VB.txt')]); pg.wait_for_timeout(600)
    J("(()=>{ const S=XPS_UI._state(); const s=S.specs.find(x=>x.sample==='Au foil'); S.sel=s.uid; S.ui.tab='proc'; XPS_UI.rerender(); })()"); pg.wait_for_timeout(300)
    pg.click('details:has(#xprFer) > summary'); pg.wait_for_timeout(100)
    pg.click('#xprFer'); pg.wait_for_timeout(500)
    F = J("(()=>{ const s=XPS_UI._cur(); const x=XPS_UI._pres()[s.uid].fermi; return {EF:x.EF.value,sEF:x.EF.sigma,fw:x.res.fwhm,sfw:x.res.sigma}; })()")
    check('Fermi edge: E_F ≈ 0 eV and resolution ≈ 0.40 eV (within 3σ)', abs(F['EF']) < max(3 * F['sEF'], 0.01) and abs(F['fw'] - FW) < max(3 * F['sfw'], 0.03), F)
    pg.click('#xprFerU'); pg.wait_for_timeout(200)
    check('Use as PSF → stored Gaussian PSF', J("XPS_UI._state().prc.psf==='stored'&&Math.abs(XPS_UI._state().prc.stored.fwhm-%.12f)<1e-9" % F['fw']))
    # 디컨볼루션: 깊은 단계의 Si 2p (Si⁰ 하나 + 작은 산화물), Gauss PSF 0.5 eV
    J("""(()=>{ const S=XPS_UI._state(); const s=S.specs.filter(x=>x.sample==='SiO2 on Si (synthetic)'&&x.region==='Si 2p').pop(); S.sel=s.uid; const p=XPS_UI._pr(s); p.lo=96; p.hi=106; S.prc.psf='gauss'; S.prc.fwhm=0.5; S.prc.method='jansson'; S.prc.mc=6; XPS_UI.rerender(); })()""")
    pg.wait_for_timeout(200); pg.click('#xprDec'); idle()
    Dd = J("(()=>{ const s=XPS_UI._cur(); const d=XPS_UI._pres()[s.uid].dec; return {k:d.k,stop:d.stop,gain:d.res.gain,ng:d.ng,area:d.area,band:!!d.sig}; })()")
    check('deconvolution (Worker): resolution gain > 1, noise gain and σ band, area conserved', Dd['gain'] > 1.05 and Dd['ng'] is not None and Dd['band'] and abs(Dd['area'] - 1) < 0.02, Dd)
    check('deconvolved overlay + "display only" label on the main plot', J("document.querySelector('#xpsPlot').textContent.indexOf('Deconvolved')>=0 && !!document.querySelector('#xpsPanel .badge.wr[title*=\"start values\"]')"))
    setv('[data-st="prc|rstep"]', 0.1); J("XPS_UI._proc.resample()"); pg.wait_for_timeout(200)
    rs = J("(()=>{ const x=XPS_UI._pres()[XPS_UI._cur().uid].rs; return {n:x.E.length,m:x.method,st:x.E[1]-x.E[0]}; })()")
    check('resample preview (0.05 → 0.1 eV, flux-conserving rebin)', rs['m'] == 'rebin' and abs(rs['st'] - 0.1) < 1e-9, rs)
    setv('[data-st="prc|der"]', 1); J("XPS_UI._proc.smooth()"); pg.wait_for_timeout(200)
    check('Savitzky–Golay 1st derivative overlay', J("document.querySelector('#xpsPlot').textContent.indexOf('derivative')>=0"))
    pg.screenshot(path=os.path.join(out, 'process.png'), full_page=False)

    # ================= 6 freeze · 390 px · reload =================
    lt = J("window.__lt")
    check('no main-thread task > 2 s', max(lt or [0]) < 2000, f'longest {max(lt or [0]):.0f} ms over {len(lt)} long tasks')
    pg.set_viewport_size({'width': 390, 'height': 844}); pg.wait_for_timeout(300)
    ov = {}
    for w in ['spectra', 'series', 'arxps', 'img']:
        pg.click(f'#xpsWs button[data-ws="{w}"]'); pg.wait_for_timeout(400)
        ov[w] = J("document.documentElement.scrollWidth - window.innerWidth")
        pg.screenshot(path=os.path.join(out, f'm390_{w}.png'), full_page=True)
    check('390 px: no horizontal overflow in any workspace', all(v <= 1 for v in ov.values()), ov)
    pg.set_viewport_size({'width': 1500, 'height': 1000})
    pg.wait_for_timeout(1200)
    pg.reload(); pg.wait_for_timeout(900); pg.click('#nav button[data-v="xps"]'); pg.wait_for_timeout(1500)
    J("XPS_UI._state().ui.ws='img'; XPS_UI.rerender()"); pg.wait_for_timeout(800)
    check('after reload: settings and image stacks restored (IndexedDB)', J("XPS_UI._ims().stacks.length") >= 4 and J("XPS_UI._state().ser.reg['O 1s'].S") == T_dp['regionInfo']['O 1s']['S'], J("XPS_UI._ims().stacks.length"))
    check('no console errors', not errs, errs[:5])
    b.close()

res['errors'] = errs
json.dump(res, open(os.path.join(out, 'ui_xps_multi.json'), 'w'), ensure_ascii=False, indent=1)
print(f'\n{len(res) - 1 - len(bad)}/{len(res) - 1} passed' + (' — FAIL: ' + ', '.join(bad) if bad else ''))
sys.exit(1 if bad else 0)
