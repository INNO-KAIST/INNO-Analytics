"""CasaXPS/Avantage 대응 기능 — 브라우저(Chromium) 실측, 항목마다 판정
사용: python3 test/ui_xps_casa.py <주입된 html> <그림 폴더>
1 기준 저장(Store as reference) ×2 → 다른 시료에서 LCF: add all → Fit → 분율·공선성 κ 표시
2 Monte Carlo(기준 잡음 포함 표시)   3 Apply model to other spectra → Apply & fit all(대전 보정이 다른 깊이 분포 4 개) → Batch results 표·CSV
4 측정 선모양(+ Shape, 위치 자유) → 이동 δ ≈ 0   5 Spline 배경: + Knot(그래프 클릭) · 끌기 · 두 번 클릭 지우기 · 모형 민감도 매듭 흔들기
6 Simplex = Marquardt (χ²ν)   7 정량 at% Monte Carlo 95 % 구간 표시·CSV 열"""
import sys, os, json, tempfile, math, random
from playwright.sync_api import sync_playwright

html, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
errs, res, bad = [], {}, []


def check(name, cond, info=''):
    res[name] = {'ok': bool(cond), 'info': info}
    if not cond:
        bad.append(name)
    print(('  ok   ' if cond else '  FAIL ') + name + ('  (' + str(info) + ')' if info != '' else ''), flush=True)


tmp = tempfile.mkdtemp()
random.seed(11)


def gls(x, c, F, m=0.2):   # 넓이 1 GL 합
    s = F / (2 * math.sqrt(2 * math.log(2)))
    g = math.exp(-0.5 * ((x - c) / s) ** 2) / (s * math.sqrt(2 * math.pi))
    l = (2 / (math.pi * F)) / (1 + 4 * ((x - c) / F) ** 2)
    return (1 - m) * g + m * l


def pois(l):
    return max(0.0, round(l + math.sqrt(max(l, 1e-9)) * random.gauss(0, 1)))


def write(name, xs, ys):
    p = os.path.join(tmp, name)
    with open(p, 'w') as fo:
        fo.write('Binding Energy (eV)\tCounts/s\n')
        for x, y in sorted(zip(xs, ys), reverse=True):
            fo.write(f'{x:.3f}\t{y:.3f}\n')
    return p


be = [96 + 0.05 * i for i in range(241)]
metal = lambda e: gls(e, 99.4, 0.9, 0.15)
oxide = lambda e: gls(e, 103.4, 1.6, 0.15)


def withShirley(P, base=900.0, k=0.02):
    Q, out_ = 0.0, []
    for i, p in enumerate(P):
        if i:
            Q += 0.5 * (P[i] + P[i - 1]) * 0.05
        out_.append(base + k * Q + p)
    return out_


files = [write('[RefMetal] Si2p.txt', be, [pois(v) for v in withShirley([400000 * metal(e) for e in be], 2000)]),
         write('[RefOxide] Si2p.txt', be, [pois(v) for v in withShirley([400000 * oxide(e) for e in be], 2000)])]
FR = [0.8, 0.6, 0.4, 0.2]
Q = [0.0, 0.7, -0.5, 1.1]   # 측정 BE = 참 + q, 시료 대전 보정 Δ = −q
for k, (fm, q) in enumerate(zip(FR, Q)):
    files.append(write(f'[Depth-{k + 1}] Si2p.txt', [e + q for e in be], [pois(v) for v in withShirley([30000 * (fm * metal(e) + (1 - fm) * oxide(e)) for e in be])]))
bo = [525 + 0.05 * i for i in range(281)]
files.append(write('[Depth-1] O1s.txt', bo, [pois(1500 + 40000 * gls(e, 532.4, 1.5)) for e in bo]))
bc = [524 + 0.05 * i for i in range(281)]
files.append(write('[Curvy] O1s.txt', bc, [pois(2500 + 500 * math.sin((e - 524) / 4.5) + 15 * (e - 524) + 15000 * gls(e, 531, 1.4, 0.0)) for e in bc]))

with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={'width': 1500, 'height': 950}, accept_downloads=True)
    pg.on('console', lambda m: errs.append('console.' + m.type + ': ' + m.text) if m.type == 'error' and 'ERR_TUNNEL' not in m.text and 'ERR_NAME' not in m.text else None)
    pg.on('pageerror', lambda e: errs.append('pageerror: ' + str(e)))
    pg.goto('file://' + os.path.abspath(html))
    pg.wait_for_timeout(800)
    pg.click('#nav button[data-v="xps"]')
    pg.wait_for_selector('.xps', timeout=10000)
    J = lambda js, *a: pg.evaluate(js, *a)
    pg.set_input_files('#xpsFin', files)
    pg.wait_for_function('window.XPS_UI && XPS_UI._state().specs.length>=8', timeout=20000)
    specs = J('XPS_UI._state().specs.map(s=>[s.sample,s.region])')
    print('   spectra:', specs)

    def select(sample, region):
        J("(a)=>{ const st=XPS_UI._state(); const x=st.specs.find(q=>q.sample.includes(a[0])&&q.region===a[1]); st.sel=x.uid; st.ui.selComp=null; st.ui.tab='fit'; XPS_UI.rerender(); }", [sample, region])
        pg.wait_for_timeout(150)

    def idle(t=180000):
        pg.wait_for_function("document.querySelector('#xpsProg').hidden", timeout=t)
        pg.wait_for_timeout(200)

    def fit():
        pg.click('#xpsFit')
        pg.wait_for_function('XPS_UI._pr().result && !XPS_UI._pr().result.stale', timeout=180000)
        idle()

    msg = lambda: pg.inner_text('#xpsMsg .cur') if pg.query_selector('#xpsMsg .cur') else ''
    open_det = lambda sel: J(f"(()=>{{ const d=document.querySelector('{sel}'); if(d&&!d.open) d.open=true; }})()")

    # 시료 대전 보정(직접): Depth-k 는 Δ = −q
    J("(Q)=>{ const S=XPS_UI._state(); S.specs.forEach(s=>{ const m=/Depth-(\\d)/.exec(s.sample); if(m){ S.samples[s.sample]=Object.assign(S.samples[s.sample]||{target:284.8},{shift:-Q[+m[1]-1],refNote:'test'}); } }); XPS_UI.rerender(); }", Q)

    # ---- 1 기준 저장 ×2
    for nm in ('RefMetal', 'RefOxide'):
        select(nm, 'Si 2p')
        pg.select_option('#xpsBg', 'shirley'); pg.wait_for_timeout(150)
        open_det('#xpsRefD'); pg.click('#xpsRefS'); pg.wait_for_timeout(200)
    refs = J('XPS_UI._refs().map(t=>({name:t.name,eref:t.eref,n:t.x.length,area:t.area0,raw:!!t.raw}))')
    check('1 Store as reference ×2 (background-subtracted, unit area, raw kept for MC)', len(refs) == 2 and abs(refs[0]['eref'] - 99.4) < 0.05 and abs(refs[1]['eref'] - 103.4) < 0.05 and all(r['raw'] for r in refs), [(r['name'], round(r['eref'], 3)) for r in refs])
    select('Depth-1', 'Si 2p')
    pg.select_option('#xpsBg', 'active'); pg.wait_for_timeout(150)
    open_det('#xpsRefD'); pg.wait_for_timeout(100)
    pg.click('#xpsLcf'); pg.wait_for_timeout(200)
    cs = J("(()=>{ const p=XPS_UI._pr(); return p.comps.map(c=>({shape:c.shape,be:c.be,fixed:(p.cons[c.id]||{}).be&&(p.cons[c.id]||{}).be.mode})); })()")
    check('1 LCF: add all → 2 fixed-position reference components', len(cs) == 2 and all(c['shape'] == 'MEAS' and c['fixed'] == 'fixed' for c in cs), cs)
    fit()
    r = J('(()=>{ const r=XPS_UI._pr().result; return {pct:r.summ.map(q=>q.areaPct), sig:r.summ.map(q=>q.sig.areaPct), lcf:r.lcf, chi:r.stats.chi2nu}; })()')
    strip = pg.inner_text('#xpsPanel .strip')
    check('1 LCF fit: metal fraction 80 % recovered, collinearity κ shown', abs(r['pct'][0] - 80) < 3 and r['lcf'] and 'Collinearity' in strip, f"{r['pct'][0]:.2f} ± {r['sig'][0]:.2f} %, κ {r['lcf']['cond']:.2f}, χ²ν {r['chi']:.2f}")
    pg.screenshot(path=os.path.join(out, 'c01_lcf_fit.png'), full_page=True)

    # ---- 2 Monte Carlo (기준 잡음 포함)
    pg.fill('#xpsMCN', '40'); pg.dispatch_event('#xpsMCN', 'change')
    pg.click('#xpsMC'); pg.wait_for_function('XPS_UI._pr().result.mc', timeout=180000); idle()
    badge = pg.inner_text('#xpsPanel .xps-h')
    tbl = pg.inner_text('#xpsPanel table.xps-tbl')
    check('2 Monte Carlo with reference noise (badge) · fraction column', 'incl. ref. noise' in badge and 'fraction' in tbl, badge.replace('\n', ' ')[:90])
    pg.screenshot(path=os.path.join(out, 'c02_lcf_mc.png'), full_page=True)

    # ---- 3 Apply model → batch fit (다른 Depth 3 개; 기준 시료 2 개는 끔)
    J("XPS_UI._state().ui.tab='fit'; XPS_UI.rerender()"); pg.wait_for_timeout(150)
    open_det('#xpsBatD'); pg.wait_for_timeout(100)
    boxes = J("[...document.querySelectorAll('[data-bat]')].map(i=>({u:i.dataset.bat,c:i.checked,l:i.parentElement.textContent.trim()}))")
    for bx in boxes:
        uid = bx['u']
        smp = J("(u)=>XPS_UI._state().specs.find(s=>s.uid===u).sample", uid)
        want = 'Depth' in smp and J("(u)=>XPS_UI._state().specs.find(s=>s.uid===u).region", uid) == 'Si 2p'
        if bool(bx['c']) != want:
            pg.click(f'[data-bat="{uid}"]'); pg.wait_for_timeout(80)
    nsel = J("document.querySelectorAll('[data-bat]:checked').length")
    pg.click('#xpsBatF')
    pg.wait_for_function("XPS_UI._state().ui.batch && XPS_UI._state().ui.tab==='result'", timeout=180000); idle()
    rows = J("""(()=>{ const S=XPS_UI._state(); const b=S.ui.batch; return [b.src].concat(b.uids).map(u=>{ const s=S.specs.find(x=>x.uid===u); const r=XPS_UI._pr(s).result; const sh=S.samples[s.sample].shift||0;
          return {smp:s.sample, pct:r&&!r.stale?r.summ[0].areaPct:null, be:r?r.summ.map(q=>q.be+sh):null}; }); })()""")
    okf = all(rw['pct'] is not None and abs(rw['pct'] - 100 * FR[int(rw['smp'].split('Depth-')[1][0]) - 1]) < 3 for rw in rows)
    okb = all(rw['be'] and abs(rw['be'][0] - 99.4) < 0.08 and abs(rw['be'][1] - 103.4) < 0.08 for rw in rows)
    check('3 Apply & fit all: 3 selected (refs excluded by user), fractions per depth ≈ truth, corrected BE aligned', nsel == 3 and len(rows) == 4 and okf and okb, [(rw['smp'], round(rw['pct'] or -1, 1)) for rw in rows])
    bt = pg.inner_text('#xpsPanel table.xps-bat')
    with pg.expect_download() as dl:
        pg.click('#xpsBatCsv')
    csv = open(dl.value.path(), encoding='utf-8-sig').read().strip().split('\n')
    check('3 Batch results table + CSV (header + 4 rows, at% column)', 'Depth-4' in bt and len(csv) == 5 and 'at_pct' in csv[0] and 'BE_corr_eV' in csv[0], csv[0][:100])
    pg.screenshot(path=os.path.join(out, 'c03_batch.png'), full_page=True)

    # ---- 4 측정 선모양(+ Shape, 위치 자유)
    select('Depth-3', 'Si 2p')
    J("(()=>{ const p=XPS_UI._pr(); p.comps=[]; p.cons={}; p.ver++; XPS_UI.rerender(); })()")
    open_det('#xpsRefD'); pg.wait_for_timeout(100)
    ids = J('XPS_UI._refs().map(t=>t.id)')
    for i in ids:
        pg.click(f'[data-refm="{i}"]'); pg.wait_for_timeout(150)
        open_det('#xpsRefD')
    fit()
    m4 = J("(()=>{ const S=XPS_UI._state(), s=XPS_UI._cur(), sh=S.samples[s.sample].shift; const p=XPS_UI._pr(); return p.result.comps.map((c,i)=>({d:c.be+sh-c.tpl.eref, s:p.result.summ[i].sig.be, pct:p.result.summ[i].areaPct})); })()")
    check('4 measured line shapes with free position: shift δ ≈ 0 (±3σ), fraction ≈ 40 %', all(abs(c['d']) < 3 * c['s'] + 0.01 for c in m4) and abs(m4[0]['pct'] - 40) < 3, [(round(c['d'], 4), round(c['s'], 4), round(c['pct'], 1)) for c in m4])
    card = pg.inner_text('#xpsComps')
    check('4 component card shows the reference name and shift', 'Measured' in J("document.querySelector('#xpsComps').innerHTML") or 'Shift δ' in card, card[:80].replace('\n', ' '))

    # ---- 5 Spline 배경
    select('Curvy', 'O 1s')
    pg.select_option('#xpsBg', 'spline'); pg.wait_for_timeout(200)

    def px(be_, yv=None):
        return J("(a)=>{ const G=XPS_UI._geo(); const svg=document.querySelector('#xpsSvg'); const r=svg.getBoundingClientRect(); const k=r.width/920; const M=G.M; let y=a[1]; if(y==null){ let i=0; for(let j=0;j<M.n;j++) if(Math.abs(M.E[j]-a[0])<Math.abs(M.E[i]-a[0])) i=j; y=M.y[i]; } return [r.left+G.X(a[0]+G.sh)*k, r.top+G.Y(y)*k]; }", [be_, yv])

    pg.click('#xpsKnA'); pg.wait_for_timeout(100)
    x, y = px(527.0)
    pg.mouse.click(x, y); pg.wait_for_timeout(250)
    kn = J('XPS_UI._pr().bg.knots')
    has = pg.query_selector('#xpsSvg rect.knot') is not None
    check('5 + Knot then click the graph → knot at the clicked BE (auto height), handle drawn', len(kn) == 1 and abs(kn[0]['e'] - 527) < 0.1 and 'y' not in kn[0] and has, kn)
    hx, hy = J("(()=>{ const r=document.querySelector('#xpsSvg rect.knot').getBoundingClientRect(); return [r.left+r.width/2, r.top+r.height/2]; })()")
    pg.mouse.move(hx, hy); pg.mouse.down(); pg.mouse.move(hx + 25, hy - 15, steps=4); pg.mouse.up(); pg.wait_for_timeout(200)
    kn2 = J('XPS_UI._pr().bg.knots')
    check('5 drag a knot → position and (manual) height change', abs(kn2[0]['e'] - 527) > 0.05 and 'y' in kn2[0], kn2)
    pg.click('#xpsKnR'); pg.wait_for_timeout(150)
    check('5 Auto heights resets dragged height', 'y' not in J('XPS_UI._pr().bg.knots')[0])
    hx, hy = J("(()=>{ const r=document.querySelector('#xpsSvg rect.knot').getBoundingClientRect(); return [r.left+r.width/2, r.top+r.height/2]; })()")
    pg.mouse.dblclick(hx, hy); pg.wait_for_timeout(200)
    check('5 double-click a knot removes it', len(J('XPS_UI._pr().bg.knots')) == 0)
    J('XPS_UI._addKnot(527.0)'); J('XPS_UI._addKnot(535.0)'); pg.wait_for_timeout(150)
    J('XPS_UI._auto()'); idle()
    J("(()=>{ const p=XPS_UI._pr(); p.comps=p.comps.slice(0,1); p.ver++; })()")
    start5 = J("JSON.stringify(XPS_UI._pr().comps)")
    fit()
    a5 = J('XPS_UI._pr().result.summ[0].area')
    check('5 spline background fit (2 knots + ends): area within 5 % of truth 15000', abs(a5 - 15000) / 15000 < 0.05, f'{a5:.0f}')
    pg.screenshot(path=os.path.join(out, 'c05_spline.png'), full_page=True)
    chi_lm = J('XPS_UI._pr().result.stats.chi2nu')
    J("XPS_UI._state().ui.tab='result'; XPS_UI.rerender()"); pg.wait_for_timeout(200)
    J("XPS_UI._state().ui.sens={mode:'oat',dw:0.5,dk:0.5,bgs:{smart:false,shirley:false,active:false,tougaard:false},ns:false}; XPS_UI.rerender()"); pg.wait_for_timeout(150)
    has_k = pg.query_selector('#xsnK') is not None
    pg.click('#xsnRun'); pg.wait_for_function('XPS_UI._pr().sens', timeout=180000); idle()
    sv = J("(()=>{ const s=XPS_UI._pr().sens; return {n:s.sum.n, kn:s.variants.filter(v=>v.f.kn).length, eff:s.sum.comps[0].area.eff.kn}; })()")
    det = J("document.querySelector('#xpsPanel').textContent")
    check('5 model sensitivity: knot ± jitter variants (4) with a knot main effect', has_k and sv['kn'] == 4 and sv['eff'] > 0 and 'knot 1' in det, sv)
    pg.screenshot(path=os.path.join(out, 'c05b_sens.png'), full_page=True)

    # ---- 6 Simplex
    J("XPS_UI._state().ui.tab='fit'; XPS_UI.rerender()"); pg.wait_for_timeout(150)
    pg.select_option('#xpsOpt', 'simplex'); pg.wait_for_timeout(150)
    J("(c)=>{ const p=XPS_UI._pr(); p.comps=JSON.parse(c); p.ver++; }", start5)   # 같은 출발점에서 (IRLS 2 차 분산이 출발점에 따라 조금 달라 χ²ν 를 비교하려면)
    fit()
    chi_nm = J('XPS_UI._pr().result.stats.chi2nu'); st6 = pg.inner_text('#xpsPanel .strip'); opt6 = J('XPS_UI._pr().result.optimizer')
    check('6 Simplex fit reaches the Marquardt minimum (χ²ν)', opt6 == 'simplex' and 'Simplex' in st6 and abs(chi_nm - chi_lm) / chi_lm < 1e-4, f'{chi_lm:.6f} vs {chi_nm:.6f}')

    # ---- 7 정량 at% MC
    select('Depth-1', 'O 1s')
    J('XPS_UI._auto()'); idle(); fit()
    J("XPS_UI._state().ui.tab='quant'; XPS_UI.rerender()"); pg.wait_for_timeout(300)
    mc = J("[...document.querySelectorAll('.xq-mc')].map(e=>e.textContent)")
    with pg.expect_download() as dl:
        pg.click('#xqCsv')
    qcsv = open(dl.value.path(), encoding='utf-8-sig').read().strip().split('\n')
    check('7 quantification: Monte Carlo 95 % at% interval per element + CSV columns', len(mc) == 2 and 'at_pct_mc95_lo' in qcsv[0], mc)
    pg.screenshot(path=os.path.join(out, 'c07_quant.png'), full_page=True)

    # ---- 8 CasaXPS 영역 규약 + GL 표기
    select('Curvy', 'O 1s')
    pg.select_option('#xpsBg', 'linear'); pg.wait_for_timeout(150)
    pg.select_option('#xpsConv', 'casa'); pg.wait_for_timeout(200)
    t8 = pg.inner_text('#xpsPanel')
    J("(()=>{ const p=XPS_UI._pr(); p.comps[0]=Object.assign(p.comps[0],{shape:'GLP',m:0.3,F:1.2}); delete p.comps[0].T; p.ver++; XPS_UI._state().ui.selComp=p.comps[0].id; XPS_UI.rerender(); })()"); pg.wait_for_timeout(200)
    t8b = pg.inner_text('#xpsComps')
    conv = J('XPS_UI._pr().bg')
    check('8 CasaXPS region convention: Av. width field, CasaXPS region Area/Max/Min shown; GLP card shows CasaXPS GL(30) + fwhm', conv.get('conv') == 'casa' and conv.get('aw') == 1 and 'CasaXPS region: Area' in t8 and 'GL(30)' in t8b, [l for l in t8.split('\n') if 'CasaXPS region' in l][:1])
    pg.screenshot(path=os.path.join(out, 'c08_casa_conv.png'), full_page=True)
    # ---- 9 Avantage 규약 Shirley + 꼬리 성분
    pg.select_option('#xpsBg', 'shirley'); pg.wait_for_timeout(150)
    pg.select_option('#xpsConv', 'avantage'); pg.wait_for_timeout(200)
    bg9 = J('XPS_UI._pr().bg'); has_it = pg.query_selector('#xpsAvIt') is not None
    cid = J('XPS_UI._pr().comps[0].id')
    pg.select_option(f'[data-shape="{cid}"]', 'AGT'); pg.wait_for_timeout(200)
    c9 = J('XPS_UI._pr().comps[0]'); html9 = J("document.querySelector('#xpsComps').innerHTML")
    fit()
    r9 = J('XPS_UI._pr().result.stats.chi2nu')
    check('9 Avantage Shirley (iterations field) + GL product with Avantage tail (ET/TM/CT, channel width from data) fits', bg9.get('conv') == 'avantage' and bg9.get('it') == 10 and has_it and c9['shape'] == 'AGT' and abs(c9['ch'] - 0.05) < 1e-9 and 'Avantage: Product + tail' in html9 and r9 > 0, f"χ²ν {r9:.3f} (curved background, Shirley misfit expected)")

    # ---- 390 px: 새 패널이 가로로 넘치지 않음
    pg.set_viewport_size({'width': 390, 'height': 844}); select('Depth-1', 'Si 2p'); open_det('#xpsRefD'); open_det('#xpsBatD'); pg.wait_for_timeout(300)
    ov = J("document.documentElement.scrollWidth-document.documentElement.clientWidth")
    check('390 px: no horizontal page overflow with reference / batch panels open', ov <= 1, ov)
    pg.screenshot(path=os.path.join(out, 'c08_390.png'), full_page=True)
    b.close()

check('no console errors', not errs, errs[:3])
json.dump(res, open(os.path.join(out, 'ui_xps_casa.json'), 'w'), ensure_ascii=False, indent=1)
print(f"\n{len(res) - len(bad)} / {len(res)} ok")
sys.exit(1 if bad else 0)
