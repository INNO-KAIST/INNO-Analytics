"""v1.9.8 audit fixes — browser (Chromium) check, one verdict per item (SPEC §17.9)
usage: python3 test/ui_xps_v198.py <injected html> <figure folder>
1 old project (saved before v1.9.2/v1.9.5/v1.9.8): VGD spectra without the 'TF removed' flag, IMFP n 0.6, quant.imfp without mat, no state
  version → opens without errors, Quantification renders (no 'reading rho' crash), defaults deep-merged, state version set; a one-line notice
  for the TF (results not changed: the flag stays unknown) and for n; 'Re-read VGD files' path (load the same .VGD again) refreshes only the
  TF metadata → notice gone; 'Use 0.7' sets n and clears the notice
2 lines by X-ray source: main line = strongest measurable line (BE < hν − 100 eV) of the source's Scofield table (Se 3p not 2p; Gd 3d under
  Al Kα but 4d under Mg Kα; Ho/Er/Tm 4d), strong lines (Si 2p and 2s), lineChoices from the Mg table for Mg Kα; quant interference warning for
  Pt 4f inside an Al 2p window (missed when Pt's 'main line' was 4d)
3 La 3d recipe with 'Fix ratio' → derived c4f¹L marked 'fixed', no ±
4 Ba 4d data shifted by +1.3 eV beyond the ±1 eV alignment → Fit check 'alignment at its search limit' warning
5 no console errors"""
import sys, os, json, tempfile, math, random, glob
from playwright.sync_api import sync_playwright

html, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
errs, res, bad = [], {}, []


def check(name, cond, info=''):
    res[name] = {'ok': bool(cond), 'info': info}
    if not cond:
        bad.append(name)
    print(('  ok   ' if cond else '  FAIL ') + name + ('  (' + str(info) + ')' if info != '' else ''), flush=True)


random.seed(11)
tmp = tempfile.mkdtemp()
LN2 = math.log(2)


def gls(x, F, m=0.3):
    g = (2 / F) * math.sqrt(LN2 / math.pi) * math.exp(-4 * LN2 * x * x / (F * F))
    l = (2 / (math.pi * F)) / (1 + 4 * x * x / (F * F))
    return (1 - m) * g + m * l


def pois(l):
    return max(0.0, round(l + math.sqrt(max(l, 1e-9)) * random.gauss(0, 1)))


def write(name, be, ys):
    fn = os.path.join(tmp, name)
    with open(fn, 'w') as fo:
        fo.write('Binding Energy (eV)\tCounts/s\n')
        for x, y in sorted(zip(be, ys), reverse=True):
            fo.write(f'{x:.3f}\t{y:.3f}\n')
    return fn


# 2: Al 2p window with Pt 4f7/2 (71.2) / 4f5/2 (74.5) inside it; a Pt 4f region of the same sample
ba = [66 + 0.05 * i for i in range(301)]
f_al = write('[V198P] Al2p.txt', ba, [pois(1500 + 8000 * gls(e - 74.4, 1.4) + 9000 * gls(e - 71.2, 1.0) + 7000 * gls(e - 74.5, 1.0)) for e in ba])
f_pt = write('[V198P] Pt4f.txt', ba, [pois(1500 + 9000 * gls(e - 71.2, 1.0) + 7000 * gls(e - 74.5, 1.0)) for e in ba])
# 3: La 3d, La₂O₃-like (c4f¹L 50 %): main 834.9, bonding +4.9, antibonding +3.7 (equal areas), ΔE 16.8, 2/3
LA = [(834.9, 1.6, 30000), (839.8, 2.0, 15000), (843.5, 2.0, 15000)]
bl = [822 + 0.05 * i for i in range(1081)]
f_la = write('[V198L] La3d.txt', bl, [pois(3000 + sum(A * (0.6 * gls(e - b0, F) + 0.4 * gls(e - b0 - 16.8, F)) for b0, F, A in LA)) for e in bl])
# 4: Ba 4d BaCO₃ at the recipe start (90.1 − 0.6 = 89.5) + 1.3 eV, ΔE 2.6, 2/3
bb = [80 + 0.05 * i for i in range(401)]
f_ba = write('[V198B] Ba4d.txt', bb, [pois(2000 + 30000 * (0.6 * gls(e - 90.8, 1.7) + 0.4 * gls(e - 93.4, 1.7))) for e in bb])

vgd = sorted(f for f in glob.glob(D + '/k2312/*.VGD') if 'Survey' not in f)

with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={'width': 1500, 'height': 950})
    pg.on('console', lambda m: errs.append('console.' + m.type + ': ' + m.text) if m.type == 'error' and 'ERR_TUNNEL' not in m.text and 'ERR_NAME' not in m.text else None)
    pg.on('pageerror', lambda e: errs.append('pageerror: ' + str(e)))
    pg.goto('file://' + os.path.abspath(html))
    pg.wait_for_timeout(800)
    pg.click('#nav button[data-v="xps"]')
    pg.wait_for_selector('.xps', timeout=10000)
    J = lambda js, *a: pg.evaluate(js, *a)

    def idle(t=180000):
        pg.wait_for_function("document.querySelector('#xpsProg').hidden", timeout=t)
        pg.wait_for_timeout(250)

    def quant_tab():
        J("document.querySelector('.xps [data-t=\"quant\"]').click()"); pg.wait_for_timeout(600)

    # ---------- 1 migration ----------
    pg.set_input_files('#xpsFin', vgd)
    pg.wait_for_function(f'XPS_UI._state().specs.length>={len(vgd)}', timeout=30000)
    snap = J("JSON.stringify({app:'INNO Analytics XPS',v:1,state:XPS_UI._state()},(k,v)=>ArrayBuffer.isView(v)?{__f64:Array.from(v)}:v)")
    o = json.loads(snap)
    st = o['state']
    for sp in st['specs']:
        sp['meta'].pop('tfRemoved', None)
        sp['meta'].pop('tfMode', None)
    st.pop('sv', None)
    st['quant']['imfp'] = {'model': 'tpp2m', 'n': 0.6}     # no 'mat' (crashed the panel: reading 'rho')
    st['quant'].pop('tf', None)
    old = os.path.join(tmp, 'old_project.json')
    open(old, 'w').write(json.dumps(o))
    pg.reload(); pg.wait_for_timeout(800); pg.click('#nav button[data-v="xps"]'); pg.wait_for_selector('.xps'); pg.wait_for_timeout(800)
    errs.clear()
    pg.set_input_files('#xpsFin', old); pg.wait_for_timeout(1500)
    quant_tab()
    s1 = J("""(()=>{ const st=XPS_UI._state(); const m=document.querySelector('#xqMig'); return {sv:st.sv, mat:st.quant.imfp.mat, tf:st.quant.tf, n:st.quant.imfp.n,
        unknown:st.specs.filter(s=>s.format==='VGD'&&s.meta.tf&&!('tfRemoved' in s.meta)).length, table:!!document.querySelector('.xq-tbl tbody tr'),
        mig:m?m.textContent:null, tfLine:!!document.querySelector('#xqMig [data-mig=tf]'), nLine:!!document.querySelector('#xqMig [data-mig=imfp]'), re:!!document.querySelector('#xqMigRe'),
        msg:(document.querySelector('#xpsMsg')||{}).textContent||''}; })()""")
    check('1a old project opens: no errors, Quantification table renders, defaults deep-merged (imfp.mat, quant.tf), state version 2',
          not errs and s1['table'] and s1['mat'] == {'rho': 5, 'Nv': 24, 'M': 150, 'Eg': 3} and s1['tf'] == {'model': 'auto', 'x': 1} and s1['sv'] == 2, [errs[:2], s1['sv'], s1['mat']])
    check('1b notices: VGD TF flag unknown (not changed silently: flag still unknown) with a re-read button, and IMFP n 0.6 vs default 0.7; one line each',
          s1['unknown'] > 0 and s1['tfLine'] and s1['re'] and s1['nLine'] and s1['n'] == 0.6 and all(len(x) < 200 for x in (s1['mig'] or '').split('Re-read')), (s1['unknown'], (s1['mig'] or '')[:160]))
    pg.set_input_files('#xpsFin', vgd); pg.wait_for_timeout(1500)
    quant_tab()
    s1c = J("""(()=>{ const st=XPS_UI._state(); return {n:st.specs.length, unknown:st.specs.filter(s=>s.format==='VGD'&&s.meta.tf&&!('tfRemoved' in s.meta)).length,
        tfLine:!!document.querySelector('#xqMig [data-mig=tf]'), rm:document.querySelectorAll('[data-tfrm]').length}; })()""")
    check('1c re-reading the same .VGD files refreshes only the TF metadata (no duplicate spectra): flag known, notice gone, "TF removed" badges shown',
          s1c['n'] == len(vgd) and s1c['unknown'] == 0 and not s1c['tfLine'] and s1c['rm'] > 0, s1c)
    J("document.querySelector('#xqMigN').click()"); pg.wait_for_timeout(400)
    s1d = J("({n:XPS_UI._state().quant.imfp.n, old:XPS_UI._state().quant.imfpOld, line:!!document.querySelector('#xqMig [data-mig=imfp]')})")
    m2 = J("""(()=>{ const r=XPS_UI._migrate({v:1,specs:[],quant:{imfp:{model:'power'}},ui:{tab:'quant'}}); return {mat:r.st.quant.imfp.mat, n:r.st.quant.imfp.n, sv:r.st.sv, notes:r.notes.length, survey:r.st.ui.survey}; })()""")
    check("1d 'Use 0.7' sets n and clears the notice; the IndexedDB path uses the same migration (missing nested defaults filled, no notice when nothing differs)",
          s1d['n'] == 0.7 and s1d['old'] is None and not s1d['line'] and m2['mat'] == {'rho': 5, 'Nv': 24, 'M': 150, 'Eg': 3} and m2['n'] == 0.7 and m2['sv'] == 2 and m2['notes'] == 0 and m2['survey'] == {}, [s1d, m2])
    pg.screenshot(path=os.path.join(out, 'v198_01_migration.png'), full_page=False)

    # ---------- 2 lines by source ----------
    s2 = J("""(()=>{ const L=XPS_UI._lines; const o=(e,h)=>{ const m=L.mainLine(e,h); return m?m.orb:null; };
        return {se:o('Se',1486.6), ho:o('Ho',1486.6), gdAl:o('Gd',1486.6), gdMg:o('Gd',1253.6), er:o('Er',1486.6), tm:o('Tm',1486.6), ba:o('Ba',1486.6), co:o('Co',1486.6),
          si:L.strongLines('Si',1486.6).map(x=>x.orb), pt:L.strongLines('Pt',1486.6).map(x=>x.orb),
          chMg:L.lineChoices('Ho','4d',1253.6), chAl:L.lineChoices('Ho','4d',1486.6), seMg:L.lineChoices('Se','3d',1253.6),
          sMg:L.lineStrength('Ba','3d',1253.6), sAl:L.lineStrength('Ba','3d',1486.6)}; })()""")
    check('2a main line = strongest measurable line of the source (BE < hν − 100 eV): Se 3p, Gd 3d (Al) / 4d (Mg), Ho/Er/Tm 4d (3d3/2 above hν − 100), Ba 3d, Co 2p; strong lines Si 2s+2p, Pt 4d+4f+4p',
          s2['se'] == '3p' and s2['ho'] == '4d' and s2['gdAl'] == '3d' and s2['gdMg'] == '4d' and s2['er'] == '4d' and s2['tm'] == '4d' and s2['ba'] == '3d' and s2['co'] == '2p' and set(s2['si']) == {'2s', '2p'} and {'4d', '4f'} <= set(s2['pt']),
          {k: s2[k] for k in ('se', 'ho', 'gdAl', 'gdMg', 'er', 'tm', 'si', 'pt')})
    check('2b lineChoices / lineStrength use the Mg Kα table for Mg Kα (no Ho 3d, no Se 2p — not in that table or above hν); Al keeps Ho 3d',
          '3d' not in s2['chMg'] and '3d' in s2['chAl'] and '2p' not in s2['seMg'] and s2['sMg'] != s2['sAl'], [s2['chMg'][:6], s2['seMg'][:6], s2['sMg'], s2['sAl']])
    pg.set_input_files('#xpsFin', [f_al, f_pt, f_la, f_ba])
    pg.wait_for_function(f'XPS_UI._state().specs.length>={len(vgd) + 4}', timeout=20000)
    regs = J("XPS_UI._state().specs.map(s=>s.region)")

    def select(sub):
        J("(r)=>{ const st=XPS_UI._state(); const x=st.specs.find(q=>(q.file||'').includes(r)); st.sel=x.uid; st.ui.tab='fit'; st.ui.selComp=null; XPS_UI.rerender(); }", sub)
        pg.wait_for_timeout(200)

    select('[V198P] Al2p')
    quant_tab()
    s2c = J("Array.from(document.querySelectorAll('.xps .note.warn')).map(x=>x.textContent).filter(t=>/Pt/.test(t))")
    check('2c quant interference: Pt 4f7/2/4f5/2 inside the Al 2p window is reported (old rule used Pt 4d only and missed it)',
          any('Pt 4f' in t and 'Al' in t for t in s2c), s2c[:2])

    # ---------- 3 La 3d Fix ratio → 'fixed' ----------
    select('[V198L]')
    J("XPS_UI._recipe('la3d',{species:['La2O3'],flags:{fixR:true,auger:false,plasmon:false}})"); pg.wait_for_timeout(300)
    J("XPS_UI._fit('one')"); idle()
    s3 = J("(()=>{ const d=document.querySelector('#xpsChk [data-der=ct]'); return d?{t:d.textContent, cls:d.className, tip:d.title}:null; })()")
    check("3 La 3d with 'Fix ratio' → c4f¹L shown as 'fixed' (no ±), tooltip says it is an input",
          s3 and 'fixed' in s3['t'] and '±' not in s3['t'] and 'fx' in s3['cls'] and 'not a fit result' in s3['tip'], s3 and s3['t'])

    # ---------- 4 Ba 4d alignment at the search limit ----------
    select('[V198B]')
    J("XPS_UI._recipe('ba4d',{species:['BaCO3'],flags:{}})"); pg.wait_for_timeout(300)
    J("XPS_UI._fit('one')"); idle()
    s4 = J("(()=>{ const it=Array.from(document.querySelectorAll('#xpsChk li')); const m=it.find(l=>/search limit/.test(l.textContent)); return {al:XPS_UI._pr().recipe&&XPS_UI._pr().recipe.align, it:m?{t:m.querySelector('.tx').textContent, warn:m.classList.contains('warn'), tip:m.querySelector('.tx').title}:null}; })()")
    check('4 Ba 4d shifted +1.3 eV: alignment stops at +1.0 eV → Fit check warning (charge reference / Ba species), label < 90 characters',
          s4['it'] and s4['it']['warn'] and abs(s4['al'] - 1.0) < 1e-9 and len(s4['it']['t']) < 90 and 'charge' in s4['it']['tip'], s4)
    pg.screenshot(path=os.path.join(out, 'v198_02_ba4d.png'), full_page=False)
    b.close()

check('5 no console errors', not errs, errs[:3])
json.dump(res, open(os.path.join(out, 'ui_xps_v198.json'), 'w'), ensure_ascii=False, indent=1)
print(f"\n{len(res) - len(bad)} / {len(res)} ok")
sys.exit(1 if bad else 0)
