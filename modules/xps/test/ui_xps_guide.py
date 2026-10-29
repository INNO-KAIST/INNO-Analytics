"""v1.9.3 fit guidance (SPEC §17) — browser (Chromium) check, one verdict per item
usage: python3 test/ui_xps_guide.py <injected html> <figure folder>
1 Ce 3d (synthetic, Ce³⁺ 20 %): satellite-region banner + Recipe entry next to Auto suggest; terse banner, details in title
2 naive Auto suggest → Fit → Fit check: warning badge, items with title tooltips, 'Apply recipe' action
3 Apply recipe (from the Fit check) → 5 pairs, in-fit Shirley, window 875–925 → Staged fit (Worker) → derived Ce³⁺ / u‴ shown, Ce³⁺ ≈ 20 %
4 C 1s with a left-out line → Fit → 'Missing intensity' + 'Add a line' → click → one more component → Fit → finding gone
5 Fit check line count stays short (labels < 90 characters); 390 px: no horizontal overflow; no console errors"""
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


random.seed(7)
tmp = tempfile.mkdtemp()
LN2 = math.log(2)


def gls(x, F, m=0.3):   # area-normalised GL sum (SPEC §1.2)
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


# Ce 3d: recipe lines (SPEC §17.3), Ce³⁺ fraction 0.2, ΔE 18.5, u/v 2/3, in-fit Shirley-type background (c0 + k∫P)
LINES = [('v0', 880.5, 2.4, 0.2 * 0.288), ('v', 882.4, 2.3, 0.8 * 0.296), ('v1', 885.0, 3.1, 0.2 * 0.712), ('v2', 888.9, 5.0, 0.8 * 0.405), ('v3', 898.2, 2.3, 0.8 * 0.299)]
be = [870 + 0.1 * i for i in range(601)]
P = [sum(150000 * fr * (0.6 * gls(e - b0, F) + 0.4 * gls(e - b0 - 18.5, F)) for _, b0, F, fr in LINES) for e in be]
Q, acc = [], 0.0
for i, p in enumerate(P):
    if i:
        acc += 0.5 * (p + P[i - 1]) * 0.1
    Q.append(acc)
f_ce = write('[GuideS] Ce3d.txt', be, [pois((2500 + 0.1 * q + p) * 2) / 2 for p, q in zip(P, Q)])
# C 1s with two lines (the second, 20 %, will be left out on purpose)
bc = [278 + 0.05 * i for i in range(361)]
f_c = write('[GuideC] C1s.txt', bc, [pois(2000 + 30000 * gls(e - 285, 1.2) + 6000 * gls(e - 287.2, 1.2)) for e in bc])

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

    def select(region):
        J("(r)=>{ const st=XPS_UI._state(); const x=st.specs.find(q=>q.region===r); st.sel=x.uid; st.ui.tab='fit'; st.ui.selComp=null; XPS_UI.rerender(); }", region)
        pg.wait_for_timeout(200)

    pg.set_input_files('#xpsFin', [f_ce, f_c])
    pg.wait_for_function('window.XPS_UI && XPS_UI._state().specs.length>=2', timeout=20000)
    regs = J("XPS_UI._state().specs.map(s=>s.region)")
    ce = [r for r in regs if 'Ce' in r][0]
    c1 = [r for r in regs if 'C 1s' in r or 'C1s' in r][0]
    select(ce)
    # 1 banner + recipe entry
    s1 = J("(()=>{ const b=document.querySelector('.xps-rcp-b'); const a=document.querySelector('#xpsRcpA'), o=document.querySelector('#xpsRcpO'); const auto=document.querySelector('#xpsAuto'); return {b:b&&b.textContent, bt:b&&b.title, a:!!a, opts:o?Array.from(o.options).map(x=>x.textContent):[], near:!!(a&&auto&&a.parentElement===auto.parentElement), at:a&&a.title}; })()")
    check('1 Ce 3d: one-line satellite-region banner (details in title) and Recipe select + button in the Auto suggest row',
          s1['b'] and len(s1['b']) < 90 and 'satellite' in s1['b'] and 'Burroughs' not in s1['b'] and len(s1['bt'] or '') > 200 and s1['a'] and s1['near'] and s1['opts'] == ['CeO₂ only (6)', 'Ce³⁺ only (4)', 'Mixed (10)'] and 'Pitfalls' in (s1['at'] or ''), s1['b'])
    # 2 naive Auto suggest → Fit → Fit check
    J("XPS_UI._auto()"); idle()
    J("XPS_UI._fit('one')"); idle()
    s2 = J("(()=>{ const c=document.querySelector('#xpsChk'); if(!c) return null; const it=Array.from(c.querySelectorAll('li')); return {badge:c.querySelector('.badge').textContent, warnBadge:c.querySelector('.badge').classList.contains('wr'), n:it.length, labels:it.map(l=>l.querySelector('.tx').textContent), titles:it.map(l=>l.querySelector('.tx').title.length), btns:it.map(l=>{ const b=l.querySelector('button'); return b?b.textContent:null; }), head:c.querySelector('.xps-chk-h span').title.length}; })()")
    check('2 naive fit → Fit check with a warning badge; every item has a tooltip; recipe offered as an action',
          s2 and s2['warnBadge'] and s2['n'] >= 2 and all(t > 60 for t in s2['titles']) and 'Apply recipe' in s2['btns'] and s2['head'] > 200, s2 and list(zip(s2['labels'], s2['btns'])))
    pg.screenshot(path=os.path.join(out, 'guide_01_naive.png'), full_page=False)
    # 3 Apply recipe from the Fit check → staged fit
    J("(()=>{ const b=Array.from(document.querySelectorAll('#xpsChk button')).find(x=>x.textContent==='Apply recipe'); b.click(); })()"); pg.wait_for_timeout(300)
    s3 = J("(()=>{ const p=XPS_UI._pr(), sp=XPS_UI._cur(); const sh=XPS_UI._state().samples; return {n:p.comps.length, names:p.comps.map(c=>c.name), bg:p.bg.type, slope:p.bg.cons&&p.bg.cons.s&&p.bg.cons.s.mode, lo:p.lo, hi:p.hi, dbl:p.comps.every(c=>c.doublet&&Math.abs(c.doublet.r-2/3)<1e-9), recipe:p.recipe&&p.recipe.label, msg:(document.querySelector('#xpsMsg .cur .tx')||{}).textContent}; })()")
    check('3a Apply recipe → 5 v/u pairs (u/v 2/3), in-fit Shirley + slope, window 875–925 eV, short note with a Staged fit button',
          s3['n'] == 5 and s3['dbl'] and s3['bg'] == 'active' and s3['slope'] == 'free' and abs(s3['lo'] - 875) < 0.06 and abs(s3['hi'] - 925) < 0.06 and s3['recipe'] == 'Mixed (10)' and 'Staged fit' in J("(document.querySelector('#xpsMsgA')||{}).textContent||''"), s3)
    J("document.querySelector('#xpsStaged').click()"); idle()
    s4 = J("(()=>{ const r=XPS_UI._pr().result; const c=document.querySelector('#xpsChk'); const d=c&&Array.from(c.querySelectorAll('[data-der]')).map(x=>({id:x.dataset.der, t:x.textContent, tip:x.title.length})); const st=c&&Array.from(c.querySelectorAll('.xps-chk-h span')).map(x=>x.textContent); return {staged:r&&r.staged&&r.staged.used, chi:r&&r.stats.chi2nu, d, st, labels:c?Array.from(c.querySelectorAll('li .tx')).map(x=>x.textContent):[], msg:(document.querySelector('#xpsMsg .cur .tx')||{}).textContent}; })()")
    ce3 = None
    for d in (s4['d'] or []):
        if d['id'] == 'ce3':
            try:
                ce3 = float(d['t'].split()[1])
            except Exception:
                pass
    check('3b Staged fit (Worker) → result kept and reported; derived Ce³⁺ and u‴ shown with tooltips; Ce³⁺ within ±5 % of the true 20 %',
          s4['staged'] in ('staged', 'normal') and s4['d'] and len(s4['d']) == 2 and all(x['tip'] > 80 for x in s4['d']) and ce3 is not None and abs(ce3 - 20) <= 5 and any('staged: kept' in t for t in s4['st']) and 'Staged fit done' in (s4['msg'] or ''),
          {'Ce3+': ce3, 'chi2nu': s4['chi'], 'kept': s4['staged'], 'items': s4['labels']})
    pg.screenshot(path=os.path.join(out, 'guide_02_recipe.png'), full_page=False)
    # 4 C 1s with a left-out line → add action
    select(c1)
    J("(()=>{ const p=XPS_UI._pr(); const C=window.XPS_CORE; const c=C.newComp('GLS',285,{F:1.2,area:30000}); c.name='C–C'; c.slot=0; p.comps=[c]; p.cons={}; p.bg={type:'active',n:5,c0:null,k:null,s:0,cons:{s:{mode:'free'}}}; p.ver++; XPS_UI.rerender(); })()")
    J("XPS_UI._fit('one')"); idle()
    s5 = J("(()=>{ const it=Array.from(document.querySelectorAll('#xpsChk li')); const m=it.find(l=>/Missing intensity/.test(l.textContent)); return {has:!!m, label:m&&m.querySelector('.tx').textContent, btn:m&&m.querySelector('button')&&m.querySelector('button').textContent}; })()")
    check("4a left-out 20 % line → 'Missing intensity near 287… eV' with an 'Add a line' button", s5['has'] and s5['btn'] and s5['btn'].startswith('Add a line') and '287' in s5['label'], s5)
    J("(()=>{ const m=Array.from(document.querySelectorAll('#xpsChk li')).find(l=>/Missing intensity/.test(l.textContent)); m.querySelector('button').click(); })()"); pg.wait_for_timeout(250)
    n_after = J("XPS_UI._pr().comps.length")
    J("XPS_UI._fit('one')"); idle()
    s6 = J("(()=>{ const p=XPS_UI._pr(); return {n:p.comps.length, be:p.comps.map(c=>+c.be.toFixed(2)), miss:Array.from(document.querySelectorAll('#xpsChk li')).some(l=>/Missing intensity/.test(l.textContent)), badge:(document.querySelector('#xpsChk .badge')||{}).textContent}; })()")
    check('4b Add a line → 2 components → Fit → the residual finding is gone', n_after == 2 and s6['n'] == 2 and not s6['miss'] and any(abs(x - 287.2) < 0.15 for x in s6['be']), s6)
    # 5 terse labels · 390 px · console
    labs = J("Array.from(document.querySelectorAll('#xpsChk li .tx')).map(x=>x.textContent)")
    select(ce)
    labs += J("Array.from(document.querySelectorAll('#xpsChk li .tx')).map(x=>x.textContent)")
    check('5a Fit check labels are one short line (< 90 characters)', labs and max(len(x) for x in labs) < 90, max((len(x) for x in labs), default=0))
    pg.set_viewport_size({'width': 390, 'height': 844}); pg.wait_for_timeout(400)
    J("XPS_UI.rerender()"); pg.wait_for_timeout(300)
    ov = J("document.documentElement.scrollWidth-document.documentElement.clientWidth")
    box = J("(()=>{ const r=document.querySelector('#xpsChk'); if(!r) return null; const b=r.getBoundingClientRect(); return {right:b.right, vw:innerWidth, sw:r.scrollWidth, cw:r.clientWidth}; })()")
    check('5b 390 px: no horizontal page overflow; Fit check box fits', ov <= 1 and box and box['right'] <= box['vw'] + 1 and box['sw'] <= box['cw'] + 1, [ov, box])
    pg.screenshot(path=os.path.join(out, 'guide_03_390.png'), full_page=True)
    b.close()

check('5c no console errors', not errs, errs[:3])
json.dump(res, open(os.path.join(out, 'ui_xps_guide.json'), 'w'), ensure_ascii=False, indent=1)
print(f"\n{len(res) - len(bad)} / {len(res)} ok")
sys.exit(1 if bad else 0)
