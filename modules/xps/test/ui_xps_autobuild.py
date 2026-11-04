"""v1.10 Auto build (SPEC §18) — browser (Chromium) check, one verdict per item
usage: python3 test/ui_xps_autobuild.py <injected html> <figure folder>
1 'Auto build' button next to Auto suggest with a tooltip (rules, thresholds, limits)
2 region without components → first proposal 'Start model: Auto suggest' — nothing changes before Accept; Accept → start model + next proposal
  with evidence chips (ΔAICc, ΔBIC, χ²ν …) whose title explains them; the components are unchanged while a proposal is shown
3 Accept → the model equals the proposal's model and the Undo history grows by exactly one entry
4 Reject → logged as rejected, its key excluded, a different proposal (or the stop) follows; the model is unchanged
5 Accept remaining (a new session if the chain already ended) → every step logged 'accepted (remaining)', one Undo entry each (+ the final fit); final summary card with an Established / Model not established
  badge (reasons in title), a lines table and σ; 'Undo' restores the model before the last step
6 step log: Markdown with the step table, the stop reason and the summary (download button present; text via the hook)
7 C 1s region → the first proposal is 'Start from the C 1s recipe' (charge-correction order C 1s → O 1s respected);
  O 1s of a sample whose C 1s has no recipe fit → accepting the recipe step stops with the C 1s-first message, nothing applied silently
8 390 px: no horizontal overflow, panel fits; no console errors"""
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


def shirley_like(be, P, step, k):
    Q, acc = [], 0.0
    for i, p in enumerate(P):
        if i:
            acc += 0.5 * (p + P[i - 1]) * step
        Q.append(acc)
    return [k * q for q in Q]


# S 2s (no recipe): three lines, the third a shoulder (sep 1.3 FWHM, 25 %)
bs = [222 + 0.05 * i for i in range(301)]
Ps = [30000 * gls(e - 228.0, 1.2) + 12000 * gls(e - 230.2, 1.25) + 7500 * gls(e - 231.75, 1.2) for e in bs]
Ks = shirley_like(bs, Ps, 0.05, 0.08)
f_s = write('[AB1] S2s.txt', bs, [pois(1500 + k + p) for p, k in zip(Ps, Ks)])
# C 1s and O 1s of another sample (C 1s not fitted with the recipe)
bc = [278 + 0.05 * i for i in range(361)]
f_c = write('[AB2] C1s.txt', bc, [pois(1200 + 30000 * gls(e - 284.8, 1.3) + 7000 * gls(e - 286.3, 1.3) + 3000 * gls(e - 288.8, 1.3)) for e in bc])
bo = [524 + 0.05 * i for i in range(281)]
Po = [30000 * gls(e - 529.3, 1.2) + 6000 * gls(e - 528.4, 1.2) + 16000 * gls(e - 531.5, 1.8) + 5000 * gls(e - 532.6, 1.8) for e in bo]
f_o = write('[AB2] O1s.txt', bo, [pois(1500 + p) for p in Po])

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

    def idle(t=300000):
        pg.wait_for_function("document.querySelector('#xpsProg').hidden", timeout=t)
        pg.wait_for_timeout(300)

    def select(region):
        J("(r)=>{ const st=XPS_UI._state(); const x=st.specs.find(q=>q.region===r); st.sel=x.uid; st.ui.tab='fit'; st.ui.selComp=null; XPS_UI.rerender(); }", region)
        pg.wait_for_timeout(200)

    def wait_prop(t=300000):   # a proposal, a stop or a final summary
        pg.wait_for_function("(()=>{ const p=XPS_UI._pr(); return p.ab&&(p.ab.prop||p.ab.final||p.ab.stop)&&document.querySelector('#xpsProg').hidden; })()", timeout=t)
        pg.wait_for_timeout(300)

    pg.set_input_files('#xpsFin', [f_s, f_c, f_o])
    pg.wait_for_function('window.XPS_UI && XPS_UI._state().specs.length>=3', timeout=20000)
    regs = J("XPS_UI._state().specs.map(s=>s.region)")
    rs = [r for r in regs if 'S' in r and '2s' in r][0]
    rc = [r for r in regs if 'C' in r and '1s' in r][0]
    ro = [r for r in regs if 'O' in r and '1s' in r][0]
    select(rs)
    s1 = J("(()=>{ const a=document.querySelector('#xpsAB'), u=document.querySelector('#xpsAuto'); return {a:!!a, t:a&&a.title, near:!!(a&&u&&a.parentElement===u.parentElement), txt:a&&a.textContent}; })()")
    check('1 Auto build button in the Auto suggest row, tooltip with the rules (ΔAICc/ΔBIC ≥ 10, N_eff, |r| > 0.9, not established)',
          s1['a'] and s1['near'] and s1['txt'] == 'Auto build' and all(k in (s1['t'] or '') for k in ['≥ 10', 'N_eff', '0.9', 'not established', 'Undo']), len(s1['t'] or ''))

    # 2 start: Auto suggest proposal; nothing applied before Accept
    J("document.querySelector('#xpsAB').click()"); wait_prop()
    s2a = J("(()=>{ const p=XPS_UI._pr(); return {n:p.comps.length, kind:p.ab.prop&&p.ab.prop.kind, lab:(document.querySelector('.xps-ab-p')||{}).textContent, btn:['#xpsABa','#xpsABr','#xpsABall','#xpsABs'].map(s=>!!document.querySelector(s))}; })()")
    J("document.querySelector('#xpsABa').click()"); wait_prop()
    s2 = J("(()=>{ const p=XPS_UI._pr(); const st=p.ab.prop; const ev=Array.from(document.querySelectorAll('.xps-ab-ev span')); return {n:p.comps.length, ids:p.comps.map(c=>c.id).join(), prop:st&&st.label, kind:st&&st.kind, mIds:st&&st.model?st.model.comps.map(c=>c.id).join():null, chips:ev.map(x=>x.textContent), tips:ev.map(x=>x.title.length), log:p.ab.steps.map(s=>s.decision+':'+s.kind)}; })()")
    check('2a no components → first proposal "Start model: Auto suggest" with Accept / Reject / Accept remaining / Stop; nothing applied yet',
          s2a['n'] == 0 and s2a['kind'] == 'suggest' and 'Auto suggest' in (s2a['lab'] or '') and all(s2a['btn']), s2a)
    check('2b Accept → start model, then a proposal with evidence chips (ΔAICc, ΔBIC, χ²ν) whose tooltips explain them; model still the start model',
          s2['n'] >= 1 and s2['prop'] and any('ΔAICc' in c for c in s2['chips']) and any('ΔBIC' in c for c in s2['chips']) and any('χ²ν' in c for c in s2['chips']) and min(s2['tips']) > 20
          and s2['log'] == ['accepted:suggest'], [s2['prop'], s2['chips']])
    pg.screenshot(path=os.path.join(out, 'ab_01_proposal.png'), full_page=False)

    # 3 Accept one step (the decision background) — exactly one Undo entry (two if the chain ended and the final fit ran)
    pre = J("(()=>{ const p=XPS_UI._pr(); return {h:(p.hist||[]).length, n:p.ab.steps.length, m:JSON.stringify(p.ab.prop.model.comps.map(c=>[c.id,+c.be.toFixed(6)])), cons:JSON.stringify(p.ab.prop.model.cons), bg:p.ab.prop.model.bg.type}; })()")
    J("document.querySelector('#xpsABa').click()"); wait_prop()
    s3 = J("""(pre)=>{ const p=XPS_UI._pr(); const fin=!!p.ab.final; const h=p.hist||[]; const e=fin?h[h.length-2]:h[h.length-1];
      return {dh:h.length-pre.h, fin, same:fin?true:JSON.stringify(p.comps.map(c=>[c.id,+c.be.toFixed(6)]))===pre.m, consSame:fin?true:JSON.stringify(p.cons)===pre.cons, bg:p.bg.type, dec:p.ab.steps[pre.n].decision, why:e&&e.why}; }""", pre)
    check('3 Accept → the model is the proposal\'s model, exactly one Undo entry for the step ("Auto build: …"), logged as accepted',
          s3['dec'] == 'accepted' and s3['same'] and s3['consSame'] and s3['bg'] == pre['bg'] and s3['dh'] == (2 if s3['fin'] else 1) and (s3['why'] or '').startswith('Auto build'), s3)

    # 4 Reject the next proposal
    s4 = None
    if J("!!XPS_UI._pr().ab.prop"):
        before = J("(()=>{ const p=XPS_UI._pr(); return {ids:p.comps.map(c=>c.id+':'+c.be.toFixed(4)).join(), h:(p.hist||[]).length, key:p.ab.prop.key, n:p.ab.steps.length}; })()")
        J("document.querySelector('#xpsABr').click()"); wait_prop()
        s4 = J("(b)=>{ const p=XPS_UI._pr(); return {ids:p.comps.map(c=>c.id+':'+c.be.toFixed(4)).join(), h:(p.hist||[]).length, ex:p.ab.excluded, key:p.ab.prop&&p.ab.prop.key, dec:p.ab.steps[b.n].decision, fin:!!p.ab.final}; }", before)
        s4['ok'] = s4['dec'] == 'rejected' and before['key'] in s4['ex'] and (s4['fin'] or (s4['ids'] == before['ids'] and s4['h'] == before['h'] and s4['key'] != before['key']))
    check('4 Reject → logged as rejected, key excluded, model and Undo history unchanged, a different proposal (or the final summary) follows',
          s4 is not None and s4['ok'], s4)

    # 5 Accept remaining → final summary
    if not J("!!XPS_UI._pr().ab.prop"):   # the chain already ended (final ran) → a new Auto build session on the same region
        J("document.querySelector('#xpsAB').click()"); wait_prop()
    h0 = J("(XPS_UI._pr().hist||[]).length")
    acc0 = J("XPS_UI._pr().ab.steps.filter(s=>/^accepted/.test(s.decision)&&s.kind!=='suggest').length")
    had_prop = J("!!(XPS_UI._pr().ab.prop&&XPS_UI._pr().ab.prop.kind!=='multi')")
    if had_prop:
        J("document.querySelector('#xpsABall').click()")
    pg.wait_for_function("(()=>{ const p=XPS_UI._pr(); return p.ab&&(p.ab.final||(p.ab.prop&&p.ab.prop.kind==='multi')); })()", timeout=600000); idle()
    if J("!!(XPS_UI._pr().ab.prop&&XPS_UI._pr().ab.prop.kind==='multi')"):
        J("document.querySelector('#xpsABa').click()"); pg.wait_for_function("!!XPS_UI._pr().ab.final", timeout=600000); idle()
    s5 = J("""(()=>{ const p=XPS_UI._pr(); const A=p.ab; const acc=A.steps.filter(s=>/^accepted/.test(s.decision)&&s.kind!=='suggest').length; const sum=document.querySelector('.xps-ab-sum'); const bd=sum&&sum.querySelector('.badge');
      return {h:(p.hist||[]).length, acc, rem:A.steps.filter(s=>s.decision==='accepted (remaining)').length, stop:A.stop, est:A.final.established, reasons:A.final.reasons, badge:bd&&bd.textContent, bt:bd&&bd.title, rows:sum?sum.querySelectorAll('tr').length:0, n:p.comps.length,
        res:!!(p.result&&!p.result.stale), lines:A.final.lines.map(l=>l.name+' '+l.be.toFixed(2)), steps:A.steps.map(s=>s.decision+':'+s.label.slice(0,40))}; })()""")
    check('5a Accept remaining → each step one Undo entry (+ final fit), stop reason logged, summary badge with title, lines table, fresh fit result',
          had_prop and s5['rem'] >= 1 and s5['h'] - h0 == s5['acc'] - acc0 + 1 and s5['stop'] and s5['est'] == (s5['badge'] == 'Established') and s5['badge'] in ('Established', 'Model not established') and len(s5['bt'] or '') > 40 and s5['rows'] == s5['n'] + 1 and s5['res'],
          [s5['badge'], s5['stop'], s5['lines'], s5['reasons']])
    # truth: 228.0 · 230.2 · 231.75 eV
    tr = [228.0, 230.2, 231.75]
    found = [any(abs(float(l.split()[-1]) - t) < 0.3 for l in s5['lines']) for t in tr]
    check('5b S 2s (three lines, shoulder at 1.3 FWHM / 25 %): all three found, no extra line', all(found) and len(s5['lines']) == 3, s5['lines'])
    pre_undo = J("XPS_UI._pr().comps.length")
    J("document.querySelector('#xpsUndo').click()"); pg.wait_for_timeout(300)
    s5c = J("(()=>{ const p=XPS_UI._pr(); return {n:p.comps.length, stale:!p.result||p.result.stale}; })()")
    check('5c Undo after Auto build restores the model before the final fit (one entry)', s5c['n'] == pre_undo, s5c)
    pg.screenshot(path=os.path.join(out, 'ab_02_summary.png'), full_page=False)

    # 6 log
    md = J("XPS_UI._ab.log()")
    check('6 step log (Markdown): step table, decisions, stop reason, result section; Log button present',
          md.startswith('# XPS Auto build log') and '| # | decision |' in md and 'rejected' in md and 'accepted' in md and 'Stopped:' in md and '## Result:' in md and J("!!document.querySelector('#xpsABlog')"), md[:200])

    # 7 recipes: C 1s first proposal; O 1s without the C 1s recipe → stop with the order message
    select(rc)
    J("document.querySelector('#xpsAB').click()"); wait_prop()
    s7a = J("(()=>{ const p=XPS_UI._pr(); return {kind:p.ab.prop&&p.ab.prop.kind, lab:p.ab.prop&&p.ab.prop.label, why:p.ab.prop&&p.ab.prop.why, n:p.comps.length}; })()")
    check('7a C 1s → first proposal "Start from the C 1s recipe" (mentions the charge correction), nothing applied yet',
          s7a['kind'] == 'recipe' and 'C 1s' in (s7a['lab'] or '') and 'charge correction' in (s7a['why'] or '') and s7a['n'] == 0, s7a)
    J("document.querySelector('#xpsABx').click()"); pg.wait_for_timeout(200)
    select(ro)
    J("document.querySelector('#xpsAB').click()"); wait_prop()
    s7b0 = J("(()=>{ const p=XPS_UI._pr(); return {kind:p.ab.prop&&p.ab.prop.kind, n:p.comps.length}; })()")
    J("document.querySelector('#xpsABa').click()"); pg.wait_for_timeout(1500); idle()
    s7b = J("(()=>{ const p=XPS_UI._pr(); return {stop:p.ab.stop, n:p.comps.length, rec:p.comps.some(c=>c.recipe)}; })()")
    check('7b O 1s whose C 1s has no recipe fit: recipe step → stops with the C 1s → charge correction → O 1s message, no recipe model applied',
          s7b0['kind'] == 'recipe' and s7b['stop'] and 'C 1s' in s7b['stop'] and not s7b['rec'], [s7b0, s7b])

    # 8 390 px
    select(rs)
    J("(()=>{ const p=XPS_UI._pr(); if(p.ab) p.ab.open=true; XPS_UI.rerender(); })()")
    pg.set_viewport_size({'width': 390, 'height': 844}); pg.wait_for_timeout(400)
    J("XPS_UI.rerender()"); pg.wait_for_timeout(300)
    ov = J("document.documentElement.scrollWidth-document.documentElement.clientWidth")
    box = J("(()=>{ const r=document.querySelector('#xpsABp'); if(!r) return null; const b=r.getBoundingClientRect(); return {right:b.right, vw:innerWidth, sw:r.scrollWidth, cw:r.clientWidth}; })()")
    check('8a 390 px: no horizontal page overflow; Auto build panel fits', ov <= 1 and box and box['right'] <= box['vw'] + 1 and box['sw'] <= box['cw'] + 1, [ov, box])
    pg.screenshot(path=os.path.join(out, 'ab_03_390.png'), full_page=True)
    b.close()

check('8b no console errors', not errs, errs[:3])
json.dump(res, open(os.path.join(out, 'ui_xps_autobuild.json'), 'w'), ensure_ascii=False, indent=1)
print(f"\n{len(res) - len(bad)} / {len(res)} ok")
sys.exit(1 if bad else 0)
