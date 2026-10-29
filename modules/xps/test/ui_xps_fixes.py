"""XPS 탭 결함 수정 회귀 (10-19 감사 10 건) — 브라우저(Chromium) 실측, 항목마다 판정
사용: python3 test/ui_xps_fixes.py <주입된 html> <그림 폴더>
1 템플릿 저장·불러오기(대전 보정 뒤) → 다시 피팅 χ²ν 유지   2 C 1s 맞추기: 클릭한 성분을 몰래 쓰지 않음, 고른 기준이 주 성분이 아니면 경고
3 피팅 중 수정 → 결과 버림   4 C 1s 맞추기 중 스펙트럼 전환   5 Auger: 직접 추가 기본 짝, 같은 스펙트럼·광전자 영역 거부
6 정량: 원소당 센 선만 기본 포함 + 창 안 다른 원소 주선 경고   7 머리줄 없는 KE 텍스트 → 감지·KE 로 다시 읽기(Info 탭)
9 쉼표 소수 · 배경 뺀 자료 · 좁은 창 · SNR 기본값·저장 · 고정 BE σ '—' · 'at start'   10 390 px: at% 열·그래프 글자"""
import sys, os, glob, json, tempfile, math
from playwright.sync_api import sync_playwright

html, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
errs, res, bad = [], {}, []


def check(name, cond, info=''):
    res[name] = {'ok': bool(cond), 'info': info}
    if not cond:
        bad.append(name)
    print(('  ok   ' if cond else '  FAIL ') + name + ('  (' + str(info) + ')' if info != '' else ''), flush=True)


tmp = tempfile.mkdtemp()


def gls(x, c, F, m=0.3):
    s = F / (2 * math.sqrt(2 * math.log(2)))
    return (1 - m) * math.exp(-0.5 * ((x - c) / s) ** 2) + m / (1 + ((x - c) / (F / 2)) ** 2)


def write(name, xs, ys, header=None, sep='\t', comma=False):
    p = os.path.join(tmp, name)
    with open(p, 'w') as fo:
        if header:
            fo.write(header + '\n')
        for x, y in zip(xs, ys):
            a, b = f'{x:.3f}', f'{y:.3f}'
            if comma:
                a, b = a.replace('.', ','), b.replace('.', ',')
            fo.write(a + sep + b + '\n')
    return p


import random
random.seed(5)
# 7: 머리줄 없는 KE 파일 (C 1s, KE = 1486.6 − 284.8 = 1201.8)
ke = [1190 + i * 0.05 for i in range(501)]
f_ke = write('[KE-1] C1s.txt', ke, [2000 + 30000 * gls(k, 1201.8, 1.2) + random.gauss(0, 40) for k in ke])
# 9: 쉼표 소수(세미콜론) · 배경 뺀 자료(음수 포함)
be = [525 + i * 0.05 for i in range(321)]
f_cd = write('[CD-1] O1s.csv', be[::-1], [3000 + 20000 * gls(b, 531.0, 1.4) for b in be[::-1]], header='Binding Energy (eV);Counts/s', sep=';', comma=True)
f_bs = write('[BS-1] O1s.txt', be[::-1], [8000 * gls(b, 531.0, 1.4) + random.gauss(0, 60) for b in be[::-1]], header='Binding Energy (eV)\tCounts/s')

with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={'width': 1500, 'height': 950})
    pg.on('console', lambda m: errs.append('console.' + m.type + ': ' + m.text) if m.type == 'error' and 'ERR_TUNNEL' not in m.text and 'ERR_NAME' not in m.text else None)
    pg.on('pageerror', lambda e: errs.append('pageerror: ' + str(e)))
    pg.goto('file://' + os.path.abspath(html))
    pg.wait_for_timeout(800)
    pg.click('#nav button[data-v="xps"]')
    pg.wait_for_selector('.xps', timeout=10000)
    pg.evaluate("localStorage.removeItem('innoanalytics_xps_tpl_v1')")
    xls = sorted(glob.glob(os.path.join(HERE, 'data', 'bct600', '*.xls')))
    pg.set_input_files('#xpsFin', xls)
    pg.wait_for_function('window.XPS_UI && XPS_UI._state().specs.length>=7', timeout=20000)
    J = lambda js: pg.evaluate(js)

    def select(region, fmt='Avantage XLS', sample=None):
        pg.evaluate("(a)=>{ const [r,f,s]=a; const st=XPS_UI._state(); const x=st.specs.find(q=>q.region===r&&q.format===f&&(!s||q.sample===s)); st.sel=x.uid; st.ui.selComp=null; XPS_UI.rerender(); }", [region, fmt, sample])
        pg.wait_for_timeout(150)

    def idle(t=120000):
        pg.wait_for_function("document.querySelector('#xpsProg').hidden", timeout=t)
        pg.wait_for_timeout(150)

    msg = lambda: pg.inner_text('#xpsMsg .cur') if pg.query_selector('#xpsMsg .cur') else ''

    def fit(kind='one'):
        J(f'XPS_UI._fit("{kind}")')
        pg.wait_for_function('XPS_UI._pr().result && !XPS_UI._pr().result.stale', timeout=120000)
        idle()

    # ---- 2 · 11: C 1s 자동 제안 → 피팅, strip / 맞추기는 선택 성분을 쓰지 않음
    select('C 1s')
    J('XPS_UI._auto()'); idle(); fit()
    strip = pg.inner_text('#xpsPanel .strip')
    check('11 strip: no "0 iter" after Auto suggest → Fit', ' 0 iter' not in strip and '\n0 iter' not in strip, strip.replace('\n', ' ')[:120])
    comps = J('XPS_UI._pr().result.comps.map(c=>({id:c.id,be:c.be,area:c.area,name:c.name}))')
    tot = sum(max(0, c['area']) for c in comps)
    auto = sorted([c for c in comps if c['area'] >= 0.1 * tot], key=lambda c: c['be'])[0]
    hi = max(comps, key=lambda c: c['be'])
    J(f"XPS_UI._state().ui.selComp='{hi['id']}'")
    J('XPS_UI._calib()'); pg.wait_for_timeout(400)
    sh = J('XPS_UI._state().samples[XPS_UI._cur().sample].shift')
    m1 = msg()
    check('2 calibrate ignores the clicked component (uses auto lowest-BE major)', abs(sh - (284.8 - auto['be'])) < 2e-3 and f"{auto['be']:.2f}" in m1, f"Δ {sh} · msg: {m1[:110]}")
    # 고른 기준(드롭다운)이 주 성분이 아니면 경고
    J("XPS_UI._state().ui.chgOpen=true; XPS_UI.rerender()"); pg.wait_for_timeout(150)
    pg.select_option('#xcComp', hi['id']); pg.wait_for_timeout(500)
    sh2 = J('XPS_UI._state().samples[XPS_UI._cur().sample].shift'); m2 = msg()
    check('2 explicit reference (dropdown) is used and warned when not the lowest-BE major C 1s', abs(sh2 - (284.8 - hi['be'])) < 2e-3 and '주의' in m2, f"Δ {sh2} · {m2[:120]}")
    pg.select_option('#xcComp', ''); pg.wait_for_timeout(500)
    sh = J('XPS_UI._state().samples[XPS_UI._cur().sample].shift')
    check('2 back to Auto reference', abs(sh - (284.8 - auto['be'])) < 2e-3, sh)

    # ---- 1: 템플릿 (Δ ≠ 0)
    chi0 = J('XPS_UI._pr().result.stats.chi2nu')
    J('XPS_UI._saveTpl()'); n_t = J('XPS_UI._tplList().length')
    J(f'XPS_UI._loadTpl({n_t - 1})'); pg.wait_for_timeout(200)
    bnd = J("(()=>{ const p=XPS_UI._pr(); return p.comps.map(c=>{ const b=(p.cons[c.id]||{}).be; return b&&b.mode==='free'&&b.lo!=null?[c.be-b.lo,b.hi-c.be]:null; }); })()")
    fit()
    chi1 = J('XPS_UI._pr().result.stats.chi2nu')
    check('1 template save → load (charge-corrected, Δ≠0) → refit keeps χ²ν', abs(chi1 - chi0) / chi0 < 0.02 and all(x is None or (x[0] >= -1e-6 and x[1] >= -1e-6) for x in bnd), f'Δ {sh:.3f} · χ²ν {chi0:.3f} → {chi1:.3f} · bounds {bnd}')

    # ---- 3: 피팅 중 수정 → 버림
    J("(()=>{ const p=XPS_UI._pr(); p.comps.forEach(c=>c._mark=1); XPS_UI._fit('multi'); setTimeout(()=>{ p.comps[0].area*=1.01; p.ver++; },30); })()")
    idle()
    kept = J('XPS_UI._pr().comps.every(c=>c._mark===1)'); m3 = msg()
    check('3 edit during fit → result discarded, edits kept', kept and '버렸습니다' in m3, m3[:90])

    # ---- 6 · 9(좁은 창): 정량 기본 포함 · 간섭 경고 · 좁은 창
    k = J('XPS_UI._cur().sample')
    rows = J(f"XPS_UI._quantRows({json.dumps(k)}).map(r=>({{el:r.el,line:r.line,inc:r.inc,region:r.sp.region}}))")
    co = {r['region']: r['inc'] for r in rows if r['el'] == 'Co'}
    check('6 quant: only the strongest Co line included by default', co.get('Co 2p') is True and co.get('Co 3p') is False, co)
    intf = J(f"XPS_UI._interf({json.dumps(k)})")
    check('6 quant: Ba 3d main line inside the Co 2p window is flagged', any('Co' in t and 'Ba' in t for t in intf), intf[:2])
    J("(()=>{ const S=XPS_UI._state(); const s=S.specs.find(x=>x.region==='Co 3p'&&x.format==='Avantage XLS'); const p=XPS_UI._pr(s); p.comps=[]; p.result=null; p.lo=s.be[10]; p.hi=s.be[13]; S.ui.tab='quant'; XPS_UI.rerender(); })()")
    pg.wait_for_timeout(300); qt = pg.inner_text('#xpsPanel')
    check('9 too-narrow fit window says so (not "No region spectra")', 'Fit range too narrow' in qt and 'No region spectra' not in qt)
    pg.screenshot(path=os.path.join(out, 'f06_quant.png'), full_page=True)

    # ---- 10 (fixed σ): Results 에 '± —'
    select('O 1s'); J("XPS_UI._state().ui.tab='fit'; XPS_UI.rerender()")
    J('XPS_UI._auto()'); idle()
    J("(()=>{ const p=XPS_UI._pr(); const c=p.comps[0]; (p.cons[c.id]=p.cons[c.id]||{}).be={mode:'fixed'}; p.ver++; })()")
    fit()
    J("XPS_UI._state().ui.tab='result'; XPS_UI.rerender()"); pg.wait_for_timeout(300)
    row0 = pg.inner_text('#xpsPanel table.xps-tbl tbody tr')
    check('9 fixed-BE component shows σ(BE) "—" in Results', '± —' in row0 and '± 0.00' not in row0.split('\t')[1], row0.replace('\t', ' | ')[:90])

    # ---- 4: C 1s 맞추기 중 다른 스펙트럼으로 전환 (K2312 VGD)
    J("XPS_UI._state().ui.tab='fit'; XPS_UI.rerender()")
    pg.set_input_files('#xpsFin', sorted(glob.glob(os.path.join(HERE, 'data', 'k2312', '*.VGD'))))
    pg.wait_for_function("XPS_UI._state().specs.filter(s=>s.format==='VGD').length>=6", timeout=20000); pg.wait_for_timeout(300)
    J("(()=>{ const S=XPS_UI._state(); const s=S.specs.filter(x=>x.format==='VGD'&&x.region==='Ba 4d').pop(); S.sel=s.uid; XPS_UI.rerender(); })()")
    pg.wait_for_timeout(200)
    pg.click('#xpsChgGo')
    pg.wait_for_timeout(60)
    J("(()=>{ const S=XPS_UI._state(); const s=S.specs.filter(x=>x.format==='VGD'&&x.region==='O 1s').pop(); S.sel=s.uid; XPS_UI.rerender(); })()")
    pg.wait_for_function("(()=>{ const S=XPS_UI._state(); const s=S.specs.find(x=>x.format==='VGD'); const sm=S.samples[s.sample]; return sm&&sm.ref; })()", timeout=90000)
    idle()
    st4 = J("(()=>{ const S=XPS_UI._state(); const v=S.specs.filter(x=>x.format==='VGD'); const g=r=>v.find(x=>x.region===r); const sm=S.samples[v[0].sample]; return {sel:XPS_UI._cur().region, refRegion:S.specs.find(x=>x.uid===sm.ref.uid).region, shift:sm.shift, o1s:XPS_UI._pr(g('O 1s')).comps.length, ba:XPS_UI._pr(g('Ba 4d')).comps.length, c1s:!!(XPS_UI._pr(g('C 1s')).result)}; })()")
    check('4 calibrate while switching spectrum: fits C 1s only, view stays where the user went', st4['refRegion'] == 'C 1s' and st4['c1s'] and st4['o1s'] == 0 and st4['ba'] == 0 and st4['sel'] == 'O 1s' and abs(st4['shift'] - 0.159) < 0.05, st4)

    # ---- 11 (v1.9.2): K2312 VGD 는 prop 415 = true (TF 를 측정 때 제거) → 'auto' 에서 TF 를 다시 나누지 않음 + 짧은 배지(설명은 title)
    J("(()=>{ const S=XPS_UI._state(); S.quant.tf={model:'auto'}; S.ui.tab='quant'; XPS_UI.rerender(); })()"); pg.wait_for_timeout(300)
    k11 = J("XPS_UI._state().specs.find(x=>x.format==='VGD').sample")
    q11 = J(f"XPS_UI._quantRows({json.dumps(k11)}).filter(r=>r.sp.format==='VGD').map(r=>({{el:r.el,tf:r.tf.model,rm:!!r.tf.removed,meta:r.sp.meta.tfRemoved}}))")
    b11 = J("(()=>{ const b=document.querySelector('[data-tfrm=\"panel\"]'); const t=[...document.querySelectorAll('#xpsPanel table.xq-tbl td.num')].map(e=>e.textContent); return b?{text:b.textContent, title:b.title, prof:!!document.querySelector('#xqProf')}:null; })()")
    J("(()=>{ const S=XPS_UI._state(); S.quant.tf={model:'power',x:1}; XPS_UI.rerender(); })()"); pg.wait_for_timeout(200)
    w11 = J("(()=>{ const b=document.querySelector('[data-tfrm=\"power\"]'); return b?b.textContent:null; })()")
    J("(()=>{ const S=XPS_UI._state(); S.quant.tf={model:'auto'}; S.ui.tab='fit'; XPS_UI.rerender(); })()"); pg.wait_for_timeout(200)
    check("11 VGD with TF removed at acquisition (prop 415): 'auto' transmission applies none (T = 1) for every VGD row; short badge, explanation in the tooltip; no 'Use VGD transmission' button; forcing (PE/KE)^x warns",
          len(q11) >= 4 and all(r['tf'] == 'none' and r['rm'] and r['meta'] for r in q11) and b11 and 'TF already removed at acquisition' in b11['text'] and len(b11['text']) < 50
          and 'prop 415' in b11['title'] and '9.8 at%' in b11['title'] and not b11['prof'] and w11 and 'twice' in w11, [q11[:2], b11 and b11['text'], w11])

    # ---- 12 (v1.9.2): Avantage 규약 Shirley 가 데이터 위로 지나는 경고 (avantage_a2_patch §4) — K2312 O 1s 전체 창: above 60 % → 강한 경고, Ta 4d: 경고 없음
    def avwarn(region, typ):
        J(f"(()=>{{ const S=XPS_UI._state(); const s=S.specs.filter(x=>x.format==='VGD'&&x.region==='{region}').pop(); S.sel=s.uid; const p=XPS_UI._pr(s); p.lo=s.be[0]; p.hi=s.be[s.n-1]; p.bg={{type:'{typ}',conv:'avantage',avgEV:0.5}}; p.ver=(p.ver||0)+1; S.ui.tab='fit'; XPS_UI.rerender(); }})()")
        pg.wait_for_timeout(300)
        return J("(()=>{ const e=document.querySelector('[data-avwarn]'); return e?{kind:e.dataset.avwarn, text:e.textContent, title:e.title}:null; })()")
    w12o, w12t, w12s = avwarn('O 1s', 'shirley'), avwarn('Ta 4d', 'shirley'), avwarn('O 1s', 'smart')
    check("12 Avantage-convention Shirley above the data over most of the window (K2312 O 1s): one-line strong warning 'use Smart', numbers in the tooltip; none for Ta 4d Shirley or O 1s Smart",
          w12o and w12o['kind'] == 'strong' and 'use Smart' in w12o['text'] and len(w12o['text']) < 70 and 'above the smoothed data on 60 %' in w12o['title'] and w12t is None and w12s is None, [w12o, w12t, w12s])
    pg.screenshot(path=os.path.join(out, 'f12_avantage_warn.png'), full_page=True)

    # ---- 7 · 9: 텍스트 파일 (KE · 쉼표 소수 · 배경 뺀 자료)
    pg.set_input_files('#xpsFin', [f_ke, f_cd, f_bs])
    pg.wait_for_function("['KE-1','CD-1','BS-1'].every(k=>XPS_UI._state().specs.some(s=>s.sample===k))", timeout=20000); pg.wait_for_timeout(300)
    allmsg = pg.inner_text('#xpsMsg')
    kes = J("(()=>{ const s=XPS_UI._state().specs.find(x=>x.sample==='KE-1'); return {uid:s.uid, lo:s.be[0], hi:s.be[s.n-1], ac:XPS_UI._axisCheck(s.uid)}; })()")
    check('7 headerless KE file is flagged (library C 1s not in BE range, KE fits)', kes['ac'] and kes['ac']['fits'] and kes['ac']['other'] == 'KE' and 'Read as KE' in allmsg, kes['ac'])
    J(f"(()=>{{ const S=XPS_UI._state(); S.sel='{kes['uid']}'; S.ui.tab='info'; XPS_UI.rerender(); }})()"); pg.wait_for_timeout(200)
    pg.select_option('#xiAx', 'KE'); pg.click('#xiAxOk'); pg.wait_for_timeout(300)
    ke2 = J("(()=>{ const s=XPS_UI._cur(); return {lo:s.be[0], hi:s.be[s.n-1], axis:s.meta.axis, ac:XPS_UI._axisCheck(s.uid), lo0:XPS_UI._pr().lo}; })()")
    check('7 Info tab BE/KE toggle re-reads the file as KE', ke2['axis'] == 'KE' and ke2['lo'] < 284.2 < ke2['hi'] and not ke2['ac'] and abs(ke2['lo0'] - ke2['lo']) < 1e-9, ke2)
    pg.screenshot(path=os.path.join(out, 'f07_info_axis.png'), full_page=True)
    J("XPS_UI._state().ui.tab='fit'; XPS_UI.rerender()")
    J('XPS_UI._auto()'); idle(); fit()
    kefit = J('XPS_UI._pr().result.summ.map(q=>q.be)')
    check('7 KE file fits at the right BE after re-reading', any(abs(x - 284.8) < 0.05 for x in kefit), kefit)
    cd = J("(()=>{ const s=XPS_UI._state().specs.find(x=>x.sample==='CD-1'); return {n:s.n, lo:s.be[0], hi:s.be[s.n-1], cd:s.meta.commaDecimal}; })()")
    check('9 comma-decimal CSV is read', cd['n'] == 321 and abs(cd['lo'] - 525) < 1e-6 and abs(cd['hi'] - 541) < 1e-6 and cd['cd'], cd)
    bs = J("(()=>{ const s=XPS_UI._state().specs.find(x=>x.sample==='BS-1'); const p=XPS_UI._pr(s); return {uid:s.uid, bgSub:s.meta.bgSub, bg:p.bg.type, w:p.weight, warn:s.warn.join(' / ')}; })()")
    check('9 background-subtracted input: warned, Background None + Constant weight', bs['bgSub'] and bs['bg'] == 'none' and bs['w'] == 'const' and 'None' in bs['warn'], bs['warn'][:90])
    J(f"(()=>{{ const S=XPS_UI._state(); S.sel='{bs['uid']}'; XPS_UI.rerender(); }})()"); pg.wait_for_timeout(150)
    J('XPS_UI._auto()'); idle(); fit()
    bsf = J('(()=>{ const r=XPS_UI._pr().result; return {chi2nu:r.stats.chi2nu, be:r.summ.map(q=>q.be), area:r.summ.map(q=>q.area)}; })()')
    check('9 background-subtracted fit (constant weight): sensible χ²ν and BE', 0.5 < bsf['chi2nu'] < 2 and any(abs(x - 531) < 0.05 for x in bsf['be']), bsf)

    # ---- 5: Auger 직접 추가 · 거부
    S9 = os.path.join(tmp, 's9'); os.makedirs(S9, exist_ok=True)
    dl, hv = 3.2, 1486.6

    def wr(name, lo, hi_, c, F, h, base):
        n = int(round((hi_ - lo) / 0.05)) + 1
        xs = [lo + i * 0.05 for i in range(n)]
        return write(name, xs, [base + h * gls(x, c, F) for x in xs], header='Binding Energy (eV)\tCounts / s')
    fa = wr('[S-9] Na1s.txt', 1060 + dl, 1082 + dl, 1071.5 + dl, 1.6, 20000, 3000)
    fb = wr('[S-9] Na KLL.txt', 485 + dl, 510 + dl, hv - 989.8 + dl, 2.4, 9000, 5000)
    pg.set_input_files('#xpsFin', [fa, fb])
    pg.wait_for_function("XPS_UI._state().specs.filter(s=>s.sample==='S-9').length>=2", timeout=20000)
    J("(()=>{ const S=XPS_UI._state(); const s=S.specs.find(x=>x.sample==='S-9'&&x.region==='Na 1s'); S.sel=s.uid; S.ui.tab='auger'; XPS_UI.rerender(); })()"); pg.wait_for_timeout(200)
    pg.click('#xaAdd'); pg.wait_for_timeout(200)
    e = J("(()=>{ const S=XPS_UI._state(); const e=S.auger.entries.filter(x=>x.sample==='S-9').pop(); const r=XPS_UI._augCalc(e); const R=u=>(S.specs.find(s=>s.uid===u)||{}).region; return {p:R(e.p.uid), a:R(e.a.uid), alpha:r.alpha, err:r.err}; })()")
    check('5 manual Auger entry defaults: photoelectron = Na 1s, Auger = Na KLL', e['p'] == 'Na 1s' and e['a'] == 'Na KLL' and e['alpha'] and abs(e['alpha'] - 2061.3) < 0.1, e)
    bad1 = J("(()=>{ const S=XPS_UI._state(); const e=JSON.parse(JSON.stringify(S.auger.entries.filter(x=>x.sample==='S-9').pop())); e.a=Object.assign({},e.p); return XPS_UI._augCalc(e); })()")
    check('5 Auger side = photoelectron region / same spectrum is refused (no α′ = hν)', bool(bad1.get('err')) and bad1.get('alpha') is None, (bad1.get('err') or '')[:80])
    pg.screenshot(path=os.path.join(out, 'f05_auger.png'), full_page=True)

    # ---- 9: survey SNR — 비우면 4, 저장
    select('Survey'); J("XPS_UI._state().ui.tab='survey'; XPS_UI.rerender()"); pg.wait_for_timeout(200)
    pg.fill('#xsSnr', ''); pg.dispatch_event('#xsSnr', 'change'); pg.wait_for_timeout(100)
    v_empty = pg.input_value('#xsSnr'); s_empty = J('XPS_UI._state().ui.survey[XPS_UI._cur().uid].snr')
    pg.fill('#xsSnr', '5.5'); pg.dispatch_event('#xsSnr', 'change'); pg.fill('#xsTol', '1.5'); pg.dispatch_event('#xsTol', 'change')
    pg.wait_for_timeout(1200)
    sv_uid = J('XPS_UI._cur().uid')
    pg.reload(); pg.wait_for_timeout(1000); pg.click('#nav button[data-v="xps"]'); pg.wait_for_timeout(1500)
    J(f"(()=>{{ const S=XPS_UI._state(); S.sel='{sv_uid}'; S.ui.tab='survey'; XPS_UI.rerender(); }})()"); pg.wait_for_timeout(300)
    check('9 survey SNR: empty → 4; SNR/tolerance persist across reload', v_empty == '4' and s_empty == 4 and pg.input_value('#xsSnr') == '5.5' and pg.input_value('#xsTol') == '1.5', [v_empty, s_empty, pg.input_value('#xsSnr'), pg.input_value('#xsTol')])

    # ---- 10: 390 px — at% 열이 보이고 그래프 글자가 읽을 만한 크기
    select('C 1s'); J("XPS_UI._state().ui.tab='quant'; XPS_UI.rerender()")
    pg.set_viewport_size({'width': 390, 'height': 844}); pg.wait_for_timeout(500)
    J('XPS_UI.rerender()'); pg.wait_for_timeout(300)
    th = J("(()=>{ const t=[...document.querySelectorAll('#xpsPanel table.xq-tbl th')].find(x=>x.textContent.trim()==='at %'); const r=t.getBoundingClientRect(); return {l:r.left, r:r.right, w:innerWidth}; })()")
    fsz = J("(()=>{ const t=[...document.querySelectorAll('#xpsSvg text')].find(x=>/^\\d/.test(x.textContent)); return t.getBoundingClientRect().height; })()")
    check('10 390 px: at% column visible without horizontal scroll', th['r'] <= th['w'] and th['l'] >= 0, th)
    check('10 390 px: plot tick labels ≥ 7 px on screen', fsz >= 7, round(fsz, 1))
    res['mobile_overflow'] = J('document.documentElement.scrollWidth>window.innerWidth+1')
    check('10 390 px: no page overflow', not res['mobile_overflow'])
    pg.screenshot(path=os.path.join(out, 'f10_mobile_quant.png'), full_page=True)
    pg.evaluate("window.scrollTo(0,0)"); pg.screenshot(path=os.path.join(out, 'f10_mobile_plot.png'))
    b.close()

check('no console errors', not errs, errs[:3])
print(json.dumps({'failed': bad, 'n': len(res)}, ensure_ascii=False))
sys.exit(1 if bad else 0)
