"""CasaXPS Cookbook 2026 선모양 (SPEC_shapes2.md) — 브라우저(Chromium) 실측, 항목마다 판정
사용: python3 test/ui_xps_shapes2.py <주입된 html> <그림 폴더>
1 성분 모양 선택에 새 모양(TLA · GL(m,n) · SGL(m,n) · TP) → TLA 로 바꿔 피팅(합성 TLA 자료, 참값 회복)
2 GL(m,n) 꼬리 선택 T/PHI → 매개변수 칸   3 수식어 ST/SB → 칸·기본 고정   4 CasaXPS 표기 입력(ST(1)LA(1,143), DS(0.3,600,40), TP(1,15,0)) · 오류 문구
5 칸 편집(μ) → 상태 반영   6 GL(30,50)T(2) · ST(1)LA(1,143) 로 피팅   7 390 px 가로 넘침 없음 · 콘솔 오류 없음
8 멈춤 규칙(SPEC §3.3b): 정지 이유가 Worker 결과에 실림 · 'stalled' 결과 → 한 줄 알림(매개변수 이름) + title 자세히 + 'Stalled' 칸"""
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


random.seed(5)
tmp = tempfile.mkdtemp()


def tla(x, F, a, mu):   # Cookbook p. 66 T(x_c: α, μ), x_c = −x/F (INNO x = BE − E0), unit area by a simple sum (fine grid)
    t = -x / F
    return (1 / (1 + 4 * t * t)) ** a * (math.pi / 2 - math.atan(2 * t) + math.pi / mu)


TR = dict(F=0.9, alpha=1.3, mu=3.0, be=118.0, area=60000.0)
grid = [-80 + 0.001 * i for i in range(160001)]
Z = sum(tla(x, TR['F'], TR['alpha'], TR['mu']) for x in grid) * 0.001   # (tail beyond ±80 eV < 1e-3 of the area, inside the noise)
be = [112 + 0.05 * i for i in range(241)]


def pois(l):
    return max(0.0, round(l + math.sqrt(max(l, 1e-9)) * random.gauss(0, 1)))


ys = [pois(400 + TR['area'] * tla(e - TR['be'], TR['F'], TR['alpha'], TR['mu']) / Z) for e in be]
fn = os.path.join(tmp, '[Alumina] Al2s.txt')
with open(fn, 'w') as fo:
    fo.write('Binding Energy (eV)\tCounts/s\n')
    for x, y in sorted(zip(be, ys), reverse=True):
        fo.write(f'{x:.3f}\t{y:.3f}\n')

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
    pg.set_input_files('#xpsFin', [fn])
    pg.wait_for_function('window.XPS_UI && XPS_UI._state().specs.length>=1', timeout=20000)
    J("(()=>{ const st=XPS_UI._state(); st.sel=st.specs[0].uid; st.ui.tab='fit'; XPS_UI.rerender(); })()"); pg.wait_for_timeout(200)

    def idle(t=180000):
        pg.wait_for_function("document.querySelector('#xpsProg').hidden", timeout=t)
        pg.wait_for_timeout(200)

    def fit():
        J("(()=>{ window.__prog=window.__prog||[]; if(!window.__progObs){ const el=document.querySelector('#xpsProg .t'); window.__progObs=new MutationObserver(()=>window.__prog.push(el.textContent)); window.__progObs.observe(el,{childList:true,characterData:true,subtree:true}); } })()")
        pg.click('#xpsFit')
        pg.wait_for_function('XPS_UI._pr().result && !XPS_UI._pr().result.stale', timeout=180000)
        idle()

    cur = lambda: J('XPS_UI._pr().comps[0]')
    cid = lambda: J('XPS_UI._pr().comps[0].id')
    msg = lambda: pg.inner_text('#xpsMsg') if pg.query_selector('#xpsMsg') else ''

    # 한 성분 모형으로 시작 (GLS) — 선택해서 카드를 연다
    pg.select_option('#xpsBg', 'linear'); pg.wait_for_timeout(150)
    J("(()=>{ const p=XPS_UI._pr(); p.comps=[]; p.cons={}; XPS_UI.rerender(); })()")
    J("XPS_UI._state().ui.addShape='GLS'")
    J("(()=>{ const p=XPS_UI._pr(); const c=XPS_CORE.newComp('GLS',118.1,{F:1.1,area:50000}); if(c){ p.comps.push(c); XPS_UI._state().ui.selComp=c.id; p.ver=(p.ver||0)+1; } XPS_UI.rerender(); })()")
    if not J('XPS_UI._pr().comps.length'):
        pg.evaluate("XPS_UI._auto()"); pg.wait_for_timeout(300)
        J("(()=>{ const p=XPS_UI._pr(); p.comps=p.comps.slice(0,1); XPS_UI._state().ui.selComp=p.comps[0].id; XPS_UI.rerender(); })()")
    pg.wait_for_timeout(200)

    # ---- 1 모양 선택에 새 모양 → TLA 로 피팅
    opts = J(f"Array.from(document.querySelectorAll('[data-shape=\"{cid()}\"] option')).map(o=>o.value)")
    addo = J("Array.from(document.querySelectorAll('#xpsAddS option')).map(o=>o.value)")
    pg.select_option(f'[data-shape="{cid()}"]', 'TLA'); pg.wait_for_timeout(250)
    c1 = cur(); has_mu = pg.query_selector(f'.comp[data-c="{c1["id"]}"] input[data-p="mu"]') is not None
    # 배경 = 성분 TP(1,1,0) (Cookbook p. 46 의 쓰임: 상수 바탕을 모형 성분으로), 영역 배경 None
    pg.select_option('#xpsBg', 'none'); pg.wait_for_timeout(150)
    J("(()=>{ const p=XPS_UI._pr(); const r=XPS_CORE.casaToComp('TP(1,1,0)',{be:118,area:9000,fwhm:40}); r.comp.name='Bkg'; p.comps.push(r.comp); p.cons[r.comp.id]={F:{mode:'fixed'},be:{mode:'fixed'}}; p.ver=(p.ver||0)+1; XPS_UI.rerender(); })()")
    fit()
    r1 = J('(()=>{ const r=XPS_UI._pr().result; const c=r.comps[0]; return {chi2nu:r.stats.chi2nu, mu:c.mu, alpha:c.alpha, F:c.F, be:c.be, area:c.area, legend:document.querySelector("#xpsComps").innerText}; })()')
    check('1 shape selectors list TLA · GLN · SGLN · TP; TLA card has μ; TLA fit on synthetic TLA data recovers α, μ, F, area (±10 %)',
          all(k in opts for k in ('TLA', 'GLN', 'SGLN', 'TP')) and all(k in addo for k in ('TLA', 'GLN', 'SGLN', 'TP')) and has_mu and r1['chi2nu'] < 1.4
          and abs(r1['mu'] / TR['mu'] - 1) < 0.1 and abs(r1['alpha'] / TR['alpha'] - 1) < 0.1 and abs(r1['F'] / TR['F'] - 1) < 0.1 and abs(r1['area'] / TR['area'] - 1) < 0.1 and abs(r1['be'] - TR['be']) < 0.02,
          {k: round(v, 4) for k, v in r1.items() if k != 'legend'})
    J("(()=>{ const p=XPS_UI._pr(); const tp=p.comps.find(c=>c.shape==='TP'); window.__tp=tp&&{area:tp.area,F:tp.F}; p.comps=p.comps.filter(c=>c.shape!=='TP'); p.ver++; XPS_UI.rerender(); })()")
    tp1 = J('window.__tp')
    check('1b TP(1,1,0) background component: area / width ≈ constant level 400 CPS', tp1 and abs(tp1['area'] / tp1['F'] / 400 - 1) < 0.03, tp1)
    pg.select_option('#xpsBg', 'linear'); pg.wait_for_timeout(150)
    pg.screenshot(path=os.path.join(out, 's2_01_tla_fit.png'), full_page=True)

    # ---- 2 GL(m,n) + 꼬리
    pg.select_option(f'[data-shape="{cid()}"]', 'GLN'); pg.wait_for_timeout(200)
    pg.select_option(f'[data-tail="{cid()}"]', 'T'); pg.wait_for_timeout(200)
    t2 = cur(); has_tk = pg.query_selector(f'input[data-p="tk"]') is not None
    pg.select_option(f'[data-tail="{cid()}"]', 'PHI'); pg.wait_for_timeout(200)
    t2b = cur(); has_tj = pg.query_selector(f'input[data-p="tj"]') is not None
    check('2 GL(m,n): tail select T(k) → k field; PHI(j,k) → j, k fields; Casa line shows the notation',
          t2['shape'] == 'GLN' and t2.get('tail') == 'T' and has_tk and t2b.get('tail') == 'PHI' and has_tj and 'PHI(' in pg.inner_text('#xpsComps'), [t2.get('tail'), t2b.get('tail'), t2b.get('tk'), t2b.get('tj')])

    # ---- 3 수식어 ST / SB
    pg.select_option(f'[data-shape="{cid()}"]', 'LA'); pg.wait_for_timeout(200)
    pg.select_option(f'[data-mod="{cid()}"]', 'ST'); pg.wait_for_timeout(200)
    c3 = cur(); stg_btn = J("(()=>{ const b=document.querySelector('[data-cm=\"stG\"]'); return b?b.textContent:null; })()")
    has_stmu = pg.query_selector('input[data-p="stMu"]') is not None
    pg.select_option(f'[data-mod="{cid()}"]', 'SB'); pg.wait_for_timeout(200)
    c3b = cur(); no_stmu = pg.query_selector('input[data-p="stMu"]') is None
    check('3 modifier ST(μ,γ) adds μ and γ fields (γ fixed by default = ST(μ)); SB(0) removes them',
          c3.get('mod') == 'ST' and has_stmu and stg_btn == 'Fixed' and c3b.get('mod') == 'SB' and no_stmu, [c3.get('mod'), stg_btn, c3b.get('mod')])
    pg.select_option(f'[data-mod="{cid()}"]', ''); pg.wait_for_timeout(150)

    # ---- 4 CasaXPS 표기 입력
    def casa(s):
        sel = f'[data-casa="{cid()}"]'; pg.fill(sel, s); pg.press(sel, 'Enter'); pg.wait_for_timeout(250)
    casa('ST(1)LA(1,143)'); c4 = cur()
    wg_dis = J("(()=>{ const i=document.querySelector('input[data-p=\"wG\"]'); return i?i.disabled:null; })()")
    nval = J(f"document.querySelector('[data-cn=\"{cid()}\"]').value")
    casa('DS(0.3,600,40)'); c4b = cur(); lval = J(f"document.querySelector('[data-dsl=\"{cid()}\"]').value")
    casa('TP(1,15,0)'); c4c = cur(); ta_btn = J("(()=>{ const b=document.querySelector('[data-cm=\"ta\"]'); return b?b.textContent:null; })()")
    before = cur(); casa('XYZ(1)'); m4 = msg(); after = cur()
    check("4 CasaXPS string input: ST(1)LA(1,143) → LA + ST, n 143 (G width from n, field locked); DS(0.3,600,40) → DS, cut l 40; TP(1,15,0) → TP, α fixed; unknown string → clear message, component unchanged",
          c4['shape'] == 'LA' and c4.get('mod') == 'ST' and c4.get('casaN') == 143 and wg_dis is True and nval == '143'
          and c4b['shape'] == 'DS' and c4b.get('dsl') == 40 and lval == '40' and c4b.get('casaN') == 600
          and c4c['shape'] == 'TP' and c4c.get('tb') == 15 and ta_btn == 'Fixed'
          and 'not recognised' in m4 and before == after, [c4.get('casaN'), wg_dis, c4b.get('dsl'), ta_btn, m4[:80]])

    # ---- 5 칸 편집
    # ---- 4c CasaXPS 줄 (SPEC_casa_v3 §3): 꼭대기 매끈(TLA) → fwhm·asym 그대로(정함), 꺾인 꼭대기 LA(2,3,0) → '≈' (경험값)
    eqt = lambda: J(f"(()=>{{ const e=document.querySelector('[data-casaeq=\"{cid()}\"]'); return e?e.textContent:null; }})()")
    casa('LA(2,3,0)'); t4la = eqt(); casa('TLA(1.3,3,0)'); t4tla = eqt()
    check("4c CasaXPS line under the component: TLA(1.3,3,0) shows fwhm and asym without '≈' (pinned); LA(2,3,0) shows 'fwhm ≈ …, asym ≈ …' (empirical)",
          t4tla and 'TLA(1.3,3,0)' in t4tla and 'asym' in t4tla and '≈' not in t4tla and t4la and 'fwhm ≈' in t4la and 'asym ≈' in t4la, [t4tla, t4la])
    pg.fill('input[data-p="mu"]', '2.5'); pg.press('input[data-p="mu"]', 'Enter'); pg.wait_for_timeout(200)
    check('5 editing the μ field updates the component', abs(cur()['mu'] - 2.5) < 1e-12, cur()['mu'])

    # ---- 6 다른 새 모양으로 피팅
    casa('GL(30,50)T(2)'); fit(); r6 = J('(()=>{ const r=XPS_UI._pr().result; return {chi2nu:r.stats.chi2nu, shape:r.comps[0].shape, n:r.comps[0].casaN, tk:r.comps[0].tk, fwhm:r.summ[0].fwhm}; })()')
    casa('ST(1)LA(1,143)'); fit(); r6b = J('(()=>{ const r=XPS_UI._pr().result; return {chi2nu:r.stats.chi2nu, mod:r.comps[0].mod, stMu:r.comps[0].stMu, fwhm:r.summ[0].fwhm}; })()')
    tbl = pg.inner_text('#xpsComps')
    prog = J("(window.__prog||[]).filter(t=>t.indexOf('반복')>=0).slice(-3)")
    wk = J("!!document.getElementById('xps-core-src')")
    check('6b fits run in the Web Worker (engine source injected) and report LM iterations in the progress bar', wk and len(prog) > 0, prog)
    check('6 fits with GL(30,50)T(2) and ST(1)LA(1,143) run (finite χ²ν, true FWHM reported); ST legend in the card',
          r6['shape'] == 'GLN' and r6['chi2nu'] > 0 and r6['fwhm'] and r6b['mod'] == 'ST' and r6b['chi2nu'] > 0 and r6b['fwhm'] and 'ST(' in tbl, [r6, r6b])
    pg.screenshot(path=os.path.join(out, 's2_06_st_fit.png'), full_page=True)

    # ---- 8 멈춤 규칙 (SPEC §3.3b): 피팅 결과에 정지 이유가 실려 오고(Worker), 'stalled' 면 한 줄 알림 + 자세한 것은 title
    st8 = J('(()=>{ const r=XPS_UI._pr().result; return {stop:r.stop, sstop:r.stats.stop}; })()')
    J("(()=>{ const p=XPS_UI._pr(), r=p.result, c=r.comps[0]; r.stop='stalled'; r.stall={gain:0.05,phi:1,flags:[{id:c.id,ci:0,k:'F',p:c.F,lo:0.05,hi:10,why:['flat']},{id:'bg',ci:-1,k:'c0',p:1,lo:-1e9,hi:1e9,why:['degenerate']}]}; XPS_UI.rerender(); })()")
    pg.wait_for_timeout(250)
    n8 = J("(()=>{ const ns=Array.from(document.querySelectorAll('.note.warn')).filter(e=>e.textContent.indexOf('진전 멈춤')>=0); const s=Array.from(document.querySelectorAll('.strip > div')).find(e=>e.textContent.indexOf('Stalled')===0); return {n:ns.length, text:ns[0]&&ns[0].textContent, lines:ns[0]&&Math.round(ns[0].getBoundingClientRect().height/parseFloat(getComputedStyle(ns[0]).lineHeight||'16')), title:ns[0]&&ns[0].title, strip:!!s, stripWarn:!!(s&&s.classList.contains('warn'))}; })()")
    check("8 stop reason reaches the UI from the Worker (result.stop / stats.stop = 'converged' after a normal fit); a 'stalled' result shows one short Korean notice naming the ill-determined parameters, details in the tooltip, and a 'Stalled' box in the strip",
          st8['stop'] == 'converged' and st8['sstop'] == 'converged' and n8['n'] == 1 and '진전 멈춤' in (n8['text'] or '') and 'FWHM' in (n8['text'] or '') and '배경' in (n8['text'] or '')
          and '구별되지 않음' in (n8['title'] or '') and 'SPEC §3.3b' in (n8['title'] or '') and n8['strip'] and n8['stripWarn'] and len(n8['text'] or '') < 120, [st8, n8])
    J("(()=>{ const n=Array.from(document.querySelectorAll('.note.warn')).find(e=>e.textContent.indexOf('진전 멈춤')>=0); if(n) n.scrollIntoView({block:'center'}); })()"); pg.wait_for_timeout(150)
    pg.screenshot(path=os.path.join(out, 's2_08_stalled.png'), full_page=False)
    J("(()=>{ const r=XPS_UI._pr().result; r.stop='converged'; r.stall=null; XPS_UI.rerender(); })()")

    # ---- 7 390 px
    pg.set_viewport_size({'width': 390, 'height': 844}); pg.wait_for_timeout(400)
    J("XPS_UI.rerender()"); pg.wait_for_timeout(300)
    ov = J("document.documentElement.scrollWidth-document.documentElement.clientWidth")
    row = J("(()=>{ const r=document.querySelector('.xps-s2'); if(!r) return null; const b=r.getBoundingClientRect(); return {w:b.width, right:b.right, vw:innerWidth, sw:r.scrollWidth, cw:r.clientWidth}; })()")
    check('7 390 px: no horizontal page overflow; shape-options row fits', ov <= 1 and row and row['right'] <= row['vw'] + 1 and row['sw'] <= row['cw'] + 1, [ov, row])
    pg.screenshot(path=os.path.join(out, 's2_07_390.png'), full_page=True)
    b.close()

check('no console errors', not errs, errs[:3])
json.dump(res, open(os.path.join(out, 'ui_xps_shapes2.json'), 'w'), ensure_ascii=False, indent=1)
print(f"\n{len(res) - len(bad)} / {len(res)} ok")
sys.exit(1 if bad else 0)
