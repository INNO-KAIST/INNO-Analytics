"""XPS 탭 브라우저 실측 (Playwright · Chromium)
사용: python3 test/ui_xps.py <주입된 html> <그림 폴더>
파일 열기 → 대전 보정(C 1s) → 자동 제안 → 피팅(Worker) → 다중 시작 → Monte Carlo → 결과·정량·서베이·정보 탭,
콘솔 오류 0 개를 확인한다."""
import sys, os, glob, time, json
from playwright.sync_api import sync_playwright

html, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
files = sorted(glob.glob(os.path.join(HERE, 'data', 'bct600', '*.xls'))) + sorted(glob.glob(os.path.join(HERE, 'data', 'k2312', '*.VGD')))
errs = []
res = {}
with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={'width': 1500, 'height': 950})
    # 외부 글꼴(Google Fonts) 요청 실패는 이 격리 환경의 네트워크 문제라 제외
    pg.on('console', lambda m: errs.append('console.' + m.type + ': ' + m.text) if m.type == 'error' and 'ERR_TUNNEL' not in m.text and 'ERR_NAME' not in m.text else None)
    pg.on('pageerror', lambda e: errs.append('pageerror: ' + str(e)))
    pg.goto('file://' + os.path.abspath(html))
    pg.wait_for_timeout(800)
    pg.click('#nav button[data-v="xps"]')
    pg.wait_for_selector('.xps', timeout=10000)
    pg.set_input_files('#xpsFin', files)
    pg.wait_for_function('window.XPS_UI && XPS_UI._state().specs.length>=12', timeout=20000)
    res['n_specs'] = pg.evaluate('XPS_UI._state().specs.length')
    res['formats'] = pg.evaluate('[...new Set(XPS_UI._state().specs.map(s=>s.format))]')
    pg.screenshot(path=os.path.join(out, '01_list.png'))

    def select(region, fmt):
        pg.evaluate("([r,f])=>{ const s=XPS_UI._state().specs.find(x=>x.region===r&&x.format===f); XPS_UI._state().sel=s.uid; XPS_UI.rerender(); }", [region, fmt])

    # C 1s: 자동 제안 → 피팅
    select('C 1s', 'Avantage XLS')
    pg.wait_for_timeout(300)
    pg.evaluate('XPS_UI._auto()')
    res['c1s_auto_n'] = pg.evaluate('XPS_UI._pr().comps.length')
    pg.evaluate('XPS_UI._fit("single")')
    pg.wait_for_function('XPS_UI._pr().result && !XPS_UI._pr().result.stale', timeout=60000)
    res['c1s_chi2nu'] = pg.evaluate('XPS_UI._pr().result.stats.chi2nu')
    res['c1s_be'] = pg.evaluate('XPS_UI._pr().result.summ.map(q=>+q.be.toFixed(3))')
    # 대전 보정: 가장 큰 C 1s 성분을 284.8 로
    pg.evaluate('XPS_UI._shiftC1s()')
    pg.wait_for_timeout(300)
    res['shift'] = pg.evaluate('(()=>{ const sp=XPS_UI._cur(); return XPS_UI._state().samples[sp.sample].shift; })()')
    pg.screenshot(path=os.path.join(out, '02_c1s_fit.png'))
    # 다중 시작
    pg.evaluate('XPS_UI._fit("multi")')
    pg.wait_for_function('XPS_UI._pr().result && XPS_UI._pr().result.multi && !XPS_UI._pr().result.stale', timeout=120000)
    res['c1s_multi'] = pg.evaluate('(()=>{ const m=XPS_UI._pr().result.multi; return {N:m.N,frac:m.frac}; })()')
    # Monte Carlo
    pg.evaluate('XPS_UI._state().ui.mcN=60')
    pg.evaluate('XPS_UI._mc()')
    pg.wait_for_function('XPS_UI._pr().result.mc', timeout=180000)
    res['c1s_mc_n'] = pg.evaluate('XPS_UI._pr().result.mc.n')
    pg.wait_for_timeout(400)
    pg.screenshot(path=os.path.join(out, '03_c1s_result.png'), full_page=True)

    # Co 2p (doublet 자동 제안) → 피팅
    select('Co 2p', 'Avantage XLS')
    pg.wait_for_timeout(300)
    pg.evaluate('XPS_UI._auto()')
    res['co2p_doublets'] = pg.evaluate('XPS_UI._pr().comps.filter(c=>c.doublet).map(c=>c.doublet)')
    pg.evaluate('XPS_UI._fit("single")')
    pg.wait_for_function('XPS_UI._pr().result && !XPS_UI._pr().result.stale', timeout=120000)
    res['co2p_chi2nu'] = pg.evaluate('XPS_UI._pr().result.stats.chi2nu')
    pg.evaluate("XPS_UI._state().ui.tab='fit'; XPS_UI.rerender()")
    pg.wait_for_timeout(300)
    pg.screenshot(path=os.path.join(out, '04_co2p.png'))

    # 정량 탭
    pg.evaluate("XPS_UI._state().ui.tab='quant'; XPS_UI.rerender()")
    pg.wait_for_timeout(500)
    pg.screenshot(path=os.path.join(out, '05_quant.png'), full_page=True)
    res['quant_text'] = pg.inner_text('#xpsPanel')[:600]

    # 서베이
    select('Survey', 'Avantage XLS')
    pg.evaluate("XPS_UI._state().ui.tab='survey'; XPS_UI.rerender()")
    pg.wait_for_timeout(300)
    pg.click('#xsRun')
    pg.wait_for_timeout(800)
    res['survey_acc'] = pg.evaluate('(()=>{ const sp=XPS_UI._cur(); const st=XPS_UI._state().ui.survey[sp.uid]; return Object.keys(st.acc); })()')
    pg.screenshot(path=os.path.join(out, '06_survey.png'), full_page=True)

    # VGD 정보 탭
    select('C 1s', 'VGD')
    pg.evaluate("XPS_UI._state().ui.tab='info'; XPS_UI.rerender()")
    pg.wait_for_timeout(400)
    res['info_text'] = pg.inner_text('#xpsPanel')[:700]
    pg.screenshot(path=os.path.join(out, '07_vgd_info.png'), full_page=True)

    # 겹쳐 보기 (O 1s 두 시료)
    pg.evaluate("(()=>{ const S=XPS_UI._state(); S.ui.overlay=S.specs.filter(s=>s.region==='O 1s').map(s=>s.uid); XPS_UI.rerender(); })()")
    pg.wait_for_timeout(400)
    pg.screenshot(path=os.path.join(out, '08_overlay.png'))

    # 첨부 시료와 무관한 원소: 합성 Fe 2p(이름으로 인식) · 이름 없는 Au 4f(정보 탭에서 원소 지정)
    import numpy as np, tempfile
    tmp = tempfile.mkdtemp()
    def gls(x, F, m):
        return (1 - m) * (2 / F) * np.sqrt(np.log(2) / np.pi) * np.exp(-4 * np.log(2) * x * x / F / F) + m * (2 / (np.pi * F)) / (1 + 4 * x * x / F / F)
    def make(path, be0, dE, r, lo, hi, A=20000.0):
        E = np.round(np.arange(lo, hi + 1e-9, 0.05), 3); P = A / (1 + r) * gls(E - be0, 1.2, 0.3) + A * r / (1 + r) * gls(E - be0 - dE, 1.2, 0.3)
        Q = np.concatenate([[0], np.cumsum(0.5 * (P[1:] + P[:-1]) * 0.05)]); mu = 3000 + 0.06 * Q + P
        y = np.random.default_rng(3).poisson(mu * 0.3) / 0.3
        open(path, 'w').write(''.join(f'{e:.2f}\t\t{v:.1f}\n' for e, v in zip(E[::-1], y[::-1])))
    make(os.path.join(tmp, '[Steel-A] Fe2p.txt'), 706.8, 13.1, 0.5, 700, 730)
    make(os.path.join(tmp, '[Steel-A] Peak3.txt'), 84.0, 3.67, 0.75, 78, 95)
    pg.set_input_files('#xpsFin', [os.path.join(tmp, '[Steel-A] Fe2p.txt'), os.path.join(tmp, '[Steel-A] Peak3.txt')])
    pg.wait_for_function('XPS_UI._state().specs.length>=14', timeout=20000)
    select('Fe 2p', '텍스트'); pg.wait_for_timeout(200)
    pg.evaluate('XPS_UI._auto()')
    res['fe_suggest'] = pg.evaluate('(()=>{ const p=XPS_UI._pr(); return p.comps.map(c=>({be:+c.be.toFixed(2),d:c.doublet,dEc:(p.cons[c.id]||{}).dE})); })()')
    pg.evaluate('XPS_UI._fit("single")')
    pg.wait_for_function('XPS_UI._pr().result && !XPS_UI._pr().result.stale', timeout=60000)
    res['fe_fit'] = pg.evaluate('(()=>{ const r=XPS_UI._pr().result; return {chi2nu:+r.stats.chi2nu.toFixed(3), be:r.summ.map(q=>+q.be.toFixed(3)), area:r.summ.map(q=>Math.round(q.area)), dE:r.comps.map(c=>c.doublet&&+c.doublet.dE.toFixed(3))}; })()')
    # 이름에서 원소를 못 읽은 스펙트럼 → 정보 탭에서 Au 4f 지정
    pg.evaluate("(()=>{ const s=XPS_UI._state().specs.find(x=>x.file&&x.file.indexOf('Peak3')>=0); XPS_UI._state().sel=s.uid; XPS_UI._state().ui.tab='info'; XPS_UI.rerender(); })()")
    pg.wait_for_timeout(300)
    res['unknown_reg'] = pg.evaluate('XPS_UI._cur().reg')
    pg.fill('#xiEl', 'au'); pg.select_option('#xiOrb', '4f'); pg.click('#xiSet'); pg.wait_for_timeout(300)
    res['assigned_reg'] = pg.evaluate('XPS_UI._cur().reg.label')
    pg.evaluate('XPS_UI._auto()')
    pg.evaluate('XPS_UI._fit("single")')
    pg.wait_for_function('XPS_UI._pr().result && !XPS_UI._pr().result.stale', timeout=60000)
    res['au_fit'] = pg.evaluate('(()=>{ const r=XPS_UI._pr().result; return {chi2nu:+r.stats.chi2nu.toFixed(3), be:r.summ.map(q=>+q.be.toFixed(3)), dE:r.comps.map(c=>c.doublet&&+c.doublet.dE.toFixed(3))}; })()')
    # 대전 보정: 선택 성분(Au 4f7/2)을 기준 BE 84.0 으로 → Δ ≈ 0
    pg.evaluate("(()=>{ const p=XPS_UI._pr(); XPS_UI._state().ui.selComp=p.comps[0].id; const sp=XPS_UI._cur(); XPS_UI._state().samples[sp.sample]=Object.assign(XPS_UI._state().samples[sp.sample]||{}, {target:84.0}); XPS_UI.rerender(); })()")
    pg.evaluate("XPS_UI._state().ui.chgOpen=true; XPS_UI.rerender()"); pg.wait_for_timeout(200)
    pg.click('#xcSel'); pg.wait_for_timeout(300)
    res['au_shift'] = pg.evaluate("XPS_UI._state().samples[XPS_UI._cur().sample].shift")
    pg.evaluate("XPS_UI._state().ui.tab='quant'; XPS_UI.rerender()"); pg.wait_for_timeout(400)
    res['steel_quant'] = pg.inner_text('#xpsPanel')[:400]
    pg.screenshot(path=os.path.join(out, '10_generic_quant.png'), full_page=True)

    # 대전 보정 한 번에: K2312 VGD 6 개(한 지점) — Ba 4d 를 보는 중에 'C 1s 로 맞추기' → C 1s 자동 제안·피팅 → Δ → Ba 4d 로 돌아옴
    pg.evaluate("XPS_UI._state().ui.tab='fit'; XPS_UI.rerender()")
    import glob
    nv0 = pg.evaluate("XPS_UI._state().specs.length")
    pg.set_input_files('#xpsFin', sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'k2312', '*.VGD'))))
    pg.wait_for_timeout(1500)
    res['dup_skipped'] = pg.evaluate("XPS_UI._state().specs.length") == nv0
    # 처음 상태로: 이 시료의 보정·C 1s 피팅을 지운다
    pg.evaluate("(()=>{ const S=XPS_UI._state(); const v=S.specs.filter(s=>s.format==='VGD'); const k=v[0].sample; S.samples[k]={shift:0,target:284.8,refNote:''}; v.forEach(s=>{ const p=XPS_UI._pr(s); p.comps=[]; p.cons={}; p.result=null; }); XPS_UI.rerender(); })()")
    res['vgd_samples'] = pg.evaluate("[...new Set(XPS_UI._state().specs.filter(s=>s.format==='VGD').map(s=>s.sample))]")
    pg.evaluate("(()=>{ const S=XPS_UI._state(); const s=S.specs.filter(x=>x.format==='VGD'&&x.region==='Ba 4d').pop(); S.sel=s.uid; XPS_UI.rerender(); })()")
    pg.wait_for_timeout(200); res['chg_before'] = pg.inner_text('#xpsChg .chip')
    pg.click('#xpsChgGo')
    pg.wait_for_function("(()=>{ const sp=XPS_UI._cur(); const sm=XPS_UI._state().samples[sp.sample]; return sm&&sm.ref; })()", timeout=60000); pg.wait_for_timeout(300)
    res['chg_after'] = {'sel': pg.evaluate('XPS_UI._cur().region'), 'shift': pg.evaluate('XPS_UI._state().samples[XPS_UI._cur().sample].shift'), 'chip': pg.inner_text('#xpsChg .chip')}
    pg.screenshot(path=os.path.join(out, '11_charge.png'), full_page=False)

    # 다중항 준위: K2312 Co 3p 자동 제안 → doublet 없음 + 안내
    pg.evaluate("(()=>{ const S=XPS_UI._state(); const s=S.specs.filter(x=>x.format==='VGD'&&x.region==='Co 3p').pop(); S.sel=s.uid; S.ui.tab='fit'; XPS_UI.rerender(); })()")
    pg.evaluate('XPS_UI._auto()'); pg.wait_for_timeout(300)
    res['co3p_mult'] = {'doublets': pg.evaluate('XPS_UI._pr().comps.filter(c=>c.doublet).length'), 'n': pg.evaluate('XPS_UI._pr().comps.length'), 'banner': '다중항 준위' in pg.inner_text('#xpsPanel')}

    # Auger parameter + 화학량론: 합성 텍스트 파일(시료 S-9, Na 1s BE 1071.50+δ, Na KLL KE 989.80 → α′ = 2061.30, 대전 δ = 3.20)
    import math
    tmpd = os.path.join(out, 'synth'); os.makedirs(tmpd, exist_ok=True)
    def gls(x, c, F, m=0.3):
        s = F / (2 * math.sqrt(2 * math.log(2))); gg = math.exp(-0.5 * ((x - c) / s) ** 2); ll = 1 / (1 + ((x - c) / (F / 2)) ** 2); return (1 - m) * gg + m * ll
    dlt = 3.2; hv = 1486.6
    def wr(name, lo, hi, c, F, h, base):
        with open(os.path.join(tmpd, name), 'w') as fo:
            fo.write('Binding Energy (eV)\tCounts / s\n')
            n = int(round((hi - lo) / 0.05)) + 1
            for i in range(n):
                x = lo + i * 0.05; fo.write(f'{x:.3f}\t{base + h * gls(x, c, F):.3f}\n')
    wr('[S-9] Na1s.txt', 1060 + dlt, 1082 + dlt, 1071.5 + dlt, 1.6, 20000, 3000)
    wr('[S-9] Na KLL.txt', 485 + dlt, 510 + dlt, hv - 989.8 + dlt, 2.4, 9000, 5000)
    wr('[S-9] F1s.txt', 675 + dlt, 695 + dlt, 684.6 + dlt, 1.7, 30000, 4000)
    pg.set_input_files('#xpsFin', [os.path.join(tmpd, x) for x in ['[S-9] Na1s.txt', '[S-9] Na KLL.txt', '[S-9] F1s.txt']])
    pg.wait_for_function("XPS_UI._state().specs.some(s=>s.sample==='S-9')", timeout=20000)
    pg.evaluate("(()=>{ const S=XPS_UI._state(); const s=S.specs.find(x=>x.sample==='S-9'&&x.region==='Na 1s'); S.sel=s.uid; S.ui.tab='auger'; XPS_UI.rerender(); })()")
    pg.click('#xaAuto'); pg.wait_for_timeout(300)
    res['auger'] = pg.evaluate("(()=>{ const A=XPS_UI._state().auger; return A.entries.filter(e=>e.sample==='S-9').map(e=>({el:e.el})); })()")
    res['auger_text'] = pg.inner_text('.aug-card .h')
    pg.fill('#xrEl', 'Na'); pg.fill('#xrLb', 'test ref'); pg.fill('#xrBE', '1072.0'); pg.fill('#xrKE', '989.0'); pg.click('#xrAdd'); pg.wait_for_timeout(200)
    res['wagner_pts'] = pg.evaluate("document.querySelectorAll('#xpsPanel svg.wag circle, #xpsPanel svg.wag rect[width=\"8\"]').length")
    pg.screenshot(path=os.path.join(out, '12_auger.png'), full_page=True)
    pg.evaluate("XPS_UI._state().ui.tab='quant'; XPS_UI.rerender()"); pg.wait_for_timeout(300)
    pg.fill('#xstF', 'NaF'); pg.click('#xstOk'); pg.wait_for_timeout(300)
    res['stoich'] = pg.evaluate("[...document.querySelectorAll('#xpsPanel table.xps-tbl')].pop().innerText")[:300]
    pg.fill('#xstF', 'NaXq2'); pg.click('#xstOk'); pg.wait_for_timeout(200)
    res['stoich_err'] = '읽지 못했습니다' in pg.inner_text('#xpsPanel')
    pg.fill('#xstF', 'NaF'); pg.click('#xstOk'); pg.wait_for_timeout(200)
    pg.screenshot(path=os.path.join(out, '13_stoich.png'), full_page=True)

    # 모형 민감도: BCT C 1s(xls, 피팅됨) — 결과 탭에서 실행
    pg.evaluate("(()=>{ const S=XPS_UI._state(); const s=S.specs.find(x=>x.format==='Avantage XLS'&&x.region==='C 1s'); S.sel=s.uid; S.ui.tab='fit'; XPS_UI.rerender(); })()")
    pg.evaluate('XPS_UI._fit("single")'); pg.wait_for_function('XPS_UI._pr().result && !XPS_UI._pr().result.stale', timeout=60000)
    pg.evaluate("XPS_UI._state().ui.tab='result'; XPS_UI.rerender()"); pg.wait_for_timeout(300)
    pg.click('#xsnRun'); pg.wait_for_function('XPS_UI._pr().sens && XPS_UI._pr().sens.sum', timeout=120000); pg.wait_for_timeout(300)
    res['sens'] = pg.evaluate("(()=>{ const s=XPS_UI._pr().sens; return {n:s.sum.n, ok:s.sum.nOk, comps:s.sum.comps.map(c=>[+c.areaPct.min.toFixed(1),+c.areaPct.max.toFixed(1)])}; })()")
    res['sens_rows'] = pg.evaluate("document.querySelectorAll('#xpsPanel table.sens tbody tr').length")
    pg.screenshot(path=os.path.join(out, '14_sens.png'), full_page=True)

    # 다시 읽기(IndexedDB 보관) 확인
    pg.reload(); pg.wait_for_timeout(1000)
    pg.click('#nav button[data-v="xps"]'); pg.wait_for_timeout(1500)
    res['after_reload_specs'] = pg.evaluate('XPS_UI._state().specs.length')

    # 휴대폰 폭
    pg.set_viewport_size({'width': 390, 'height': 844}); pg.wait_for_timeout(500)
    pg.screenshot(path=os.path.join(out, '09_mobile.png'), full_page=False)
    res['mobile_overflow'] = pg.evaluate('document.documentElement.scrollWidth>window.innerWidth+1')
    b.close()

res['errors'] = errs
print(json.dumps(res, ensure_ascii=False, indent=1))
sys.exit(1 if errs else 0)
