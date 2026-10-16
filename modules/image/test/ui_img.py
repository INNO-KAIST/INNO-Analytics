"""Image Analysis 탭 브라우저 실측 (Playwright · Chromium)
사용: python3 test/ui_img.py <주입된 html> <그림 폴더>
파일 열기(여러 형식·옆 파일) → Worker 로 필터(ImageJ 결과와 화소 대조) → Threshold → Analyze Particles(ImageJ 와 개수·면적 대조)
→ 마우스로 ROI·Measure·Line profile → Find Maxima → Macro 기록·batch·.ijm → 되돌리기 → 새로 고침 후 복원 → 모바일, 콘솔 오류 0"""
import sys, os, json, numpy as np
from playwright.sync_api import sync_playwright
html, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
H = os.path.dirname(os.path.abspath(__file__)); REF = os.path.join(H, '..', 'ref', 'data'); IOD = os.path.join(H, 'io')
man = json.load(open(os.path.join(REF, 'manifest.json')))
ref = lambda n: np.fromfile(os.path.join(REF, 'out', n + '.f32'), np.float32)
files = [os.path.join(REF, 'in', f) for f in ['g8.tif', 'bin8.tif', 'stk.tif', 'g16.tif']] + [os.path.join(IOD, f) for f in ['tfei.tif', 'afm.ibw', 'scan.sxm', 'trgb.tif', 'thitachi.tif', 'thitachi.txt', 'p16.png', 'tij_hyper.tif']]
errs = []; res = {}; fails = []
def chk(c, name, info=''):
    (res.setdefault('checks', [])).append([bool(c), name, str(info)])
    if not c: fails.append(name + ' — ' + str(info))
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={'width': 1560, 'height': 1000})
    pg.on('console', lambda m: errs.append('console.' + m.type + ': ' + m.text) if m.type == 'error' and 'ERR_TUNNEL' not in m.text and 'ERR_NAME' not in m.text else None)
    pg.on('pageerror', lambda e: errs.append('pageerror: ' + str(e)))
    pg.goto('file://' + os.path.abspath(html)); pg.wait_for_timeout(800)
    pg.click('#nav button[data-v="imgx"]'); pg.wait_for_selector('.imx', timeout=10000)
    pg.screenshot(path=os.path.join(out, '00_empty.png'))
    pg.set_input_files('#imFin', files)
    pg.wait_for_function('window.IMG_UI && IMG_UI._state().images.length>=16', timeout=30000)
    ims = pg.evaluate('IMG_UI._state().images.map(i=>({n:i.name,w:i.w,h:i.h,bit:i.bit,rgb:i.rgb,nc:i.nc,nz:i.nz,nt:i.nt,cal:i.cal?[+i.cal.pw.toPrecision(6),i.cal.unit]:null}))')
    res['images'] = ims
    by = {i['n']: i for i in ims}
    chk(by['tfei.tif']['cal'] == [2.5, 'nm'], 'FEI 태그 보정 2.5 nm', by['tfei.tif']['cal'])
    chk(by['thitachi.tif']['cal'] and abs(by['thitachi.tif']['cal'][0] - 1.98438) < 1e-4, 'Hitachi 옆 파일 보정', by['thitachi.tif']['cal'])
    chk(by['tij_hyper.tif']['nc'] == 2 and by['tij_hyper.tif']['nz'] == 2 and by['tij_hyper.tif']['nt'] == 3, 'ImageJ hyperstack 차원', by['tij_hyper.tif'])
    chk(sum(1 for i in ims if i['n'].startswith('afm.ibw')) == 3 and sum(1 for i in ims if i['n'].startswith('scan.sxm')) == 4, 'IBW 3 채널 · SXM 4 영상')
    chk(by['trgb.tif']['rgb'], 'RGB TIFF')
    pg.wait_for_timeout(400); pg.screenshot(path=os.path.join(out, '01_loaded.png'))
    sel = lambda name: pg.evaluate("n=>{ const S=IMG_UI._state(); S.sel=S.images.find(i=>i.name===n).id; IMG_UI.rerender(); }", name)
    getp = lambda: np.array(pg.evaluate("(()=>{ const im=IMG_UI._cur(); return Array.from(IMG_UI._px.get(im.id)[0]); })()"), np.float32)
    # 1) Worker 필터 = ImageJ
    for op, prm, rn in [('gaussian', {'sigma': 2}, 'gauss_g8_2'), ('median', {'radius': 3.5}, 'median_g8_3.5'), ('subtractBg', {'radius': 20}, 'rb_g8_20'), ('bandpass', {'large': 40, 'small': 3, 'stripes': 'none', 'tolerance': 5, 'autoscale': True, 'saturate': True}, 'bp_g8')]:
        sel('g8.tif'); pg.evaluate("([o,p])=>IMG_UI._run(o,p)", [op, prm]); pg.wait_for_function('!document.querySelector("#imProg")||document.querySelector("#imProg").hidden', timeout=60000); pg.wait_for_timeout(200)
        a = getp(); r = ref(rn); eq = float(np.mean(a == r)); mx = float(np.max(np.abs(a - r)))
        chk(eq >= 0.999 and mx <= 1, f'Worker {op} = ImageJ {rn}', f'일치 {eq*100:.3f} %, 최대차 {mx}')
        pg.evaluate('document.querySelector("#imUndo").click()'); pg.wait_for_timeout(150)
    a = getp(); chk(np.array_equal(a, np.fromfile(os.path.join(REF, 'in', 'g8.f32'), np.float32)), '되돌리기 4 번 → 원본')
    # 2) Process 패널 + Preview 화면
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.panel='process'; S.ui.proc='unsharp'; S.ui.preview=true; IMG_UI.rerender(); })()"); pg.wait_for_timeout(900)
    res['preview_on'] = pg.evaluate('document.querySelector("#imStatus").textContent.includes("Preview")'); chk(res['preview_on'], 'Preview 표시')
    pg.screenshot(path=os.path.join(out, '02_process_preview.png'))
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.preview=false; S.ui.panel='threshold'; S.ui.thr.method='Otsu'; S.ui.thr.dark=true; S.ui.thr.manual=false; IMG_UI.rerender(); })()"); pg.wait_for_timeout(700)
    thr = pg.evaluate('[IMG_UI._state().ui.thr.lower,IMG_UI._state().ui.thr.upper]'); rt = [c for c in man['cases'] if c['name'] == 'thr_g8'][0]['thr']['Otsu|dark']
    chk(thr == rt, 'Threshold 패널 Otsu dark = ImageJ', f'{thr} vs {rt}')
    pg.click('#imTall'); pg.wait_for_selector('.im-try .t', timeout=20000); pg.wait_for_timeout(300)
    pg.screenshot(path=os.path.join(out, '03_threshold_tryall.png'))
    pg.click('#imTapply'); pg.wait_for_function('IMG_UI._state().images.some(i=>i.name==="Mask of g8.tif")', timeout=20000)
    m = getp(); r = ref('bin_g8_mask'); chk(np.array_equal(m, r), 'Convert to Mask = ImageJ', f'{np.mean(m==r)*100:.3f} %')
    # 3) Analyze Particles on bin8 = ImageJ 'pa'
    sel('bin8.tif'); pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.panel='analyze'; Object.assign(S.ui.pa,{minSize:0,maxSize:'',unitPx:true,minCirc:0,maxCirc:1,exclude:false,holes:false,redirect:'',show:'outlines'}); IMG_UI.rerender(); })()")
    pg.wait_for_timeout(300); pg.click('#paRun'); pg.wait_for_function('(()=>{ const im=IMG_UI._cur(); return IMG_UI._state().particles[im.id]; })()', timeout=30000)
    pa = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return IMG_UI._state().particles[im.id].rows.map(r=>[r.Area,r["Perim."],r.Feret]); })()')
    rows = [c for c in man['cases'] if c['name'] == 'pa'][0]['rows']
    ok_pa = len(pa) == len(rows) and all(abs(a[0] - r['Area']) < 1e-6 and abs(a[1] - r['Perim.']) < 1e-6 and abs(a[2] - r['Feret']) < 1e-6 for a, r in zip(pa, rows))
    chk(ok_pa, 'Analyze Particles (UI) = ImageJ', f'{len(pa)} vs {len(rows)}')
    pg.wait_for_timeout(300); pg.screenshot(path=os.path.join(out, '04_particles.png'), full_page=True)
    pg.evaluate("(()=>{ IMG_UI._state().ui.bottom='dist'; IMG_UI.rerender(); })()"); pg.wait_for_timeout(300)
    chk(pg.evaluate('!!document.querySelector(".im-dist svg .fit")'), 'Size distribution + log-normal 곡선')
    pg.screenshot(path=os.path.join(out, '05_distribution.png'), full_page=True)
    # 4) 마우스 ROI → Measure, Line profile
    sel('g8.tif'); pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.tool='rect'; S.ui.panel='measure'; IMG_UI.rerender(); })()"); pg.wait_for_timeout(300)
    box = pg.locator('#imOverlay').bounding_box(); x0, y0 = box['x'] + box['width'] * 0.35, box['y'] + box['height'] * 0.35
    pg.mouse.move(x0, y0); pg.mouse.down(); pg.mouse.move(x0 + 120, y0 + 80, steps=6); pg.mouse.up(); pg.wait_for_timeout(200)
    roi = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return IMG_UI._state().roiSel[im.id]; })()'); chk(roi and roi['type'] == 'rect', '마우스로 사각형 ROI', roi)
    n0 = pg.evaluate('IMG_UI._state().results.length'); pg.keyboard.press('Control+m'); pg.wait_for_timeout(300)
    rr = pg.evaluate('IMG_UI._state().results.slice(-1)[0]'); chk(pg.evaluate('IMG_UI._state().results.length') == n0 + 1, 'Ctrl+M → Results 한 줄')
    # 같은 사각형을 core 로 직접 재어 비교
    b_ = pg.evaluate('(()=>{ const im=IMG_UI._cur(); const r=IMG_UI._state().roiSel[im.id]; const xs=r.pts.map(p=>p[0]), ys=r.pts.map(p=>p[1]); return [Math.min(...xs),Math.min(...ys),Math.max(...xs)-Math.min(...xs),Math.max(...ys)-Math.min(...ys)]; })()')
    g8 = np.fromfile(os.path.join(REF, 'in', 'g8.f32'), np.float32).reshape(180, 240); sub = g8[b_[1]:b_[1] + b_[3], b_[0]:b_[0] + b_[2]]
    chk(abs(rr['Mean'] - float(sub.astype(np.float64).mean())) < 1e-6 and rr['Area'] == sub.size, 'ROI Mean·Area = numpy', f"{rr['Mean']} vs {sub.astype(np.float64).mean()}")
    pg.keyboard.press('t'); chk(pg.evaluate('(()=>{ const im=IMG_UI._cur(); return IMG_UI._state().rois[im.id].length; })()') == 1, 'T → ROI Manager')
    pg.evaluate("(()=>{ IMG_UI._state().ui.tool='line'; IMG_UI.rerender(); })()")
    pg.mouse.move(x0, y0 + 150); pg.mouse.down(); pg.mouse.move(x0 + 200, y0 + 160, steps=5); pg.mouse.up(); pg.wait_for_timeout(300)
    chk(pg.evaluate('!document.querySelector("#imProfile").hidden && !!document.querySelector("#imProfile svg")'), 'Line profile 그래프')
    pg.screenshot(path=os.path.join(out, '06_roi_profile.png'), full_page=True)
    # 5) Find Maxima = ImageJ
    pg.evaluate("IMG_UI._run('findMaxima',{prominence:10,light:false,exclude:false,strict:false,output:'points'})"); pg.wait_for_timeout(800)
    npt = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return IMG_UI._state().roiSel[im.id].pts.length; })()'); nref = len([c for c in man['cases'] if c['name'] == 'max_g8_10_list'][0]['rows'])
    chk(npt == nref, 'Find Maxima 개수 = ImageJ', f'{npt} vs {nref}')
    # 6) Macro: 기록 → 다른 영상(g16)에 batch → .ijm
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.macro.steps=[]; S.macro.rec=true; })()")
    sel('g8.tif'); pg.evaluate("IMG_UI._run('gaussian',{sigma:2})"); pg.wait_for_timeout(600)
    pg.evaluate("IMG_UI._run('threshold',{method:'Otsu',dark:true,manual:false,lower:0,upper:255,stack:false,newImage:true})"); pg.wait_for_timeout(800)
    pg.evaluate("IMG_UI._run('analyze',{minSize:0,maxSize:'',unitPx:true,minCirc:0,maxCirc:1,exclude:true,holes:false,redirect:'__parent',show:'outlines',thr:null})"); pg.wait_for_timeout(1200)
    steps = pg.evaluate('IMG_UI._steps().map(s=>s.op)'); chk(steps == ['gaussian', 'threshold', 'analyze'], '매크로 자동 기록', steps)
    ns0 = pg.evaluate('IMG_UI._state().summary.length')
    pg.evaluate("(async()=>{ const S=IMG_UI._state(); const im=S.images.find(i=>i.name==='g16.tif'); await IMG_UI._runSteps(S.macro.steps,im); })()"); pg.wait_for_timeout(2500)
    ns1 = pg.evaluate('IMG_UI._state().summary.length'); last = pg.evaluate('IMG_UI._state().summary.slice(-1)[0]')
    chk(ns1 == ns0 + 1 and last['Image'].startswith('Mask of g16'), '매크로를 다른 영상에 적용 → Summary', last)
    ijm = pg.evaluate('IMG_UI._toIJM()'); open(os.path.join(out, 'macro.ijm'), 'w').write(ijm)
    chk('Gaussian Blur' in ijm and 'Analyze Particles' in ijm, '.ijm 내보내기')
    pg.evaluate("(()=>{ IMG_UI._state().ui.panel='macro'; IMG_UI.rerender(); })()"); pg.wait_for_timeout(300); pg.screenshot(path=os.path.join(out, '07_macro.png'), full_page=True)
    # 7) 다른 패널 화면
    for k in ['adjust', 'binary', 'info']:
        sel('afm.ibw · HeightRetrace' if k == 'adjust' else 'bin8.tif'); pg.evaluate(f"(()=>{{ IMG_UI._state().ui.panel='{k}'; IMG_UI.rerender(); }})()"); pg.wait_for_timeout(400)
        pg.screenshot(path=os.path.join(out, f'08_{k}.png'))
    sel('tij_hyper.tif'); pg.wait_for_timeout(300); pg.screenshot(path=os.path.join(out, '09_hyperstack.png'))
    chk(pg.evaluate('!document.querySelector("#imStack").hidden'), '스택 막대 (c/z/t)')
    # 7b) 검토에서 나온 결함 회귀
    sel('g8.tif'); pg.evaluate("(()=>{ const S=IMG_UI._state(); S.macro.steps=[]; })()"); ops0 = pg.evaluate('IMG_UI._cur().ops.length')
    pg.evaluate("IMG_UI._run('gaussian',{sigma:1.5})"); pg.wait_for_timeout(500)
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.panel='adjust'; IMG_UI.rerender(); })()"); pg.wait_for_timeout(200)
    pg.fill('#imScPx', '100'); pg.fill('#imScD', '50'); pg.select_option('#imScU', 'nm'); pg.click('#imScSet'); pg.wait_for_timeout(200)
    pg.evaluate('document.querySelector("#imUndo").click()'); pg.wait_for_timeout(300)
    st = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return {cal:im.cal&&im.cal.pw, steps:IMG_UI._steps().map(s=>s.op), ops:im.ops.length}; })()')
    chk(st['cal'] == 0.5 and st['steps'] == ['setScale'] and st['ops'] == ops0, '되돌리기: 보정 유지·해당 매크로 단계만 지움', st)
    # Preview 계산 중에 Apply
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.panel='process'; S.ui.proc='median'; S.ui.preview=true; IMG_UI.rerender(); })()"); pg.wait_for_timeout(260)
    v0 = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return im.ops.length; })()')
    pg.evaluate("document.querySelector('#imApply').click()"); pg.wait_for_function('!document.querySelector("#imProg")||document.querySelector("#imProg").hidden', timeout=30000); pg.wait_for_timeout(500)
    v1 = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return {ops:im.ops.length, pv:IMG_UI._state().ui.preview, status:document.querySelector("#imStatus").textContent}; })()')
    chk(v1['ops'] == v0 + 1 and not v1['pv'] and 'Preview' not in v1['status'], 'Preview 중 Apply → 적용되고 Preview 꺼짐', v1)
    # Z Project: hyperstack 은 t 를 유지, RGB 합은 3 채널
    sel('tij_hyper.tif'); pg.evaluate("IMG_UI._run('zproject',{method:'max',float:false})"); pg.wait_for_timeout(800)
    zp = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return [im.name,im.nc,im.nz,im.nt,IMG_UI._px.get(im.id).length]; })()')
    chk(zp[1:] == [2, 1, 3, 6], 'Z Project (hyperstack) → 채널·시간 유지', zp)
    sel('trgb.tif'); pg.evaluate("IMG_UI._run('duplicate',{slice:false})"); pg.wait_for_timeout(300)
    # 3 장 RGB 스택을 만들어 합
    pg.evaluate("(()=>{ const S=IMG_UI._state(); const im=IMG_UI._cur(); const P=IMG_UI._px.get(im.id); IMG_UI._px.set(im.id,P.concat(P.map(p=>Float32Array.from(p)),P.map(p=>Float32Array.from(p)))); im.nz=3; im.ver++; IMG_UI.rerender(); })()")
    pg.evaluate("IMG_UI._run('zproject',{method:'sum',float:false})"); pg.wait_for_timeout(800)
    zr = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return [im.rgb,im.nc,im.bit,IMG_UI._px.get(im.id).length]; })()')
    chk(zr == [False, 3, 32, 3], 'Z Project RGB Sum → R·G·B 3 채널 32-bit', zr)
    # 7c) 현미경 형식 nd2·czi (lossless zlib·zstd1 — 브라우저 해독), 위치 많은 파일 → 고르기 창
    n0 = pg.evaluate('IMG_UI._state().images.length')
    pg.set_input_files('#imFin', [os.path.join(IOD, f) for f in ['n_t3z4c2.nd2', 'n_lossless.nd2', 'c_zstd1.czi', 'c_mosaic.czi', 'n_p8z4.nd2']])
    pg.wait_for_function('!document.querySelector("#imModal").hidden && document.querySelector("#ipOk")', timeout=30000)
    pg.screenshot(path=os.path.join(out, '11_import_choose.png'))
    pg.evaluate("document.querySelectorAll('#imModal .ipP').forEach((x,i)=>{ x.checked=i<2; x.dispatchEvent(new Event('input',{bubbles:true})); })")
    pg.fill('#ipZs', '2'); pg.click('#ipOk')
    pg.wait_for_function(f'IMG_UI._state().images.length>={n0}+6', timeout=30000); pg.wait_for_timeout(300)
    TR = json.load(open(os.path.join(IOD, 'truth.json')))
    def same_as_truth(name, key):
        a = np.array(pg.evaluate("n=>{ const S=IMG_UI._state(); const im=S.images.find(i=>i.name===n); return [].concat(...IMG_UI._px.get(im.id).map(p=>Array.from(p))); }", name), np.float32)
        t = np.fromfile(os.path.join(IOD, key + '.truth'), np.float32); return a.shape == t.shape and np.array_equal(a, t)
    mi = pg.evaluate("IMG_UI._state().images.filter(i=>/\\.(nd2|czi)/.test(i.name)).map(i=>({n:i.name,nc:i.nc,nz:i.nz,nt:i.nt,ch:i.chNames,lut:i.disp.map(d=>d.lut),z:i.zStep,cal:i.cal&&[+i.cal.pw.toPrecision(6),i.cal.unit]}))")
    by2 = {i['n']: i for i in mi}; res['micro'] = mi
    m1 = by2.get('n_t3z4c2.nd2', {})
    chk(m1.get('nc') == 2 and m1.get('nz') == 4 and m1.get('nt') == 3 and m1.get('ch') == ['FITC', 'TRITC'] and m1.get('z') == 0.5 and m1.get('cal') == [0.25, 'µm'], 'ND2: 차원·채널 이름·z step·보정', m1)
    chk(same_as_truth('n_t3z4c2.nd2', 'n_t3z4c2.nd2') and same_as_truth('n_lossless.nd2', 'n_lossless.nd2'), 'ND2 화소 = nd2 패키지 (비압축·lossless)')
    chk(same_as_truth('c_zstd1.czi', 'c_zstd1.czi') and same_as_truth('c_mosaic.czi', 'c_mosaic.czi'), 'CZI 화소 = libCZI/czifile (zstd1·mosaic 합성)')
    p8 = [i for i in mi if i['n'].startswith('n_p8z4.nd2')]
    chk(len(p8) == 2 and all(i['nz'] == 2 and abs(i['z'] - 0.5) < 1e-12 for i in p8), '고르기 창: 위치 2 개 · Z 2 장마다 → z step 0.5 µm', p8)
    # 7c-2) 고르기 창이 떠 있는 동안 다른 파일을 또 열면 → 차례로 (앞 파일이 사라지지 않음)
    n1 = pg.evaluate('IMG_UI._state().images.length')
    pg.set_input_files('#imFin', [os.path.join(IOD, 'n_p8z4.nd2')]); pg.wait_for_function('!document.querySelector("#imModal").hidden && document.querySelector("#ipOk")', timeout=20000)
    pg.set_input_files('#imFin', [os.path.join(IOD, 'c_f32.czi')]); pg.wait_for_timeout(600); pg.click('#ipOk')
    pg.wait_for_function(f'IMG_UI._state().images.length>={n1}+9', timeout=30000); pg.wait_for_timeout(300)
    nn = pg.evaluate('IMG_UI._state().images.slice(%d).map(i=>i.name)' % n1)
    chk(sum(1 for x in nn if x.startswith('n_p8z4')) == 8 and 'c_f32.czi' in nn, '고르기 창 중 또 열기 → 두 파일 모두 열림', nn)
    # 7c-3) AFM 형식 (Nanoscope·GWY·ARDF) — 채널 수·값 단위·보정·화소 = 정답(Gwyddion·gwyfile·magic-afm)
    AF = os.path.join(IOD, 'afm'); TA = json.load(open(os.path.join(AF, 'truth', 'truth.json'))); n2 = pg.evaluate('IMG_UI._state().images.length')
    pg.set_input_files('#imFin', [os.path.join(AF, f) for f in ['ns_v531.spm', 'g_multi.gwy', 'a_fmap.ardf']])
    pg.wait_for_function(f'IMG_UI._state().images.length>={n2}+12', timeout=30000); pg.wait_for_timeout(300)
    am = pg.evaluate('IMG_UI._state().images.slice(%d).map(i=>({n:i.name,w:i.w,h:i.h,nz:i.nz,bit:i.bit,u:i.valueUnit,cal:i.cal&&[i.cal.pw,i.cal.ph,i.cal.unit],z:i.zStep||null}))' % n2); res['afm'] = am
    chk([x['n'].split(' · ')[0] for x in am] == ['ns_v531.spm'] * 3 + ['g_multi.gwy'] * 6 + ['a_fmap.ardf'] * 3, 'AFM: 채널마다 영상 (Nanoscope 3 · GWY 4+mask+volume · ARDF 3)', [x['n'] for x in am])
    chk([x['u'] for x in am] == ['nm', 'deg', 'V', 'nm', '', 'V', 'MPa', 'N/m', 'nm', 'nm', 'nN', 'MPa'], 'AFM: 값 단위 (길이 nm, Pa → MPa, N → nN)', [x['u'] for x in am])
    bk = am[8]; chk(bk['nz'] == 5 and abs(bk['z'] - 0.5) < 1e-9, 'GWY volume → z-stack, z step 0.5 µm', bk)
    def afm_same(k, fname, j):
        ti = TA[fname]['images'][j]; im_id = pg.evaluate('k=>IMG_UI._state().images[k].id', n2 + k)
        a = np.array(pg.evaluate("id=>[].concat(...IMG_UI._px.get(id).map(p=>Array.from(p)))", im_id), np.float64)
        f = pg.evaluate('k=>{ const im=IMG_UI._state().images[k]; return im.siScale||null; }', n2 + k)
        t = np.fromfile(os.path.join(AF, 'truth', ti['file']), np.float64)
        e = t * (f if ti['kind'] != 'mask' else 1)   # float32 로 한 번 반올림한 차이만 허용 (상대 1.5e-7, test_afm.js 와 같음)
        return f is not None and a.shape == t.shape and np.allclose(a, e, rtol=1.5e-7, atol=0, equal_nan=True)
    chk(afm_same(0, 'ns_v531.spm', 0) and afm_same(1, 'ns_v531.spm', 1) and afm_same(3, 'g_multi.gwy', 0) and afm_same(4, 'g_multi.gwy', 1) and afm_same(8, 'g_multi.gwy', 5) and afm_same(9, 'a_fmap.ardf', 0) and afm_same(11, 'a_fmap.ardf', 2),
        'AFM 화소 = 정답 × 표시 배율 (float32 반올림 안)')
    chk(abs(am[0]['cal'][0] - 2000 / 32) < 1e-9 and am[0]['cal'][2] == 'nm', 'Nanoscope 보정: 2 ~m / 32 = 62.5 nm/px', am[0]['cal'])
    pg.screenshot(path=os.path.join(out, '11b_afm.png'))
    # 값 단위: 같은 물리량을 남기는 연산(Gaussian)은 nm 유지, 8-bit 눈금 변환은 지움, 되돌리면 되살아남, Duplicate 는 이어받음
    vu = lambda: pg.evaluate('(()=>{ const im=IMG_UI._cur(); return [im.valueUnit, im.siScale]; })()')
    pg.evaluate('k=>{ IMG_UI._state().sel=IMG_UI._state().images[k].id; }', n2)
    pg.evaluate("IMG_UI._run('gaussian',{sigma:1})"); pg.wait_for_timeout(500); u1 = vu()
    pg.evaluate("IMG_UI._run('type',{to:'8',scale:true})"); pg.wait_for_timeout(500); u2 = vu()
    pg.evaluate("IMG_UI._undo()"); pg.wait_for_timeout(300); u3 = vu()
    pg.evaluate("IMG_UI._run('duplicate',{})"); pg.wait_for_timeout(400); u4 = vu()
    chk(u1 == ['nm', 1e9] and u2 == ['', None] and u3 == ['nm', 1e9] and u4 == ['nm', 1e9], '값 단위: Gaussian 유지 · 8-bit 변환 지움 · 되돌리기 복원 · Duplicate 이어받음', [u1, u2, u3, u4])
    # 7c-4) Color Threshold (RGB) — 패널 조작 → mask = 핵심 함수(ImageJ 와 2²⁴ 색 대조됨), Sample·Auto, 매크로 .ijm 을 ImageJ 로 실행한 결과와 같음
    sel('trgb.tif'); pg.evaluate("IMG_UI._state().ui.panel='threshold'; IMG_UI.rerender()"); pg.wait_for_selector('#ctSp', timeout=10000)
    pg.wait_for_selector('#ctApply:not([disabled])', timeout=20000)
    pg.select_option('#ctSp', 'HSB'); pg.wait_for_selector('#ctApply:not([disabled])', timeout=20000)
    for k, lo, hi in [(0, 30, 170), (1, 40, 255), (2, 20, 240)]:
        pg.fill(f'.ctMinN[data-k="{k}"]', str(lo)); pg.dispatch_event(f'.ctMinN[data-k="{k}"]', 'change')
        pg.fill(f'.ctMaxN[data-k="{k}"]', str(hi)); pg.dispatch_event(f'.ctMaxN[data-k="{k}"]', 'change')
    pg.click('.ctPass[data-k="0"]')   # Hue stop
    pg.wait_for_timeout(300); pg.screenshot(path=os.path.join(out, '11c_color_threshold.png'))
    cs = pg.evaluate('JSON.parse(JSON.stringify(IMG_UI._state().ui.cthr))')
    chk(cs['space'] == 'HSB' and cs['min'] == [30, 40, 20] and cs['max'] == [170, 255, 240] and cs['pass'] == [False, True, True], 'Color Threshold 패널 입력 → 상태', cs)
    nb = pg.evaluate('IMG_UI._state().images.length'); pg.click('#ctApply'); pg.wait_for_function(f'IMG_UI._state().images.length>{nb}', timeout=20000); pg.wait_for_timeout(200)
    mk = pg.evaluate("""(()=>{ const S=IMG_UI._state(); const m=S.images[S.images.length-1]; const src=S.images.find(i=>i.name==='trgb.tif'); const P=IMG_UI._px.get(src.id);
        const C=S.ui.cthr; const IC=self.IMG_CORE; const ref=IC.colorThreshold(P[0],P[1],P[2],src.w,src.h,{space:C.space,min:C.min,max:C.max,pass:C.pass}).data; const a=IMG_UI._px.get(m.id)[0];
        let bad=0,n=0; for(let i=0;i<a.length;i++){ if(a[i]!==ref[i]) bad++; if(a[i]) n++; } return {name:m.name,bit:m.bit,bad,n,N:a.length}; })()""")
    chk(mk['name'] == 'Mask of trgb.tif' and mk['bit'] == 8 and mk['bad'] == 0 and 0 < mk['n'] < mk['N'], 'Color Threshold → 8-bit mask = 핵심 함수', mk)
    # ImageJ 로 .ijm 실행 → 같은 mask
    import subprocess, tempfile
    ijm = pg.evaluate("(()=>{ const C=IMG_UI._state().ui.cthr; return IMG_UI._ops.colorThreshold.ij({space:C.space,min:C.min,max:C.max,pass:C.pass,newImage:true}); })()")
    td = tempfile.mkdtemp(); mf = os.path.join(td, 'ct.ijm'); rf = os.path.join(td, 'ct.raw')
    open(mf, 'w').write(f'setBatchMode(true); setOption("BlackBackground", true); open("{os.path.join(IOD, "trgb.tif")}");\ni=7;\n' + ijm + f'\nprint("I "+i+" "+getTitle()); saveAs("Raw Data", "{rf}");\n')
    lg = subprocess.run(f'JAVA_TOOL_OPTIONS= xvfb-run -a java -Xmx1g -cp /tmp/ij.jar ij.ImageJ -batch {mf}', shell=True, capture_output=True, text=True).stdout
    same = os.path.exists(rf) and np.array_equal(np.fromfile(rf, np.uint8).astype(np.float32), np.array(pg.evaluate("(()=>{ const S=IMG_UI._state(); return Array.from(IMG_UI._px.get(S.images[S.images.length-1].id)[0]); })()"), np.float32))
    chk(same and 'I 7 Mask of trgb.tif' in lg, 'Color Threshold .ijm 을 ImageJ 로 실행 → 같은 mask, 바깥 변수 i 그대로, 이름 Mask of …', lg.strip()[-200:])
    # 스택·hyperstack: '모든 장' + 새 영상 / 제자리 — 앱 mask = 핵심 함수, .ijm(장마다 고리) 을 ImageJ 로 실행한 결과와 같은 화소·차원
    for fn in ['rgb_stack.tif', 'rgb_hyper.tif']:
        for newimg in [True, False]:
            nb = pg.evaluate('IMG_UI._state().images.length'); pg.set_input_files('#imFin', [os.path.join(IOD, fn)]); pg.wait_for_function(f'IMG_UI._state().images.length>{nb}', timeout=20000); pg.wait_for_timeout(200)
            prm = "{space:'Lab',min:[30,100,0],max:[220,170,140],pass:[true,true,false],stack:true,newImage:%s}" % ('true' if newimg else 'false')
            info = pg.evaluate("""(async()=>{ const S=IMG_UI._state(); const src=S.images[S.images.length-1]; const P=IMG_UI._px.get(src.id).slice(); const IC=self.IMG_CORE; const nz=src.nz, nt=src.nt;
                const ref=[]; for(let t=0;t<nt;t++) for(let z=0;z<nz;z++){ const k=3*(z+nz*t); ref.push(Array.from(IC.colorThreshold(P[k],P[k+1],P[k+2],src.w,src.h,%s).data)); }
                await IMG_UI._run('colorThreshold',%s); const m=S.images.find(i=>i.id===S.sel); const got=IMG_UI._px.get(m.id).map(a=>Array.from(a));
                return {same:JSON.stringify(got)===JSON.stringify(ref),nz:m.nz,nt:m.nt,rgb:m.rgb,name:m.name,ijm:IMG_UI._ops.colorThreshold.ij(%s)}; })()""" % (prm, prm, prm))
            rf = os.path.join(td, 'cts.raw'); mf = os.path.join(td, 'cts.ijm')
            if os.path.exists(rf): os.remove(rf)
            open(mf, 'w').write(f'setBatchMode(true); setOption("BlackBackground", true); open("{os.path.join(IOD, fn)}");\ni=7;\n' + info['ijm'] + f'\ngetDimensions(w,h,c,z,t); print("I "+i+" "+getTitle()+" "+z+" "+t+" "+bitDepth()); saveAs("Raw Data", "{rf}");\n')
            lg = subprocess.run(f'JAVA_TOOL_OPTIONS= xvfb-run -a java -Xmx1g -cp /tmp/ij.jar ij.ImageJ -batch {mf}', shell=True, capture_output=True, text=True).stdout
            got = np.array(pg.evaluate("(()=>{ const S=IMG_UI._state(); return [].concat(...IMG_UI._px.get(S.sel).map(a=>Array.from(a))); })()"), np.float32)
            ijraw = np.fromfile(rf, np.uint8).astype(np.float32) if os.path.exists(rf) else None
            want = f"I 7 {'Mask of ' if newimg else ''}{fn} {info['nz']} {info['nt']} 8"
            chk(info['same'] and ijraw is not None and np.array_equal(ijraw, got) and want in lg, f'Color Threshold 스택 {fn} {"새 영상" if newimg else "제자리"}: 앱 = 핵심, .ijm(장마다) = ImageJ 화소·차원', [info['same'], want, lg.strip()[-160:]])
    # 7c-5) Analyze Skeleton — 앱 결과 = Fiji AnalyzeSkeleton (같은 골격을 SkelDriver 로), .ijm (Skeletonize + Analyze Skeleton) 을 Fiji 플러그인을 넣은 ImageJ 로 실행해 같은 표
    if os.environ.get('SKIP_SKEL') != '1' and os.path.exists('/tmp/askel_out/SkelDriver.class') and os.path.exists('/tmp/ijplug/plugins/AnalyzeSkeleton_.jar'):
        sel('bin8.tif'); n0s = pg.evaluate('IMG_UI._state().skelTrees.length')
        pg.evaluate("IMG_UI._run('skeleton',{skel:true,lsp:true,tagged:true,lspImg:true})"); pg.wait_for_function(f'IMG_UI._state().skelTrees.length>{n0s}', timeout=30000); pg.wait_for_timeout(300)
        sk = pg.evaluate("""(()=>{ const S=IMG_UI._state(); const tg=S.images.find(i=>i.name==='Tagged skeleton of bin8.tif'); const lp=S.images.find(i=>i.name==='Longest shortest paths of bin8.tif');
            return {w:tg.w,h:tg.h,tag:Array.from(IMG_UI._px.get(tg.id)[0]),lsp:Array.from(IMG_UI._px.get(lp.id)[0]),trees:S.skelTrees.filter(r=>r.Image==='bin8.tif'),br:S.skelBranches.filter(r=>r.Image==='bin8.tif'),sel:S.images.find(i=>i.id===S.sel).name,lut:tg.disp[0].lut}; })()""")
        import tifffile, subprocess, tempfile
        tdk = tempfile.mkdtemp(); skf = os.path.join(tdk, 'sk.tif'); tag = np.array(sk['tag'], np.float32).reshape(sk['h'], sk['w'])
        tifffile.imwrite(skf, ((tag > 0) * 255).astype(np.uint8))
        lg = subprocess.run(f'JAVA_TOOL_OPTIONS= java -cp /tmp/ij.jar:/tmp/askel_out SkelDriver {skf} 0 0 1 1 1 1', shell=True, capture_output=True, text=True).stdout.splitlines()
        T = [list(map(float, l.split()[2:])) for l in lg if l.startswith('T ')]; E = [list(map(float, l.split()[1:])) for l in lg if l.startswith('E ')]
        TAG = dict(tuple(map(int, x.split(':'))) for l in lg if l.startswith('TAG') for x in l.split()[1:])
        app_t = [[r['# Branches'], r['# Junctions'], r['# End-point voxels'], r['# Junction voxels'], r['# Slab voxels'], r['Average Branch Length'], r['# Triple points'], r['# Quadruple points'], r['Maximum Branch Length'], r['Longest Shortest Path'], r['spx'], r['spy'], r['spz']] for r in sk['trees']]
        same_t = len(app_t) == len(T) and all(np.allclose(a, b, rtol=1e-12, atol=1e-12) for a, b in zip(app_t, T))
        tagok = all(int(v) == TAG.get(i, 0) for i, v in enumerate(sk['tag']))
        # Branch information: 나무마다 길이 내림차순 (Fiji 표 순서) — Fiji 가지를 같은 규칙으로 정렬해 비교
        fe = []
        for t in range(len(T)):
            L = [e for e in E if int(e[0]) == t]; L = sorted(enumerate(L), key=lambda x: (-x[1][1], x[0])); fe += [x[1] for x in L]
        same_b = len(fe) == len(sk['br']) and all(abs(a[1] - b['Branch length']) < 1e-12 and [a[2], a[3], a[5], a[6]] == [b['V1 x'], b['V1 y'], b['V2 x'], b['V2 y']] for a, b in zip(fe, sk['br']))
        chk(same_t and same_b and tagok and sk['lut'] == 'Fire' and len(T) > 3, f'Analyze Skeleton (bin8 → Skeletonize): 골격 {len(T)} · 가지 {len(E)} — 표·Branch information·태그 = Fiji', [same_t, same_b, tagok, sk['lut']])
        lspN = sum(1 for v in sk['lsp'] if v == 96); chk(lspN > 0 and sum(1 for v in sk['lsp'] if v > 0) == sum(1 for v in sk['tag'] if v > 0), 'Longest shortest path 영상: 골격 255 위에 경로 96', lspN)
        ijm = pg.evaluate("IMG_UI._ops.skeleton.ij({skel:true,lsp:true,tagged:true})"); of = os.path.join(tdk, 'out.txt'); mf = os.path.join(tdk, 'm.ijm')
        open(mf, 'w').write(f'open("{os.path.join(REF, "in", "bin8.tif")}"); setOption("BlackBackground", true);\n' + ijm + f'\ns="N "+nResults; for(i=0;i<nResults;i++) s=s+"\\nR "+getResult("# Branches",i)+" "+d2s(getResult("Average Branch Length",i),9)+" "+d2s(getResult("Longest Shortest Path",i),9);\nFile.saveString(s,"{of}"); eval("script","System.exit(0);");\n')
        subprocess.run(f'JAVA_TOOL_OPTIONS= timeout 120 xvfb-run -a java -Dplugins.dir=/tmp/ijplug -cp /tmp/ij.jar ij.ImageJ -macro {mf}', shell=True, capture_output=True)
        rr = open(of).read().splitlines() if os.path.exists(of) else []
        ijr = [list(map(float, l.split()[1:])) for l in rr if l.startswith('R ')]
        chk(len(ijr) == len(T) and all(abs(a[0] - b[0]) < 1e-9 and abs(a[1] - b[5]) < 1e-6 and abs(a[2] - b[9]) < 1e-6 for a, b in zip(ijr, T)), '.ijm (Skeletonize → Analyze Skeleton) 을 Fiji 플러그인 ImageJ 로 실행 → 같은 표', rr[:3])
        # 가지치기: (a) shortest branch + Prune ends (b) lowest intensity branch, 원본 = g8.tif → SkelDriver 와 같은 표, (b) 는 .ijm 을 Fiji 로 실행해도 같은 표
        gid = pg.evaluate("IMG_UI._state().images.find(i=>i.name==='g8.tif').id")
        for nm, prm, drv in [('shortest branch + Prune ends', {'prune': 'shortest branch', 'pruneEnds': True}, '1 1'), ('lowest intensity branch (원본 g8)', {'prune': 'lowest intensity branch', 'orig': gid, 'pruneEnds': False}, '3 0')]:
            sel('bin8.tif'); n0 = pg.evaluate('IMG_UI._state().skelTrees.length')
            P = dict(skel=True, lsp=True, tagged=False, lspImg=False, **prm)
            pg.evaluate("P=>IMG_UI._run('skeleton',P)", P); pg.wait_for_function(f'IMG_UI._state().skelTrees.length>{n0}', timeout=30000); pg.wait_for_timeout(200)
            tr = pg.evaluate(f"IMG_UI._state().skelTrees.slice({n0})")
            lg = subprocess.run(f'JAVA_TOOL_OPTIONS= java -cp /tmp/ij.jar:/tmp/askel_out SkelDriver {skf} {drv} 1 1 1 1 {os.path.join(REF, "in", "g8.tif")}', shell=True, capture_output=True, text=True).stdout.splitlines()
            T2 = [list(map(float, l.split()[2:])) for l in lg if l.startswith('T ')]
            a2 = [[r['# Branches'], r['# Junctions'], r['# End-point voxels'], r['# Junction voxels'], r['# Slab voxels'], r['Average Branch Length'], r['# Triple points'], r['# Quadruple points'], r['Maximum Branch Length'], r['Longest Shortest Path'], r['spx'], r['spy'], r['spz']] for r in tr]
            okp = len(a2) == len(T2) and all(np.allclose(x, y, rtol=1e-12, atol=1e-12) for x, y in zip(a2, T2))
            chk(okp and len(T2) > 3, f'Analyze Skeleton 가지치기 {nm}: 표 = Fiji ({len(T2)} 골격, 가지 {int(sum(t[0] for t in T2))})', [len(a2), len(T2)])
            if 'orig' in prm:
                ijm2 = pg.evaluate("P=>IMG_UI._ops.skeleton.ij(P)", dict(P, origName='g8.tif')); of2 = os.path.join(tdk, 'out2.txt'); mf2 = os.path.join(tdk, 'm2.ijm')
                open(mf2, 'w').write(f'open("{os.path.join(REF, "in", "g8.tif")}"); open("{os.path.join(REF, "in", "bin8.tif")}"); setOption("BlackBackground", true);\n' + ijm2 + f'\ns="N "+nResults; for(i=0;i<nResults;i++) s=s+"\\nR "+getResult("# Branches",i)+" "+d2s(getResult("Average Branch Length",i),9); s=s+"\\nT "+getTitle()+" "+getValue("Mean");\nFile.saveString(s,"{of2}"); eval("script","System.exit(0);");\n')
                subprocess.run(f'JAVA_TOOL_OPTIONS= timeout 120 xvfb-run -a java -Dplugins.dir=/tmp/ijplug -cp /tmp/ij.jar ij.ImageJ -macro {mf2}', shell=True, capture_output=True)
                rr2 = open(of2).read().splitlines() if os.path.exists(of2) else []
                ij2 = [list(map(float, l.split()[1:])) for l in rr2 if l.startswith('R ')]
                tl = [l for l in rr2 if l.startswith('T ')]
                chk(len(ij2) == len(T2) and all(abs(x[0] - y[0]) < 1e-9 and abs(x[1] - y[5]) < 1e-6 for x, y in zip(ij2, T2)) and tl and tl[0].split()[1] == 'bin8.tif',
                    '.ijm 가지치기(original_image) 를 Fiji 로 실행 → 같은 표 · 끝나면 원본 영상 다시 선택', [len(ij2), len(T2), tl])
        # z step 없는 3D 스택: 1 화소로 세고 경고 (ok 알림 뒤 = 맨 위), LSP 영상만 골라도 LSP 계산, Slice 열
        sk2 = pg.evaluate("""(async()=>{ const w=12,h=10,d=3, pl=[]; for(let z=0;z<d;z++){ const a=new Float32Array(w*h); if(z===0) for(let x=1;x<11;x++) a[5*w+x]=255; else a[5*w+10]=255; pl.push(a); } /* x 1..10 (z0) + (10,5) 에서 z 로 2 화소 — 26-이웃이라 (9,5,0)·(10,5,0)·(10,5,1) 이 접점: Fiji SkelDriver LSP = 10.414213562373096 */
            IMG_UI._addImage({name:'skz.tif',w,h,bit:8,nc:1,nz:d,nt:1,planes:pl,cal:null}); const n0=IMG_UI._state().skelTrees.length;
            await IMG_UI._run('skeleton',{skel:false,lsp:false,tagged:false,lspImg:true}); const S=IMG_UI._state(); const t=S.skelTrees.slice(n0);
            return {n:t.length, lsp:t.map(r=>r['Longest Shortest Path']), msgs:IMG_UI._msgs().slice(0,3).map(m=>m.k+':'+m.m), lspimg:!!S.images.find(i=>i.name==='Longest shortest paths of skz.tif')}; })()""")
        chk(sk2['n'] == 1 and abs(sk2['lsp'][0] - 10.414213562373096) < 1e-12 and sk2['lspimg'] and sk2['msgs'][0].startswith('warn:z 간격을 모릅니다') and sk2['msgs'][1].startswith('ok:'),
            'Analyze Skeleton 3D (z step 없음): 1 화소로 세고 경고는 ok 뒤 · LSP 영상만 골라도 LSP 계산', sk2)
    # Wand (Worker): 이진 영상의 입자 몇 개를 클릭 → 외곽선 = ImageJ doWand(8-connected, tolerance 0) 의 꼭짓점과 같음
    sel('bin8.tif'); wb = pg.evaluate("""(async()=>{ const im=IMG_UI._cur(); const d=IMG_UI._px.get(im.id)[0]; const pts=[]; for(let i=0;i<d.length&&pts.length<4;i+=997) if(d[i]>0) pts.push([i%im.w,Math.floor(i/im.w)]);
        const out=[]; for(const [x,y] of pts){ await IMG_UI._wand(im,x,y); const S=IMG_UI._state(); const a=S.roiSel[im.id]; out.push({x,y,pts:a&&a.pts}); } return out; })()""")
    if not wb or any(w['pts'] is None for w in wb):
        wb2 = pg.evaluate("(()=>{ const S=IMG_UI._state(); return Object.keys(S).filter(k=>/act|roi/i.test(k)); })()"); chk(False, 'Wand: ROI 를 찾지 못함', wb2)
    else:
        tdw = tempfile.mkdtemp(); mf = os.path.join(tdw, 'w.ijm'); of = os.path.join(tdw, 'w.txt')
        open(mf, 'w').write(f'open("{os.path.join(REF, "in", "bin8.tif")}"); s="";\n' + ''.join(f'doWand({w["x"]},{w["y"]},0,"8-connected"); getSelectionCoordinates(xs,ys); s=s+"W"; for(i=0;i<xs.length;i++) s=s+" "+xs[i]+","+ys[i]; s=s+"\\n";\n' for w in wb) + f'File.saveString(s,"{of}");\n')
        subprocess.run(f'JAVA_TOOL_OPTIONS= xvfb-run -a java -Xmx1g -cp /tmp/ij.jar ij.ImageJ -batch {mf}', shell=True, capture_output=True)
        ijw = [sorted(tuple(map(int, t.split(','))) for t in l.split()[1:]) for l in open(of).read().splitlines() if l.startswith('W')] if os.path.exists(of) else []
        appw = [sorted(tuple(int(round(v)) for v in p) for p in w['pts']) for w in wb]
        chk(len(ijw) == len(appw) == len(wb) and all(a == b for a, b in zip(appw, ijw)), f'Wand (Worker) {len(wb)} 입자: 외곽선 꼭짓점 = ImageJ doWand', [len(a) for a in appw])
    # Sample: 사각형 ROI → 범위, Auto (HSB = Brightness 만)
    sel('trgb.tif'); pg.evaluate("IMG_UI._state().ui.panel='threshold'; IMG_UI._setActive({type:'rect',pts:[[5,5],[25,5],[25,20],[5,20]]}); IMG_UI.rerender()"); pg.wait_for_selector('#ctSample', timeout=10000); pg.wait_for_selector('#ctApply:not([disabled])', timeout=20000)
    pg.click('#ctSample'); pg.wait_for_timeout(300)
    smp = pg.evaluate("""(()=>{ const S=IMG_UI._state(); const src=S.images.find(i=>i.name==='trgb.tif'); const P=IMG_UI._px.get(src.id); const IC=self.IMG_CORE; const ch=IC.colorChannels(P[0],P[1],P[2],'HSB');
        const idx=[]; for(let y=5;y<20;y++) for(let x=5;x<25;x++) idx.push(y*src.w+x); const s=IC.colorSample(ch,idx,'HSB'); const C=S.ui.cthr; return {ok:JSON.stringify([s.min,s.max,s.pass0])===JSON.stringify([C.min,C.max,C.pass[0]]),C,s}; })()""")
    chk(smp['ok'], 'Sample (사각형 ROI) → 채널 범위 = ImageJ sample 규칙', smp)
    pg.click('#ctReset'); pg.select_option('#ctM', 'Otsu'); pg.click('#ctAuto'); pg.wait_for_timeout(200)
    au = pg.evaluate('JSON.parse(JSON.stringify(IMG_UI._state().ui.cthr))'); chk(au['min'][:2] == [0, 0] and au['max'][:2] == [255, 255] and au['max'][2] < 255, 'Auto (HSB): Brightness 만 바꿈', au)
    pg.evaluate("IMG_UI._setActive(null)")
    sel('n_t3z4c2.nd2'); pg.evaluate("IMG_UI._run('threshold',{method:'Default',dark:true,stack:true,newImage:false})"); pg.wait_for_timeout(600)
    chk(pg.evaluate('(()=>{ const im=IMG_UI._cur(); return im.nc===1 && !im.chNames; })()'), '채널 수가 바뀌면 채널 이름 지움', pg.evaluate('(()=>{ const im=IMG_UI._cur(); return [im.nc, im.chNames]; })()'))
    # 7d) Register (StackReg rigid) — 정답을 아는 drift stack, 결과 화소 = pystackreg
    pg.set_input_files('#imFin', [os.path.join(IOD, 'reg_rigid.tif')]); pg.wait_for_function("IMG_UI._state().images.some(i=>i.name==='reg_rigid.tif')", timeout=20000)
    sel('reg_rigid.tif'); pg.evaluate("IMG_UI._run('register',{method:'rigid',ref:'previous',nFirst:1,channel:1})")
    pg.wait_for_function("IMG_UI._state().lastReg && IMG_UI._state().lastReg.rows.length===6", timeout=120000); pg.wait_for_timeout(300)
    RG = json.load(open(os.path.join(IOD, 'reg_rigid.json'))); rows = pg.evaluate("IMG_UI._state().lastReg.rows")
    ea = max(abs(r_['angle'] - a_) for r_, a_ in zip(rows, RG['angle'])); ed = max(np.hypot(r_['dx'] - x_, r_['dy'] - y_) for r_, x_, y_ in zip(rows, RG['dx'], RG['dy']))
    chk(ea < 0.05 and ed < 0.05, 'Register rigid: 정답 회전(반시계 +)·중심 drift 회복', f'각도 최대 오차 {ea:.4f}°, drift 최대 오차 {ed:.4f} px')
    a = np.array(pg.evaluate("(()=>{ const im=IMG_UI._cur(); return [].concat(...IMG_UI._px.get(im.id).map(p=>Array.from(p))); })()"), np.float32)
    e = np.fromfile(os.path.join(IOD, 'reg_rigid.expected.f32'), np.float32); rel = float(np.max(np.abs(a - e)) / (e.max() - e.min()))
    chk(a.shape == e.shape and rel < 1e-5, 'Register rigid 결과 화소 = pystackreg transform_stack', f'최대 상대차 {rel:.2e}')
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.bottom='plots'; IMG_UI.rerender(); })()"); pg.wait_for_timeout(400); pg.screenshot(path=os.path.join(out, '12_register.png'))
    # 7e) 전수 검수 회귀: 되돌리기 → 보정 복원, ±∞ 영상 표시, 16-bit 표시 = ImageJ, 아주 큰 영상 확대 = 원래 해상도, Esc 로 창 닫기
    pg.evaluate("(()=>{ const im=IMG_UI._addImage({name:'cal_undo',w:100,h:50,bit:8,planes:[new Float32Array(5000).fill(7)],cal:{pw:2,ph:2,unit:'nm',src:'t'}}); IMG_UI.rerender(); })()")
    pg.evaluate("IMG_UI._run('resize',{factor:0.5,interp:'bilinear'})"); pg.wait_for_timeout(500); pg.evaluate('document.querySelector("#imUndo").click()'); pg.wait_for_timeout(300)
    cu = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return [im.w,im.h,im.cal&&im.cal.pw]; })()'); chk(cu == [100, 50, 2], 'Scale 뒤 되돌리기 → 크기·보정 모두 복원', cu)
    pg.evaluate("(()=>{ const d=new Float32Array(64*64); for(let i=0;i<d.length;i++) d[i]=i%64; d[5]=Infinity; d[6]=-Infinity; d[7]=NaN; IMG_UI._addImage({name:'inf32',w:64,h:64,bit:32,planes:[d]}); IMG_UI.rerender(); })()"); pg.wait_for_timeout(400)
    dr = pg.evaluate('(()=>{ const im=IMG_UI._cur(); return [im.disp[0].min, im.disp[0].max]; })()'); chk(dr == [0, 63], '±∞ 가 있는 32-bit: 표시 범위는 유한값으로 (ImageJ 와 같음)', dr)
    # 16-bit 표시: (int)((v−min)·256/(max−min+1)+0.5)
    pg.evaluate("(()=>{ const d=new Float32Array(256*4); for(let i=0;i<d.length;i++) d[i]=1000+i*3; const im=IMG_UI._addImage({name:'lut16',w:256,h:4,bit:16,planes:[d]}); im.disp[0].min=1100; im.disp[0].max=3000; IMG_UI._state().ui.views[im.id]={zoom:1,x:0,y:0}; IMG_UI.rerender(); })()"); pg.wait_for_timeout(400)
    px = pg.evaluate("(()=>{ const cv=document.querySelector('#imCanvas'); const g=cv.getContext('2d'); const dpr=window.devicePixelRatio||1; const out=[]; for(let x=0;x<256;x+=5){ out.push(g.getImageData(Math.floor((x+0.5)*dpr),Math.floor(1.5*dpr),1,1).data[0]); } return out; })()")
    exp = [min(255, int(max(0, 1000 + (256 + x) * 3 - 1100) * 256 / (3000 - 1100 + 1) + 0.5)) for x in range(0, 256, 5)]
    chk(px == exp, '16-bit 표시 = ImageJ ShortProcessor 규칙', f'{sum(a != b for a, b in zip(px, exp))} 곳 다름 {px[:8]} {exp[:8]}')
    # 아주 큰 영상(> 2.4×10⁷ 화소): 확대하면 원래 해상도
    pg.evaluate("(()=>{ const W=5200,H=5000; const d=new Float32Array(W*H); for(let y=0;y<H;y++) for(let x=0;x<W;x++) d[y*W+x]=((x+y)&1)?255:0; const im=IMG_UI._addImage({name:'huge',w:W,h:H,bit:8,planes:[d]}); IMG_UI._state().ui.views[im.id]={zoom:8,x:1001,y:1001}; IMG_UI.rerender(); })()"); pg.wait_for_timeout(1500)
    hv = pg.evaluate("(()=>{ const cv=document.querySelector('#imCanvas'); const g=cv.getContext('2d'); const dpr=window.devicePixelRatio||1; const v=[]; for(let k=0;k<6;k++) v.push(g.getImageData(Math.floor((k*8+4)*dpr),Math.floor(4*dpr),1,1).data[0]); return v; })()")
    chk(hv == [255, 0, 255, 0, 255, 0] or hv == [0, 255, 0, 255, 0, 255], '2.6×10⁷ 화소 영상 확대 → 원래 해상도 (바둑판 그대로)', hv)
    pg.evaluate("(()=>{ const S=IMG_UI._state(); const im=IMG_UI._cur(); IMG_UI._px.delete(im.id); S.images.pop(); S.sel=S.images[S.images.length-1].id; IMG_UI.rerender(); })()")
    pg.click('#imFmts'); pg.wait_for_timeout(200); pg.keyboard.press('Escape'); pg.wait_for_timeout(200)
    chk(pg.evaluate('document.querySelector("#imModal").hidden'), 'Esc 로 창 닫기')
    # 8) 새로 고침 → 복원
    pg.evaluate('IMG_UI._flush()'); pg.wait_for_timeout(800)
    sel('g8.tif'); pg.evaluate("IMG_UI._run('enhance',{saturated:0.35,normalize:false})"); pg.wait_for_timeout(500)   # 표시만 바꾸는 연산 뒤 저장 (예전: ver 만 올라 새로 고침 때 영상이 버려짐)
    n_before = pg.evaluate('IMG_UI._state().images.length'); pg.evaluate('IMG_UI._flush()'); pg.wait_for_timeout(1500)
    pg.reload(); pg.wait_for_timeout(800); pg.click('#nav button[data-v="imgx"]')
    pg.wait_for_function('window.IMG_UI && IMG_UI._state().images.length>0', timeout=30000); pg.wait_for_timeout(500)
    n_after = pg.evaluate('IMG_UI._state().images.length'); chk(n_after == n_before, '새로 고침 후 영상 복원', f'{n_after}/{n_before}')
    chk(pg.evaluate('IMG_UI._px.size') == n_after, '화소 자료 복원')
    # 9) 모바일
    pg.set_viewport_size({'width': 400, 'height': 860}); pg.wait_for_timeout(500); pg.screenshot(path=os.path.join(out, '10_mobile.png'), full_page=True)
    ovf = pg.evaluate('document.documentElement.scrollWidth - document.documentElement.clientWidth'); chk(ovf <= 2, '모바일 가로 넘침 없음', ovf)
    b.close()
res['errors'] = errs; chk(not errs, '콘솔 오류 0', errs[:5])
json.dump(res, open(os.path.join(out, 'ui_result.json'), 'w'), ensure_ascii=False, indent=1)
print('\n'.join(('ok  ' if c[0] else 'FAIL') + ' ' + c[1] + ('  (' + c[2] + ')' if c[2] and not c[0] else '') for c in res['checks']))
print('FAILS', len(fails))
