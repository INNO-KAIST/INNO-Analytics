"""Analyze Particles 보호 장치 · ROI 안만 세기 · 문턱 극성 안내 (Playwright · Chromium) — 10-19
사용: python3 test/ui_pa_guard.py <주입된 html> <결과 폴더>
합성 TEM 명시야 영상(밝은 바탕, 어두운 원판 400 개 — 붙은 쌍 포함, 잡음 σ 25, 왼쪽 아래 축척 막대·'100 nm' 비슷한 글자)을 만들어
① 기본 문턱(Dark background 꺼짐)·흐림 없이 → 잡음 점 수만 개 → 상한에서 몇 초 안에 멈추고 안내
② Dark background 켜면 → Threshold 패널에 '절반 넘게 선택' 경고, Analyze 도 상한에서 멈춤
③ Gaussian σ2 → Convert to Mask → Fill Holes → Watershed → 사각형 ROI(축척 막대 뺌) 안 Analyze → 개수 ≈ 정답
④ Max count 올리면 그대로 측정, 표는 앞 3000 줄만"""
import sys, os, time, json
import numpy as np
from PIL import Image
from playwright.sync_api import sync_playwright
html, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(7); W = H = 1024; R = 11.0
img = np.full((H, W), 175.0); cs = []
while len(cs) < 400:
    x, y = rng.uniform(30, W - 30), rng.uniform(30, H - 120)
    if all((x - a) ** 2 + (y - b) ** 2 > (2 * R + 3) ** 2 for a, b in cs): cs.append((x, y))
# 붙은 쌍 40 개 (거리 2R − 1)
pairs = []
for a, b in cs[:40]:
    t = rng.uniform(0, 2 * np.pi); x2, y2 = a + (2 * R - 1) * np.cos(t), b + (2 * R - 1) * np.sin(t)
    if 30 < x2 < W - 30 and 30 < y2 < H - 120 and all((x2 - c) ** 2 + (y2 - d) ** 2 > (2 * R - 2) ** 2 for c, d in cs if (c, d) != (a, b)): pairs.append((x2, y2))
allc = cs + pairs
yy, xx = np.mgrid[0:H, 0:W]
for a, b in allc: img[(xx - a) ** 2 + (yy - b) ** 2 <= R * R] = 70
img[H - 60:H - 54, 40:200] = 20  # 축척 막대
for k, x0 in enumerate([60, 95, 130]):  # '100' 글자 비슷한 고리·막대
    if k == 0: img[H - 100:H - 70, x0:x0 + 6] = 20
    else: img[(xx - x0 - 10) ** 2 / 100 + (yy - (H - 85)) ** 2 / 225 <= 1] = 20; img[(xx - x0 - 10) ** 2 / 36 + (yy - (H - 85)) ** 2 / 121 <= 1] = 175
img = np.clip(img + rng.normal(0, 25, img.shape), 0, 255).astype(np.uint8)
fn = os.path.join(out, 'pa_synth.tif'); Image.fromarray(img).save(fn)
truth = len(allc); print('정답 입자', truth, '(붙은 쌍', len(pairs), ')')
errs = []; ok = []
def chk(c, n, i=''):
    ok.append(bool(c)); print(('  ✓ ' if c else '  ✗ ') + n + (f' — {i}' if i != '' else ''))
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={'width': 1500, 'height': 1000})
    pg.on('pageerror', lambda e: errs.append(str(e)))
    pg.goto('file://' + os.path.abspath(html)); pg.wait_for_timeout(800)
    pg.click('#nav button[data-v="imgx"]'); pg.wait_for_selector('.imx')
    pg.set_input_files('#imFin', [fn]); pg.wait_for_function('IMG_UI._state().images.length>=1', timeout=30000); pg.wait_for_timeout(500)
    idle = lambda: pg.wait_for_function('!IMG_UI._busy()', timeout=60000)
    chk(pg.evaluate('IMG_UI._state().ui.thr.dark') is False and pg.evaluate('[IMG_UI._state().ui.pa.minCirc,IMG_UI._state().ui.pa.exclude]')==[0.5,True], '새 상태 기본값 = Dark background 꺼짐 · Circularity 0.5–1 · Exclude on edges')
    def analyze():
        pg.evaluate("(()=>{const S=IMG_UI._state(); S.ui.panel='analyze'; IMG_UI.rerender();})()"); pg.wait_for_selector('#paRun')
        t = time.time(); pg.click('#paRun'); idle(); pg.wait_for_timeout(200); return time.time() - t
    dt = analyze(); m = pg.evaluate('IMG_UI._msgs()[0].m'); pa = pg.evaluate('!!IMG_UI._state().particles[IMG_UI._cur().id]')
    chk(dt < 15 and '상한' in m and 'Gaussian' in m and not pa, f'흐림 없이 → 잡음 점 상한에서 멈춤 ({dt:.1f} s), 안내', m[:110])
    pg.evaluate("(()=>{const S=IMG_UI._state(); S.ui.panel='threshold'; S.ui.thr.dark=true; IMG_UI.rerender();})()"); pg.wait_for_timeout(900)
    chk('절반이 넘게 선택' in pg.inner_text('#imPanel') and 'Dark background' in pg.inner_text('#imPanel .note.warn'), 'Dark background 켜면 Threshold 패널 경고')
    pg.screenshot(path=os.path.join(out, 'pa_thr_warn.png'))
    dt = analyze(); m = pg.evaluate('IMG_UI._msgs()[0].m'); chk('바탕' in m and 'Dark background' in m, f'Analyze 도 극성 안내 ({dt:.1f} s)', m[:90])
    pg.evaluate("(()=>{const S=IMG_UI._state(); S.ui.thr.dark=false; IMG_UI.rerender();})()")
    run = lambda op, P: (pg.evaluate('([o,P])=>IMG_UI._run(o,P)', [op, P]), idle(), pg.wait_for_timeout(250))
    run('gaussian', {'sigma': 2, 'scaled': False, 'stack': False}); run('threshold', {'method': 'Default', 'dark': False, 'manual': False})
    run('binary', {'fn': 'fill', 'iterations': 1, 'count': 1, 'edm32': False, 'stack': False}); run('binary', {'fn': 'watershed', 'iterations': 1, 'count': 1, 'edm32': False, 'stack': False})
    pg.evaluate("(()=>{const A=IMG_UI._state().ui.pa; Object.assign(A,{minSize:100,maxSize:'',unitPx:true,minCirc:0.5,maxCirc:1,exclude:false});})()")
    analyze(); n_all = pg.evaluate('IMG_UI._state().particles[IMG_UI._cur().id].rows.length')
    pg.evaluate("r=>IMG_UI._setActive(r)", {'id': 'R', 'type': 'rect', 'pts': [[0, 0], [W, H - 115]]}); analyze()
    pa = pg.evaluate('(()=>{const p=IMG_UI._state().particles[IMG_UI._cur().id]; return {n:p.rows.length,thr:p.summary.Threshold}})()')
    chk(n_all > pa['n'], 'ROI 없이는 축척 막대·글자도 입자로 셈 → 사각형 ROI 안만 세면 빠짐', f"{n_all} → {pa['n']}")
    chk(abs(pa['n'] - truth) <= max(4, 0.02 * truth) and 'ROI' in pa['thr'], f'Gaussian → Mask → Fill → Watershed → ROI 안 = 정답 {truth} 개 근처', pa)
    rec = pg.evaluate('IMG_UI._steps().filter(s=>s.op==="analyze").slice(-1)[0].params.roi'); chk(rec and rec['type'] == 'rect', '매크로 기록에 ROI 포함', rec)
    pg.screenshot(path=os.path.join(out, 'pa_ok.png'))
    # Max count 올리면 그대로 (Undo 4 번 → 원본, 흐림 없이, 상한 200000) — 표는 앞 3000 줄
    for _ in range(4): pg.evaluate('IMG_UI._undo()'); pg.wait_for_timeout(150)
    pg.evaluate("(()=>{ const S=IMG_UI._state(); IMG_UI._setActive(null); Object.assign(S.ui.pa,{minSize:0,minCirc:0,maxCount:200000}); })()")
    dt = analyze(); n = pg.evaluate('(IMG_UI._state().particles[IMG_UI._cur().id]||{rows:[]}).rows.length')
    pg.evaluate("(()=>{const S=IMG_UI._state(); S.ui.bottom='particles'; S.ui.bottomOpen=true; IMG_UI.rerender();})()"); pg.wait_for_timeout(1500)
    tr = pg.locator('#imBottom tbody tr').count(); txt = pg.inner_text('#imBottom')
    chk(n > 3000 and tr == 3000 and '앞 3,000 줄만' in txt, f'Max count 올리면 측정({n:,} 개, {dt:.1f} s), 표는 앞 3000 줄', tr)
    # ⑤ Auto (recommended) 단추: 원본에서 한 번에 — 원본은 그대로, ROI 이어받음
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.sel=S.images[0].id; IMG_UI.rerender(); })()"); pg.wait_for_timeout(300)
    v0 = pg.evaluate('IMG_UI._cur().ver'); nimg = pg.evaluate('IMG_UI._state().images.length')
    pg.evaluate("r=>IMG_UI._setActive(r)", {'id': 'R2', 'type': 'rect', 'pts': [[0, 0], [W, H - 115]]})
    pg.evaluate("(()=>{const S=IMG_UI._state(); S.ui.panel='analyze'; IMG_UI.rerender();})()"); pg.wait_for_selector('#paAuto')
    t = time.time(); pg.click('#paAuto')
    pg.wait_for_function('IMG_UI._state().images.length>%d && !IMG_UI._busy() && IMG_UI._state().particles[IMG_UI._cur().id]' % nimg, timeout=120000); pg.wait_for_timeout(500)
    r = pg.evaluate('(()=>{const im=IMG_UI._cur(); const p=IMG_UI._state().particles[im.id]; const o=IMG_UI._state().images[0]; return {n:p.rows.length,name:im.name,thr:p.summary.Threshold,orig:o.ver,mean:p.rows[0].Mean,min:IMG_UI._state().ui.pa.minSize}})()')
    chk(abs(r['n'] - truth) <= max(4, 0.02 * truth) and r['orig'] == v0 and 'ROI' in r['thr'], f'Auto (recommended): {time.time()-t:.1f} s, 정답 {truth} 근처 · 원본 그대로 · ROI 이어받음', r)
    chk(r['mean'] is not None and r['mean'] < 150, '세기(Mean)는 원본에서 (입자는 어두움 ≈ 70)', r['mean'])
    pg.screenshot(path=os.path.join(out, 'pa_auto.png'))
    # ⑥ ROI 없이 Auto: 축척 막대·글자를 스스로 빼고, 크기 분포 표를 연다
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.sel=S.images[0].id; IMG_UI._setActive(null); IMG_UI.rerender(); })()"); pg.wait_for_timeout(300)
    nimg = pg.evaluate('IMG_UI._state().images.length'); pg.evaluate("(()=>{const S=IMG_UI._state(); S.ui.panel='analyze'; IMG_UI.rerender();})()"); pg.click('#paAuto')
    pg.wait_for_function('IMG_UI._state().images.length>%d && !IMG_UI._busy() && IMG_UI._state().particles[IMG_UI._cur().id]' % nimg, timeout=120000); pg.wait_for_timeout(500)
    r = pg.evaluate('(()=>{const im=IMG_UI._cur(); const p=IMG_UI._state().particles[im.id]; return {n:p.rows.length,bottom:IMG_UI._state().ui.bottom,msg:IMG_UI._msgs()[0].m}})()')
    chk(abs(r['n'] - truth) <= max(4, 0.02 * truth) and '축척 막대' in r['msg'] and r['bottom'] == 'dist', f'ROI 없이 Auto: 축척 막대 스스로 뺌 · 정답 {truth} 근처 · Size distribution 열림', r)
    # ⑦ 겹친 입자·크기 이상치는 찾되 분포에서 뺌 · 이미지에서 눌러 고른 입자만 분포
    pg.evaluate("(()=>{const S=IMG_UI._state(); S.ui.panel='analyze'; IMG_UI.rerender();})()"); pg.wait_for_selector('#pcPick')
    r = pg.evaluate('(()=>{const im=IMG_UI._cur(); const p=IMG_UI._state().particles[im.id]; return {n:p.rows.length,ovl:p.ovl?p.ovl.reduce((a,b)=>a+b,0):-1,card:document.querySelector("#imPanel .im-card").innerText}})()')
    chk(r['ovl'] > 0 and 'overlapping' in r['card'] and 'used' in r['card'], f"겹친 입자 표시 {r['ovl']} 개 · 패널에 used/overlapping", r['card'][:120])
    opts = pg.evaluate('[...document.querySelectorAll("#dsUse option")].map(o=>o.textContent)')
    nused = pg.evaluate('(()=>{const p=IMG_UI._state().particles[IMG_UI._cur().id]; return p.rows.length})()')
    chk(len(opts) == 3 and opts[0].startswith('Used') and f'All found ({nused})' in opts[1], '분포: Used · All found · Picked 선택', opts)
    nU = int(opts[0].split('(')[1].rstrip(')')); chk(nU <= nused - r['ovl'], f'Used = 찾은 {nused} − 겹침 {r["ovl"]} − 이상치 → {nU}')
    pg.click('#pcOvl'); pg.wait_for_timeout(300); o2 = pg.evaluate('document.querySelector("#dsUse option").textContent'); chk(int(o2.split('(')[1].rstrip(')')) > nU, 'Exclude overlapping 끄면 Used 늘어남', o2)
    pg.click('#pcOvl'); pg.wait_for_timeout(300)
    pg.fill('#pcK', '0'); pg.dispatch_event('#pcK', 'change'); pg.wait_for_timeout(300); o3 = pg.evaluate('document.querySelector("#dsUse option").textContent')
    chk(int(o3.split('(')[1].rstrip(')')) == nused - r['ovl'], 'Size outliers 0 = 끔 → Used = 찾은 − 겹침', o3)
    pg.fill('#pcK', '3'); pg.dispatch_event('#pcK', 'change'); pg.wait_for_timeout(300)
    # 고르기: 처음 세 입자 중심을 화면에서 누름
    pg.click('#pcPick'); pg.wait_for_timeout(200)
    cents = pg.evaluate('IMG_UI._state().particles[IMG_UI._cur().id].cent.slice(0,40).map(c=>IMG_UI._toScr(c[0],c[1]))')
    box = pg.locator('#imOverlay').bounding_box(); clicked = 0
    for sx, sy in cents:
        if 5 < sx < box['width'] - 5 and 5 < sy < box['height'] - 5 and clicked < 3:
            pg.mouse.click(box['x'] + sx, box['y'] + sy); pg.wait_for_timeout(150); clicked += 1
    pk = pg.evaluate('IMG_UI._state().particles[IMG_UI._cur().id].pick.length')
    chk(clicked == 3 and pk == 3, f'Pick: 이미지에서 3 개 누름 → {pk} 개 골라짐')
    pg.mouse.click(box['x'] + cents[0][0], box['y'] + cents[0][1]) if 5 < cents[0][0] < box['width'] - 5 and 5 < cents[0][1] < box['height'] - 5 else None
    pg.wait_for_timeout(150); pk2 = pg.evaluate('IMG_UI._state().particles[IMG_UI._cur().id].pick.length'); chk(pk2 == 2, '다시 누르면 빠짐', pk2)
    nroi = pg.evaluate('(IMG_UI._state().rois[IMG_UI._cur().id]||[]).length'); chk(nroi <= 1, '고르는 동안 ROI 를 그리지 않음', nroi)
    pg.click('#pcPick'); pg.wait_for_timeout(200); pg.click('#pcPickDist'); pg.wait_for_timeout(600)
    txt = pg.inner_text('#imBottom'); chk(pg.evaluate('document.querySelector("#dsUse").value') == 'picked' and '고른 2 개' in txt, 'Distribution of picked → 고른 2 개만', txt[:200])
    pg.evaluate("(()=>{const S=IMG_UI._state(); S.ui.bottom='particles'; IMG_UI.rerender();})()"); pg.wait_for_timeout(800)
    th = pg.evaluate('[...document.querySelectorAll("#imBottom thead th")].map(t=>t.textContent).slice(0,3)'); chk(th[:3] == ['#', 'Status', 'Picked'], 'Particles 표에 Status · Picked 열', th)
    pg.screenshot(path=os.path.join(out, 'pa_pick.png'))
    b.close()
chk(not errs, '페이지 오류 0', errs[:3])
print(f'\n{sum(ok)}/{len(ok)} 통과'); sys.exit(0 if all(ok) else 1)
