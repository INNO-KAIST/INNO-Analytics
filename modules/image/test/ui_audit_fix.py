"""전수조사(10-19) 고친 것 확인 (Playwright · Chromium)
사용: python3 test/ui_audit_fix.py <주입된 html> <결과 폴더>
① Auto 가 손으로 하는 Analyze 설정(S.ui.pa)을 바꾸지 않음 ② Convert to Mask 사본에서 Analyze → 세기는 원본(Mean ≠ 255)
③ 굵은 마스크로 Analyze Skeleton → 먼저 Skeletonize(가지 수 적음) ④ 밝은 바탕에 Subtract Background(Light 끔) → 경고·되돌려 다시
⑤ 32-bit 음수 이미지에 NaN Background 기본 0–255 → 멈춤 ⑥ Summary 열 단위"""
import sys, os
import numpy as np
from PIL import Image
from playwright.sync_api import sync_playwright
html, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(3); W = H = 512; R = 10.0
img = np.full((H, W), 175.0); cs = []
while len(cs) < 120:
    x, y = rng.uniform(20, W - 20), rng.uniform(20, H - 20)
    if all((x - a) ** 2 + (y - b) ** 2 > (2 * R + 4) ** 2 for a, b in cs): cs.append((x, y))
yy, xx = np.mgrid[0:H, 0:W]
for a, b in cs: img[(xx - a) ** 2 + (yy - b) ** 2 <= R * R] = 70
img = np.clip(img + rng.normal(0, 15, img.shape), 0, 255).astype(np.uint8)
fn = os.path.join(out, 'af_synth.tif'); Image.fromarray(img).save(fn)
errs = []; ok = []
def chk(c, n, i=''):
    ok.append(bool(c)); print(('  ✓ ' if c else '  ✗ ') + n + (f' — {i}' if i != '' else ''))
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={'width': 1500, 'height': 1000})
    pg.on('pageerror', lambda e: errs.append(str(e)))
    pg.goto('file://' + os.path.abspath(html)); pg.wait_for_timeout(800)
    pg.click('#nav button[data-v="imgx"]'); pg.wait_for_selector('.imx')
    pg.set_input_files('#imFin', [fn]); pg.wait_for_function('IMG_UI._state().images.length>=1', timeout=30000); pg.wait_for_timeout(400)
    idle = lambda: pg.wait_for_function('!IMG_UI._busy()', timeout=120000)
    run = lambda op, P: (pg.evaluate('([o,P])=>IMG_UI._run(o,P)', [op, P]), idle(), pg.wait_for_timeout(250))
    pa0 = pg.evaluate('JSON.stringify(IMG_UI._state().ui.pa)')
    pg.evaluate("(()=>{const S=IMG_UI._state(); S.ui.panel='analyze'; IMG_UI.rerender();})()"); pg.wait_for_selector('#paAuto')
    n = pg.evaluate('IMG_UI._state().images.length'); pg.click('#paAuto')
    pg.wait_for_function('IMG_UI._state().images.length>%d && !IMG_UI._busy() && IMG_UI._state().particles[IMG_UI._cur().id]' % n, timeout=120000); pg.wait_for_timeout(300)
    na = pg.evaluate('IMG_UI._state().particles[IMG_UI._cur().id].rows.length')
    chk(pg.evaluate('JSON.stringify(IMG_UI._state().ui.pa)') == pa0 and abs(na - 120) <= 3, f'① Auto 뒤에도 손 설정 그대로 · {na} 개', pa0[:80])
    # ② Convert to Mask 사본 → Analyze: 세기는 원본
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.sel=S.images[0].id; IMG_UI.rerender(); })()"); pg.wait_for_timeout(200)
    run('gaussian', {'sigma': 1, 'scaled': False, 'stack': False}); run('threshold', {'method': 'Default', 'dark': False, 'manual': False, 'newImage': True})
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.panel='analyze'; IMG_UI.rerender(); })()"); pg.wait_for_selector('#paRun'); pg.click('#paRun'); idle(); pg.wait_for_timeout(300)
    r = pg.evaluate('(()=>{ const p=IMG_UI._state().particles[IMG_UI._cur().id]; return p?{n:p.rows.length,mean:p.rows.reduce((a,x)=>a+x.Mean,0)/p.rows.length,red:p.redirect}:null })()')
    chk(r and r['mean'] < 150, '② Convert to Mask 사본 Analyze → Mean 은 원본에서(≈ 70, 255 아님)', r)
    # ③ 굵은 마스크 → Analyze Skeleton
    n = pg.evaluate('IMG_UI._state().images.length'); run('skeleton', {'skel': False, 'prune': 'none', 'pruneEnds': False, 'lsp': True, 'tagged': True, 'lspImg': False})
    m = pg.evaluate('IMG_UI._msgs()[0].m'); nb = pg.evaluate('IMG_UI._state().skelBranches.length')
    chk('먼저 Skeletonize' in m and nb < 2000, f'③ 굵은 마스크 → 먼저 Skeletonize, 가지 {nb}', m[:90])
    # ④ Subtract Background (Light 끔) on 밝은 바탕
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.sel=S.images[0].id; IMG_UI.rerender(); })()"); pg.wait_for_timeout(200)
    run('subtractBg', {'radius': 50, 'light': False, 'create': False, 'nosmooth': False, 'stack': False})
    m = pg.evaluate('IMG_UI._msgs()[0]'); chk(m['k'] == 'warn' and 'Light background' in m['m'] and m['act'], '④ 밝은 바탕 + Light 끔 → 경고 · 되돌려 다시 단추', m['m'][:70])
    pg.click('#imMsgA'); idle(); pg.wait_for_timeout(400)
    st = pg.evaluate('IMG_UI._steps().slice(-1)[0]'); chk(st['op'] == 'subtractBg' and st['params']['light'] is True, '④ 단추 → Light background 로 다시', st['params'])
    # ⑤ NaN Background
    run('type', {'to': '32'}); run('math', {'fn': 'sub', 'value': 200, 'stack': False})
    pg.evaluate("(()=>{ try{ IMG_UI._run('nanbg',{lower:0,upper:255}); }catch(e){} })()"); idle(); pg.wait_for_timeout(300)
    m = pg.evaluate('IMG_UI._msgs()[0].m'); chk('NaN' in m or '엉뚱한' in m, '⑤ 음수 32-bit + 0–255 → 멈춤', m[:80])
    # ⑥ Summary 단위
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.bottom='summary'; S.ui.bottomOpen=true; IMG_UI.rerender(); })()"); pg.wait_for_timeout(500)
    th = pg.evaluate('[...document.querySelectorAll("#imBottom thead th")].map(t=>t.textContent)'); chk(any('Density' in t and '/px²' in t for t in th), '⑥ Summary 밀도 열 단위', [t for t in th if 'Density' in t])
    b.close()
chk(not errs, '페이지 오류 0', errs[:3])
print(f'\n{sum(ok)}/{len(ok)} 통과'); sys.exit(0 if all(ok) else 1)
