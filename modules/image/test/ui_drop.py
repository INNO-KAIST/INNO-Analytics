"""Image Analysis ▸ Drop 패널 브라우저 실측 (Playwright · Chromium)
사용: python3 test/ui_drop.py <주입된 html> <결과 폴더>
합성 방울(ref/gen_drop.py 정답)을 화면에서 열어 접촉각·표면장력을 재고, 사각형 ROI · 선 ROI 기준선 · 스택 시계열 · 표 CSV ·
새로 고침 후 복원 · 휴대폰 폭 · 콘솔 오류 0 을 본다. 실측 OpenDrop 그림(test/drop_real, 저장소 제외)이 있으면 그것도."""
import sys, os, json
import numpy as np
from PIL import Image
from playwright.sync_api import sync_playwright
html, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
H = os.path.dirname(os.path.abspath(__file__)); DD = os.path.join(H, 'drop'); RD = os.path.join(H, 'drop_real')
T = {t['file']: t for t in json.load(open(os.path.join(DD, 'truth.json')))}
errs = []; fails = []; checks = []
def chk(c, name, info=''):
    checks.append([bool(c), name, str(info)]); print(('  ✓ ' if c else '  ✗ ') + name + (f' — {info}' if info != '' else ''))
    if not c: fails.append(name + ' — ' + str(info))
# 스택: pendant_03 을 0 · 1 · 2 px 옮긴 3 장 (같은 γ 이어야)
a = np.asarray(Image.open(os.path.join(DD, 'pendant_03.png')).convert('L'))
stk = os.path.join(out, 'pend_stack.tif'); fr = [Image.fromarray(np.roll(a, k, axis=1)) for k in range(3)]
fr[0].save(stk, save_all=True, append_images=fr[1:])
a2 = np.asarray(Image.open(os.path.join(DD, 'sessile_13.png')).convert('L')); sstk = os.path.join(out, 'sess_stack.tif'); f2 = [Image.fromarray(np.roll(a2, k, axis=1)) for k in range(3)]; f2[0].save(sstk, save_all=True, append_images=f2[1:])
files = [sstk, os.path.join(DD, 'sessile_13.png'), os.path.join(DD, 'sessile_15.png'), os.path.join(DD, 'pendant_03.png'), stk]
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={'width': 1500, 'height': 980})
    pg.on('console', lambda m: errs.append('console.' + m.type + ': ' + m.text) if m.type == 'error' and 'ERR_TUNNEL' not in m.text and 'ERR_NAME' not in m.text else None)
    pg.on('pageerror', lambda e: errs.append('pageerror: ' + str(e)))
    pg.goto('file://' + os.path.abspath(html)); pg.wait_for_timeout(800)
    pg.click('#nav button[data-v="imgx"]'); pg.wait_for_selector('.imx', timeout=10000)
    pg.set_input_files('#imFin', files)
    pg.wait_for_function('window.IMG_UI && IMG_UI._state().images.length>=5', timeout=30000)
    sel = lambda name: pg.evaluate("n=>{ const S=IMG_UI._state(); S.sel=S.images.find(i=>i.name===n).id; IMG_UI.rerender(); }", name)
    idle = lambda: pg.wait_for_function('!IMG_UI._busy()', timeout=60000)
    def run():
        n0 = pg.evaluate('IMG_UI._state().drops.length'); pg.click('#dpRun'); idle(); pg.wait_for_timeout(250); return pg.evaluate('IMG_UI._state().drops').__getitem__(slice(n0, None))
    # 1) 접촉각: 패널 탭 → 영상 전체 · 자동 기준선
    sel('sessile_13.png'); pg.click('#imTabs button[data-p="drop"]'); pg.wait_for_selector('#dpRun')
    chk(pg.locator('#dpMode button.active').inner_text() == 'Contact angle', '기본 모드 Contact angle')
    rows = run(); t = T['sessile_13.png']; by = {r['Method']: r for r in rows}
    chk(len(rows) == 4, '방법 4 개 = 표 4 줄', [r['Method'] for r in rows])
    yl = by.get('Young–Laplace', {}); chk(abs(yl.get('θ mean', 999) - t['theta']) < 1.0, f"YL θ ≈ 정답 {t['theta']}°", yl.get('θ mean'))
    chk(abs(by['Circle (spherical cap)']['θ mean'] - t['theta']) < 1.0, '원 θ (Bo ≈ 0) ≈ 정답', by['Circle (spherical cap)']['θ mean'])
    chk(all(r.get('Baseline', '').startswith('auto') for r in rows), '기준선 auto 기록')
    card = pg.inner_text('#imPanel .im-card'); chk('Young–Laplace' in card and '°' in card, '결과 카드(방법별 표)', card[:80].replace('\n', ' '))
    ov = pg.evaluate('(()=>{ const F=IMG_UI._state().dropFit[IMG_UI._cur().id]; return F&&{t:F.ov.tangents.length,c:F.ov.curves.length,b:!!F.ov.baseline}; })()')
    chk(ov and ov['t'] == 8 and ov['c'] == 2 and ov['b'], '덧그림: 접선 8 · 곡선 2 · 기준선', ov)
    pg.screenshot(path=os.path.join(out, '01_sessile.png'))
    # 덧그림 화소: 캔버스 위에 분홍(YL) 선이 그려졌는지
    pink = pg.evaluate('''(()=>{ const c=document.querySelector("#imOverlay"); const g=c.getContext("2d"); const d=g.getImageData(0,0,c.width,c.height).data; let n=0; for(let i=0;i<d.length;i+=4) if(d[i]>220&&d[i+1]<90&&d[i+2]>170&&d[i+3]>200) n++; return n; })()''')
    chk(pink > 200, '캔버스에 YL 곡선(분홍) 화소', pink)
    # 2) 방법 끄기 → 덧그림에서 사라짐, 다음 측정은 남은 방법만
    pg.click('[data-dm="circle"]'); pg.click('[data-dm="ellipse"]'); pg.wait_for_timeout(150)
    rows = run(); chk(sorted(r['Method'] for r in rows) == ['Tangent (polynomial)', 'Young–Laplace'], '고른 방법만', [r['Method'] for r in rows])
    pg.click('[data-dm="circle"]'); pg.click('[data-dm="ellipse"]')
    # 3) 사각형 ROI + 선 ROI 기준선 (Bo 0.6, θ 135 — 원·타원은 맞지 않고 YL 만 맞아야)
    sel('sessile_15.png'); t = T['sessile_15.png']; (p1, p2) = t['baseline']
    W, Hh = t['W'], t['H']
    pg.evaluate("r=>IMG_UI._setActive(r)", {'id': 'L1', 'type': 'line', 'pts': [[p1[0] - 60, p1[1] - (p2[1] - p1[1]) / (p2[0] - p1[0]) * 60], [p2[0] + 60, p2[1] + (p2[1] - p1[1]) / (p2[0] - p1[0]) * 60]]})
    pg.wait_for_timeout(150); pg.click('#dpBl'); pg.wait_for_timeout(150)
    chk('line' in pg.inner_text('#imPanel'), '기준선 = 선 ROI 로 정함')
    pg.evaluate("r=>IMG_UI._setActive(r)", {'id': 'R1', 'type': 'rect', 'pts': [[30, 20], [W - 30, Hh - 40]]}); pg.wait_for_timeout(150)
    chk('rect' in pg.inner_text('#imPanel'), '영역 = 사각형 ROI')
    rows = run(); by = {r['Method']: r for r in rows}
    chk(all(r.get('Baseline', '').startswith('line') for r in rows), '기준선 line 기록')
    chk(abs(by['Young–Laplace']['θ mean'] - t['theta']) < 1.0, f"YL θ (Bo {t['Bo']}) ≈ 정답 {t['theta']}°", by['Young–Laplace']['θ mean'])
    chk(by['Circle (spherical cap)']['θ mean'] < t['theta'] - 5, '원은 중력 납작 방울에서 작게 나옴(정보)', by['Circle (spherical cap)']['θ mean'])
    warn = pg.evaluate('IMG_UI._state().dropFit[IMG_UI._cur().id].warn'); chk(any('방법 사이' in w for w in warn), '방법 차이 경고', warn)
    pg.screenshot(path=os.path.join(out, '02_sessile_roi_line.png'))
    pg.click('#dpBlA'); pg.wait_for_timeout(100); chk('auto' in pg.inner_text('#imPanel'), 'Auto 로 되돌림')
    # 4) 펜던트: 바늘 지름(정답)으로 γ
    sel('pendant_03.png'); t = T['pendant_03.png']; pg.click('#dpMode button[data-m="pendant"]'); pg.wait_for_selector('#dpNd')
    pg.fill('#dpNd', str(t['needle_mm'])); pg.dispatch_event('#dpNd', 'change')
    pg.click('#dpPhys summary'); pg.fill('#dpRi', '1000'); pg.dispatch_event('#dpRi', 'change'); pg.fill('#dpRo', '0'); pg.dispatch_event('#dpRo', 'change')
    rows = run(); r = rows[0]; dg = (r['γ'] / t['gamma_mNm'] - 1) * 100
    chk(abs(dg) < 0.5, f"γ ≈ 정답 {t['gamma_mNm']:.2f} mN/m", f"{r['γ']:.3f} ({dg:+.2f} %)")
    chk(abs(r['V (mm³)'] / t['V_mm3'] - 1) < 0.01, '부피 ≈ 정답 (1 %)', f"{r['V (mm³)']:.4f} / {t['V_mm3']:.4f}")
    chk(abs(r['Tilt (°)'] - t['tilt_deg']) < 0.1, '기울기 ≈ 정답', f"{r['Tilt (°)']:.3f} / {t['tilt_deg']:.3f}")
    chk(r['Scale'] == 'needle', '눈금 = 바늘')
    st = pg.inner_text('#imPanel .im-card'); chk('mN/m' in st and 'Worthington' in st, '펜던트 결과 카드', st[:60].replace('\n', ' '))
    pg.screenshot(path=os.path.join(out, '03_pendant.png'))
    # 5) 스택 시계열 (3 장, 같은 γ)
    sel('pend_stack.tif'); pg.wait_for_selector('#dpFr'); pg.select_option('#dpFr', 'all')
    pg.fill('#dpNd', str(t['needle_mm'])); pg.dispatch_event('#dpNd', 'change')
    rows = run(); g = [q['γ'] for q in rows if q.get('γ') is not None]
    chk(len(rows) == 3 and len(g) == 3 and max(g) - min(g) < 0.05 * np.mean(g) / 10, '스택 3 장 · γ 같음', [round(x, 3) for x in g])
    plots = pg.evaluate('IMG_UI._state().plots.map(p=>p.title)'); chk(any('Surface tension' in x for x in plots), '시계열 그래프 Plots 에', plots[-1:] if plots else plots)
    # 5b) 세실 스택 시계열 → 그래프 (주 방법 이름으로)
    sel('sess_stack.tif'); pg.click('#dpMode button[data-m="sessile"]'); pg.wait_for_selector('#dpFr'); pg.select_option('#dpFr', 'all')
    rows = run(); yl3 = [q['θ mean'] for q in rows if q.get('Method') == 'Young–Laplace']
    chk(len(yl3) == 3 and max(yl3) - min(yl3) < 0.3, '세실 스택 3 장 YL θ 같음', [round(x, 2) for x in yl3])
    plots = pg.evaluate('IMG_UI._state().plots.map(p=>p.title)'); chk(any('Contact angle (Young–Laplace)' in x for x in plots), '세실 시계열 그래프', plots[-1:])
    card = pg.inner_text('#imPanel .im-card'); chk(card.count('Young–Laplace') == 1, '카드는 이번 측정 한 장 것만', card[:80].replace('\n', ' '))
    # 화소가 바뀌면(필터) 옛 덧그림 안 그림
    pg.select_option('#dpFr', 'current'); sel('sessile_13.png'); pg.click('#dpRun'); idle(); pg.wait_for_timeout(200)
    pg.evaluate("IMG_UI._run('gaussian',{sigma:1})"); idle(); pg.wait_for_timeout(300)
    pink2 = pg.evaluate('''(()=>{ const c=document.querySelector("#imOverlay"); const g=c.getContext("2d"); const d=g.getImageData(0,0,c.width,c.height).data; let n=0; for(let i=0;i<d.length;i+=4) if(d[i]>220&&d[i+1]<90&&d[i+2]>170&&d[i+3]>200) n++; return n; })()''')
    chk(pink2 < 30 and '화소가 바뀌었습니다' in pg.inner_text('#imPanel'), '필터 뒤 옛 덧그림 숨김 + 알림', pink2)
    pg.evaluate('IMG_UI._undo()'); pg.wait_for_timeout(200)
    # 6) 아래 Drop 표 · CSV
    pg.click('#imBottom [data-b="drop"]'); pg.wait_for_selector('#drCsv')
    nrow = pg.locator('#imBbody tbody tr').count(); tot = pg.evaluate('IMG_UI._state().drops.length'); chk(nrow == tot and tot >= 20, 'Drop 표 줄 수', f'{nrow} / {tot}')
    with pg.expect_download() as d: pg.click('#drCsv')
    csv = open(d.value.path(), encoding='utf-8-sig').read(); chk(csv.startswith('Image,Frame,Mode,Method') and 'γ (mN/m)' in csv, 'CSV 머리', csv.split('\n')[0][:120])
    pg.screenshot(path=os.path.join(out, '04_table.png'))
    # 7) 매크로 기록
    steps = pg.evaluate('IMG_UI._steps().filter(s=>s.op==="drop").length'); chk(steps >= 7, 'Macro 에 Drop 단계 기록', steps)
    # 8) 새로 고침 후 복원
    pg.evaluate('IMG_UI._flush()'); pg.wait_for_timeout(500); pg.reload(); pg.wait_for_timeout(900)
    pg.click('#nav button[data-v="imgx"]'); pg.wait_for_function('window.IMG_UI && IMG_UI._state().images.length>=5', timeout=30000); pg.wait_for_timeout(500)
    keep = pg.evaluate('[IMG_UI._state().drops.length, IMG_UI._state().ui.drop.mode, IMG_UI._state().ui.drop.needleMM]')
    chk(keep[0] == tot and keep[1] == 'sessile' and abs(float(keep[2]) - t['needle_mm']) < 1e-9, '새로 고침 뒤 표·설정 그대로', keep)
    # 9) 실측 OpenDrop (있으면)
    if os.path.exists(os.path.join(RD, 'water_in_air001.png')):
        pg.set_input_files('#imFin', [os.path.join(RD, 'water_in_air001.png'), os.path.join(RD, 'drop_on_surface.png')])
        pg.wait_for_function('IMG_UI._state().images.length>=7', timeout=30000)
        sel('water_in_air001.png'); pg.click('#imTabs button[data-p="drop"]'); pg.click('#dpMode button[data-m="pendant"]'); pg.wait_for_selector('#dpNd')
        pg.fill('#dpNd', '0.7176'); pg.dispatch_event('#dpNd', 'change'); pg.click('#dpPhys summary') if not pg.is_visible('#dpG') else None
        pg.fill('#dpG', '9.8'); pg.dispatch_event('#dpG', 'change'); pg.fill('#dpRi', '1000'); pg.dispatch_event('#dpRi', 'change'); pg.fill('#dpRo', '0'); pg.dispatch_event('#dpRo', 'change')
        r = run()[0]; chk(abs(r['γ'] / 73.31 - 1) < 0.015, 'OpenDrop water_in_air001 γ ↔ 73.31 mN/m (±1.5 %)', f"{r['γ']:.2f}")
        pg.screenshot(path=os.path.join(out, '05_opendrop_pendant.png'))
        sel('drop_on_surface.png'); pg.click('#dpMode button[data-m="sessile"]'); pg.wait_for_selector('#dpRun')
        pg.evaluate("r=>IMG_UI._setActive(r)", {'id': 'R2', 'type': 'rect', 'pts': [[100, 330], [621, 761]]}); pg.wait_for_timeout(150)
        rows = run(); tp = [q for q in rows if q['Method'] == 'Tangent (polynomial)'][0]
        chk(abs(tp['θ left'] - 120.0) < 4 and abs(tp['θ right'] - 124.2) < 4, 'OpenDrop drop_on_surface 접선 ↔ 120.0 / 124.2° (±4°)', f"{tp['θ left']:.1f} / {tp['θ right']:.1f}")
        pg.screenshot(path=os.path.join(out, '06_opendrop_sessile.png'))
    # 10) 휴대폰 폭: 가로 넘침 없음
    pg.set_viewport_size({'width': 390, 'height': 844}); pg.wait_for_timeout(400)
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.ui.panel='drop'; IMG_UI.rerender(); })()"); pg.wait_for_timeout(300)
    ovf = pg.evaluate('document.documentElement.scrollWidth - window.innerWidth'); chk(ovf <= 1, '390 px 가로 넘침 없음', ovf)
    pg.screenshot(path=os.path.join(out, '07_mobile.png'), full_page=True)
    b.close()
chk(not errs, '콘솔 오류 0', errs[:5])
json.dump({'checks': checks, 'errors': errs}, open(os.path.join(out, 'ui_drop.json'), 'w'), ensure_ascii=False, indent=1)
print(f"\n{sum(1 for c in checks if c[0])}/{len(checks)} 통과")
sys.exit(1 if fails else 0)
