"""탭 막대 맞춤(nav_fit.js) · EM Image 안 Metadata 보기 — Chromium 시험
사용: python3 modules/nav/ui_nav.py "INNO Analytics.html" <결과 폴더> [dm3 파일]"""
import sys, os
from playwright.sync_api import sync_playwright
html, out = sys.argv[1], sys.argv[2]; dm3 = sys.argv[3] if len(sys.argv) > 3 else None; os.makedirs(out, exist_ok=True)
errs = []; ok = []
def chk(c, n, i=''):
    ok.append(bool(c)); print(('  ✓ ' if c else '  ✗ ') + n + (f' — {i}' if i != '' else ''))
JS = '''(()=>{const n=document.querySelector('#nav'); const vis=[...n.querySelectorAll('button')].filter(b=>getComputedStyle(b).display!=='none'); return {rows:new Set(vis.map(b=>Math.round(b.getBoundingClientRect().top))).size,short:n.classList.contains('nav-short'),more:!!n.querySelector('.nav-more'),hid:[...n.querySelectorAll('.nav-hid')].map(b=>b.dataset.v),over:n.scrollWidth-n.clientWidth,keys:[...n.querySelectorAll('button[data-v]')].map(b=>b.dataset.v),active:(n.querySelector('button.active')||{}).dataset}})()'''
with sync_playwright() as p:
    b = p.chromium.launch()
    for w in [360, 390, 720, 721, 820, 1024, 1280, 1440, 1920]:
        pg = b.new_page(viewport={'width': w, 'height': 800}, has_touch=w <= 720, is_mobile=w <= 720)
        pg.on('pageerror', lambda e: errs.append(str(e)))
        pg.goto('file://' + os.path.abspath(html)); pg.wait_for_timeout(500)
        for v in ['interf', 'dre']:
            pg.click(f'#nav button[data-v="{v}"]'); pg.wait_for_timeout(200); i = pg.evaluate(JS)
            chk(i['rows'] == 1 and i['over'] <= 1, f'{w} px · {v}: 한 줄, 넘침 없음', i)
        if w == 360: chk('meta' not in i['keys'] and len(i['keys']) == 9, '탭 9 개 (EM Metadata 빠짐)', i['keys'])
        pg.close()
    # 아주 좁은 막대(확대 200 % 등) → More 메뉴
    pg = b.new_page(viewport={'width': 1280, 'height': 800}); pg.on('pageerror', lambda e: errs.append(str(e)))
    pg.goto('file://' + os.path.abspath(html)); pg.wait_for_timeout(500); pg.click('#nav button[data-v="dre"]'); pg.wait_for_timeout(200)
    pg.evaluate("document.querySelector('#nav').style.maxWidth='420px'"); pg.wait_for_timeout(300); i = pg.evaluate(JS)
    chk(i['more'] and i['short'] and i['rows'] == 1 and i['over'] <= 1 and 'dre' not in i['hid'] and len(i['hid']) >= 2, '좁으면 More ▾ · 지금 탭(DRE)은 보임', i)
    pg.click('#nav .nav-more'); pg.wait_for_timeout(150)
    items = pg.evaluate("[...document.querySelectorAll('.nav-menu button')].map(b=>b.dataset.v)"); chk(items == i['hid'], 'More 메뉴 = 접힌 탭', items)
    pg.screenshot(path=os.path.join(out, 'nav_more.png'), clip={'x': 300, 'y': 60, 'width': 700, 'height': 360})
    pg.keyboard.press('Escape'); chk(pg.locator('.nav-menu').count() == 0, 'Esc 로 닫힘')
    pg.click('#nav .nav-more'); pg.click('.nav-menu button[data-v="interf"]'); pg.wait_for_timeout(250); i2 = pg.evaluate(JS)
    chk(i2['active'].get('v') == 'interf' and 'interf' not in i2['hid'], '메뉴에서 고른 탭으로 가고 그 탭은 막대에 보임', i2)
    pg.evaluate("document.querySelector('#nav').style.maxWidth=''"); pg.wait_for_timeout(300); i3 = pg.evaluate(JS)
    chk(not i3['more'] and not i3['hid'], '넓어지면 More 사라짐', i3)
    # EM Image 안 Metadata
    if dm3:
        pg.set_input_files('#fin', [dm3]); pg.wait_for_function('document.querySelector("#nav button[data-v=image]") && !document.querySelector("#nav button[data-v=image]").disabled', timeout=30000); pg.wait_for_timeout(600)
        pg.click('#nav button[data-v="image"]'); pg.wait_for_selector('.emsub')
        chk(pg.locator('.emsub button.active').inner_text() == 'Image', 'EM Image 위에 Image · Metadata 전환, 처음은 Image')
        pg.click('.emsub button[data-sv="meta"]'); pg.wait_for_selector('#mdBody')
        n = pg.locator('#mdBody tr').count(); act = pg.evaluate("document.querySelector('#nav button.active').dataset.v")
        chk(n > 5 and act == 'image', 'Metadata 보기: 표 · 탭 막대는 EM Image 가 선택', f'{n} 줄, active {act}')
        pg.fill('#mdQ', 'Pixel'); pg.wait_for_timeout(100); chk(0 < pg.locator('#mdBody tr').count() < n, '메타데이터 찾기')
        pg.screenshot(path=os.path.join(out, 'em_meta.png'))
        pg.click('.emsub button[data-sv="image"]'); pg.wait_for_timeout(400); chk(pg.locator('#mdBody').count() == 0 and pg.locator('.emsub button.active').inner_text() == 'Image', '다시 Image')
    pg.close(); b.close()
chk(not errs, '페이지 오류 0', errs[:3])
print(f'\n{sum(ok)}/{len(ok)} 통과'); sys.exit(0 if all(ok) else 1)
