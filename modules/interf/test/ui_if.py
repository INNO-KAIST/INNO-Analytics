import asyncio, json
from playwright.async_api import async_playwright
T='/home/claude/modules_if/test/'
Q="""(()=>{const R=IF_SUG._last(); if(!R) return null; return {mode:R.mode, fit:R.fit?R.fit.cols.map(c=>[c.el+' '+c.fam, Math.round(c.A)]):null,
 out:R.out.map(o=>({el:o.el, best:o.best?[o.best.line.nm,o.best.line.fam,+o.best.line.e.toFixed(3),+o.best.err.toFixed(4),+o.best.purity.toFixed(3),o.best.top&&o.best.top.line?o.best.top.line.el+' '+o.best.top.line.nm:null]:null,
 alt:o.alt?[o.alt.line.nm,+o.alt.err.toFixed(4),+o.alt.purity.toFixed(3)]:null, fu:o.famUnsure, ux:!!(o.best&&o.best.unexpl)}))}})()"""
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); ctx=await b.new_context(viewport={'width':1440,'height':1000}); pg=await ctx.new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e))); pg.on('console',lambda m: m.type=='error' and 'ERR_' not in m.text and errs.append(m.text))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(800); await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(800)
        # 1) 파일 없이 — 모형 모드
        await pg.click('#nav button[data-v="interf"]'); await pg.wait_for_timeout(300)
        await pg.fill('#ifEls','Na, F, Lu, Yb, Tm, Cu, C, O'); await pg.fill('#ifComp','Na 1, F 4, Lu 0.6, Yb 0.38, Tm 0.02, Cu 0.05, C 0.5, O 0.1'); await pg.click('#ifRun'); await pg.wait_for_timeout(400)
        r=await pg.evaluate(Q); print('MODEL', json.dumps(r,ensure_ascii=False))
        print('details closed:', await pg.evaluate("(()=>{const d=document.querySelector('details.if-all'); return d?[d.open, d.querySelector('summary').textContent.trim().slice(0,40)]:null})()"))
        await pg.screenshot(path=T+'model.png',full_page=False)
        # 2) 파일 열기 — 실측 모드
        await pg.set_input_files('#fin',T+'synthetic_eds.emd'); await pg.wait_for_timeout(3000)
        await pg.click('#nav button[data-v="interf"]'); await pg.wait_for_timeout(300)
        await pg.fill('#ifEls','Na, F, Lu, Yb, Tm, Cu, C, O'); await pg.click('#ifRun'); await pg.wait_for_timeout(1500)
        r=await pg.evaluate(Q); print('MEAS', json.dumps(r,ensure_ascii=False))
        el=await pg.query_selector('#ifSug'); await el.screenshot(path=T+'meas.png')
        # 3) Cu 빼고 — 미설명 신호 경고?
        await pg.fill('#ifEls','Na, F, Lu, Yb, Tm, C, O'); await pg.click('#ifRun'); await pg.wait_for_timeout(1200)
        r=await pg.evaluate(Q); print('NO-CU', json.dumps([o for o in r['out'] if o['ux']],ensure_ascii=False))
        # 4) 원소 하나
        await pg.fill('#ifEls','Yb'); await pg.click('#ifRun'); await pg.wait_for_timeout(800); r=await pg.evaluate(Q); print('ONE', json.dumps(r['out'],ensure_ascii=False))
        # 5) 다크 + 모바일
        await pg.set_viewport_size({'width':390,'height':900}); await pg.fill('#ifEls','Na, F, Lu, Yb, Tm'); await pg.evaluate("document.getElementById('ifRun').click()"); await pg.wait_for_timeout(800)
        print('mobile hscroll', await pg.evaluate("document.documentElement.scrollWidth>window.innerWidth+1"))
        print('errors', errs[:5]); await b.close()
asyncio.run(main())
