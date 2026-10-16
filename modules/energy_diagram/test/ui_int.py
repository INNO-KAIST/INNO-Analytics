import asyncio, json
from playwright.async_api import async_playwright
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await (await b.new_context(viewport={'width':1440,'height':1100})).new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e)))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(700); await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(700)
        await pg.click('#nav button[data-v="energy"]'); await pg.wait_for_timeout(300)
        await pg.fill('#edF','β-NaYF4:25%Yb,0.3%Tm'); await pg.click('#edApply'); await pg.wait_for_timeout(600)
        print('hidden?', await pg.evaluate("document.getElementById('edIntCard').hidden"))
        print(await pg.evaluate("document.querySelectorAll('#edIntCard .ed-int-dir').length"), (await pg.inner_text('#edIntCard'))[:900])
        await pg.click('#edIntCard button[data-int]'); await pg.wait_for_timeout(200)
        print('SEL', (await pg.inner_text('#edIntCard .ed-sug'))[:400])
        await pg.click('#edIntCard details summary'); await pg.fill('#edCda','1e-40'); await pg.fill('#edTauD','1'); await pg.wait_for_timeout(200)
        print('FR', await pg.inner_text('#edFr'))
        n0=await pg.evaluate("document.querySelectorAll('#edProcCard tbody tr').length"); await pg.click('#edIntCard button[data-intadd]'); await pg.wait_for_timeout(200)
        n1=await pg.evaluate("document.querySelectorAll('#edProcCard tbody tr').length"); print('proc rows',n0,n1)
        el=await pg.query_selector('#edIntCard'); await el.screenshot(path='/home/claude/ed/test/int.png')
        await pg.fill('#edF','β-NaYF4:2%Er'); await pg.click('#edApply'); await pg.wait_for_timeout(400); print('single ion hidden', await pg.evaluate("document.getElementById('edIntCard').hidden"))
        await pg.set_viewport_size({'width':390,'height':900}); await pg.fill('#edF','β-NaYF4:20%Yb,2%Er'); await pg.evaluate("document.getElementById('edApply').click()"); await pg.wait_for_timeout(500)
        print('mobile hscroll', await pg.evaluate("document.documentElement.scrollWidth>window.innerWidth+1"))
        print('errors',errs[:5]); await b.close()
asyncio.run(main())
