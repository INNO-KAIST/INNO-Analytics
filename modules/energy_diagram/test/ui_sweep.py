import asyncio
from playwright.async_api import async_playwright
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await b.new_page(viewport={'width':1440,'height':950})
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e)))
        await pg.goto('file:///home/claude/ed/out.html'); await pg.wait_for_timeout(400)
        await pg.click('#nav button[data-v="energy"]')
        ions=await pg.evaluate("Object.keys(ED_CORE.IONS)")
        for w in (1440,768,390):
            await pg.set_viewport_size({'width':w,'height':900})
            for theme in (0,1):
                for i in range(0,len(ions),3):
                    grp=ions[i:i+3]
                    await pg.evaluate("(g)=>{const S=ED_UI._state(); S.cols=[]; S.procs=[]; S.formula=''; S.lasers=[365,808,980,1064]; S.lasersUser=true; ED_UI.rerender(); g.forEach(k=>{ S.cols.push({uid:Math.random().toString(36).slice(2),ion:k,conc:[null],layers:['-'],notes:[]}); }); ED_UI.rerender(); }",grp)
                    await pg.click('#edSugRun'); await pg.wait_for_timeout(30)
                    ov=await pg.evaluate("document.documentElement.scrollWidth>window.innerWidth+1")
                    txt=await pg.evaluate("document.getElementById('view').innerText.includes('렌더 오류')")
                    if ov or txt: errs.append(f'{w} {theme} {grp} overflow={ov} rendererr={txt}')
                await pg.click('#chipTheme')
        for i in range(10):
            await pg.select_option('#edPreset',str(i)); await pg.wait_for_timeout(30)
        print('errors',errs or 'none')
        await b.close()
asyncio.run(main())
