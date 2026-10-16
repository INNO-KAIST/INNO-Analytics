import asyncio, json, time
from playwright.async_api import async_playwright
D='/home/claude/dre/test/dy033/'
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await (await b.new_context(viewport={'width':1440,'height':1000})).new_page()
        pg.on('dialog',lambda d: asyncio.ensure_future(d.accept()))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(700); await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(700)
        await pg.click('#nav button[data-v="dre"]'); await pg.set_input_files('#dreImp',D+'model.json'); await pg.wait_for_timeout(800)
        await pg.click('#dreTabs button[data-t="data"]')
        for _ in range(4): await pg.click('button[data-fix$="|del"]'); await pg.wait_for_timeout(150)
        await pg.click('button[data-fix$="|split"]'); await pg.wait_for_timeout(200)
        for _ in range(2):
            bb=await pg.query_selector('button[data-fix$="|bg"]')
            if bb: await bb.click(); await pg.wait_for_timeout(200)
        # up branch 만 남긴다
        await pg.evaluate("(()=>{const S=DRE_UI._state(); S.samples[0].data=S.samples[0].data.filter(d=>d.dir>0); S.samples[0].data[0].relErr=0.02; S.fit.starts=24; S.fit.profile=false; S.fit.cpu=100; S.result=null; DRE_UI.rerender();})()")
        await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(200); t=time.time(); await pg.click('#fRun')
        for i in range(1800):
            await pg.wait_for_timeout(1000)
            if await pg.evaluate("!!DRE_UI._state().result"): break
        r=await pg.evaluate("""(()=>{const R=DRE_UI._state().result; if(!R) return null; return {s2:R.info.s2, time:R.time, adq:R.adequacy.items.map(q=>[q.key,q.level,q.summary]), free:R.free.map((f,j)=>[f.key,+R.best.theta[j].toFixed(2),f.lo,f.hi]), ms:R.ms.map(m=>+m.chi2.toFixed(0)).sort((a,b)=>a-b).slice(0,8)}})()""")
        print(json.dumps(r,ensure_ascii=False,indent=0)); print('wall',time.time()-t); await b.close()
asyncio.run(main())
