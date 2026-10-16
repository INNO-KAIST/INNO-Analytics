import asyncio, json
from playwright.async_api import async_playwright
async def run(pg,w):
    await pg.evaluate(f"(()=>{{const S=DRE_UI._state(); S.fit.workers={w}; S.fit.starts=16; S.fit.npts=7; S.fit.seed=3; S.result=null; DRE_UI.rerender();}})()")
    await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(200)
    t=await pg.evaluate("Date.now()"); await pg.click('#fRun')
    for i in range(600):
        await pg.wait_for_timeout(500)
        if await pg.evaluate("!!DRE_UI._state().result"): break
    r=await pg.evaluate("(()=>{const R=DRE_UI._state().result; return {chi:R.ms.map(x=>+x.chi2.toFixed(6)), best:R.best.theta.map(v=>+v.toFixed(6)), prof:R.prof.map(p=>p?[+p.lo.toFixed(5),+p.hi.toFixed(5)]:null), time:R.time, workers:R.workers, adq:!!R.adequacy, der:(R.derived||[]).length}})()")
    return r
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(args=['--enable-features=SharedArrayBuffer']); ctx=await b.new_context(viewport={'width':1440,'height':950}); pg=await ctx.new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e))); pg.on('dialog',lambda d: asyncio.ensure_future(d.accept()))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(600)
        await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(600)
        await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(300)
        await pg.click('#dreTabs button[data-t="data"]'); await pg.click('#dreSynth'); await pg.click('#sOK'); await pg.wait_for_timeout(300)
        print('hardwareConcurrency', await pg.evaluate("navigator.hardwareConcurrency"))
        await pg.evaluate("(()=>{const S=DRE_UI._state(); S.fit.workers=0; S.fit.starts=16; S.fit.npts=7; S.fit.seed=3;})()")
        await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(200)
        r0=await run(pg,0)
        cal=await pg.evaluate("JSON.parse(localStorage.getItem('innoanalytics_dre_hw_v1')||'null')")
        print('calibration',json.dumps({k:cal[k] for k in ['hc','best','phys','smt']} if cal else None), [ (r['K'],round(r['speedup'],2)) for r in (cal['rows'] if cal else [])])
        print('auto run workers',r0['workers'],'time',r0['time'])
        r1=await run(pg,1); r4=await run(pg,4)
        print('auto==w1', r0['chi']==r1['chi'] and r0['prof']==r1['prof'])
        print('w1 time',r1['time'],'workers',r1['workers'],'| w4 time',r4['time'],'workers',r4['workers'])
        print('ms chi2 identical', r1['chi']==r4['chi'], 'best identical', r1['best']==r4['best'], 'profiles identical', r1['prof']==r4['prof'], 'adq/derived', r4['adq'], r4['der'])
        # 중단 버튼
        await pg.click('#fRun'); await pg.wait_for_timeout(700); await pg.click('#fStop'); await pg.wait_for_timeout(300)
        print('after stop running?', await pg.evaluate("!!document.querySelector('#fStop:not([disabled])')"))
        await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(300)
        await pg.evaluate("document.querySelectorAll('#dreDiag details.dre-ins').forEach(d=>d.open=true)"); await pg.click('#slGo')
        for i in range(120):
            await pg.wait_for_timeout(500)
            if await pg.evaluate("!!DRE_UI._state().insSlice"): break
        print('slice ok', await pg.evaluate("(()=>{const q=DRE_UI._state().insSlice; return q&&q.z.length===17&&q.z.every(r=>r.length===17)})()"))
        el=await pg.query_selector('#fHw'); 
        if el: await el.screenshot(path='/home/claude/dre/shots/hw.png')
        print('errors',errs[:5]); await b.close()
asyncio.run(main())
