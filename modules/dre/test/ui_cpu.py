import asyncio, json
from playwright.async_api import async_playwright
Q="(()=>{const R=DRE_UI._state().result; return {chi:R.ms.map(x=>+x.chi2.toFixed(6)), prof:R.prof.map(p=>p?[+p.lo.toFixed(5),+p.hi.toFixed(5)]:null), workers:R.workers}})()"
async def wait(pg):
    for i in range(600):
        await pg.wait_for_timeout(300)
        if await pg.evaluate("!!DRE_UI._state().result"): return
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); ctx=await b.new_context(viewport={'width':1440,'height':950})
        await ctx.add_init_script("Object.defineProperty(navigator,'hardwareConcurrency',{get:()=>8})")
        pg=await ctx.new_page(); errs=[]; pg.on('pageerror',lambda e:errs.append(str(e))); pg.on('dialog',lambda d: asyncio.ensure_future(d.accept()))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(600); await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(600)
        # 가짜 측정 결과: 8스레드에서 best 8
        await pg.evaluate("""localStorage.setItem('innoanalytics_dre_hw_v1',JSON.stringify({hc:8,ua:(navigator.userAgentData&&navigator.userAgentData.platform)||navigator.platform||'',mem:null,cores:null,best:8,phys:4,smt:1.2,date:new Date().toISOString(),top:10,rows:[{K:1,thr:2,speedup:1,eff:1},{K:2,thr:4,speedup:2,eff:1},{K:4,thr:7.6,speedup:3.8,eff:.95},{K:6,thr:8.6,speedup:4.3,eff:.72},{K:8,thr:10,speedup:5,eff:.62}]}))""")
        await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(300)
        await pg.click('#dreTabs button[data-t="data"]'); await pg.click('#dreSynth'); await pg.click('#sOK'); await pg.wait_for_timeout(300)
        await pg.evaluate("(()=>{const S=DRE_UI._state(); S.fit.workers=0; S.fit.starts=24; S.fit.npts=7; S.fit.seed=3; S.result=null; delete S.fit.cpu; DRE_UI.rerender();})()")
        await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(200)
        print('default:', await pg.eval_on_selector('#fCpu','e=>e.value'), '|', await pg.inner_text('#fCpuN'))
        for v in ['100','90','50','25']:
            await pg.select_option('#fCpu',v); await pg.wait_for_timeout(50); print(v,'->',await pg.inner_text('#fCpuN'))
        # 기준: 워커 1개
        await pg.evaluate("(()=>{const S=DRE_UI._state(); S.fit.workers=1; S.result=null; DRE_UI.rerender();})()"); await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(200)
        await pg.click('#fRun'); await wait(pg); r1=await pg.evaluate(Q)
        # 자동 + 실행 중 한도 변경 25 → 100 → 50
        await pg.evaluate("(()=>{const S=DRE_UI._state(); S.fit.workers=0; S.result=null; DRE_UI.rerender();})()"); await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(200)
        await pg.select_option('#fCpu','25'); await pg.click('#fRun'); await pg.wait_for_timeout(1500)
        print('running txt', await pg.inner_text('#fTxt'))
        await pg.select_option('#fCpu','100'); await pg.wait_for_timeout(1500); print('after→100', await pg.inner_text('#fTxt'), 'pool', await pg.evaluate("document.querySelector('#fTxt')&&0"))
        await pg.select_option('#fCpu','50'); await wait(pg); r2=await pg.evaluate(Q)
        print('identical chi2', r1['chi']==r2['chi'], 'identical profiles', r1['prof']==r2['prof'], 'workers r1', r1['workers'], 'r2 max', r2['workers'])
        print('saved cpu', await pg.evaluate("DRE_UI._state().fit.cpu"))
        print('errors',errs[:5]); await b.close()
asyncio.run(main())
