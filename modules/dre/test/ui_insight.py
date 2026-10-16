import asyncio, json, math
from playwright.async_api import async_playwright
OUT='/home/claude/dre/shots'
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); ctx=await b.new_context(viewport={'width':1440,'height':950}); pg=await ctx.new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append('pageerror '+str(e)))
        pg.on('console',lambda m: errs.append('console '+m.text) if m.type=='error' and 'ERR_TUNNEL' not in m.text else None)
        pg.on('dialog',lambda d: asyncio.ensure_future(d.accept()))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(500)
        await pg.evaluate("localStorage.removeItem('innoanalytics_dre_v1')"); await pg.reload(); await pg.wait_for_timeout(400)
        await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(400)
        await pg.click('#dreTabs button[data-t="data"]')
        await pg.click('#dreSynth'); await pg.click('#sOK'); await pg.wait_for_timeout(300)
        # 트레이스 붙여넣기 (진단용): 100 kW/cm² 에서 단조 상승 + 인위적 과도 상승
        rows=['t_ms,y']+[f'{t*0.05:.3f},{1000*(1-math.exp(-t*0.05/0.3))*(1+0.12*math.exp(-t*0.05/8)*(1-math.exp(-t*0.05/1))):.2f}' for t in range(0,1200)]
        await pg.select_option('#dreDType','trace'); await pg.fill('#drePaste','\n'.join(rows)); await pg.click('#dreParse')
        await pg.select_option('#mTU','1e-3'); await pg.fill('#mTP','100'); await pg.select_option('#mPU','kW/cm²'); await pg.click('#mOK'); await pg.wait_for_timeout(300)
        await pg.screenshot(path=f'{OUT}/i0_data.png',full_page=True)
        await pg.click('#dreTabs button[data-t="fit"]')
        await pg.click('#fAdq'); await pg.wait_for_timeout(2500)
        n_adq=await pg.evaluate("document.querySelectorAll('#dreDiag .dre-adq tr').length"); print('adequacy-now rows',n_adq)
        await pg.fill('#fS','8'); await pg.fill('#fN','7'); await pg.press('#fN','Tab')
        await pg.click('#fRun')
        for i in range(300):
            await pg.wait_for_timeout(1000)
            if await pg.evaluate("!!(DRE_UI._state().result)"): break
        print('fit done after',i,'s')
        await pg.wait_for_timeout(600)
        info=await pg.evaluate("(()=>{const R=DRE_UI._state().result; return {adq:!!R.adequacy, worst:R.adequacy&&R.adequacy.worst, items:R.adequacy&&R.adequacy.items.map(q=>q.key+':'+q.level), ins:!!R.insight, secs:document.querySelectorAll('#dreDiag details.dre-ins').length}})()")
        print(json.dumps(info,ensure_ascii=False))
        await pg.evaluate("document.querySelectorAll('#dreDiag details.dre-ins').forEach(d=>d.open=true)"); await pg.wait_for_timeout(300)
        await pg.select_option('#insProfJ','1'); await pg.wait_for_timeout(300)
        await pg.evaluate("document.querySelectorAll('#dreDiag details.dre-ins').forEach(d=>d.open=true)")
        await pg.click('#slGo')
        for i in range(120):
            await pg.wait_for_timeout(1000)
            if await pg.evaluate("!!DRE_UI._state().insSlice"): break
        print('slice after',i,'s')
        await pg.evaluate("document.querySelectorAll('#dreDiag details.dre-ins').forEach(d=>d.open=true)"); await pg.wait_for_timeout(300)
        await pg.screenshot(path=f'{OUT}/i1_diag.png',full_page=True)
        el=await pg.query_selector('details[data-k="pred"]')
        print('derived rows', await pg.evaluate("(DRE_UI._state().result.derived||[]).length"))
        if el: await el.screenshot(path=f'{OUT}/pred.png')
        ov=await pg.evaluate("(()=>{const o=[];document.querySelectorAll('#dreDiag *').forEach(e=>{const r=e.getBoundingClientRect(); const p=document.querySelector('#dreDiag').getBoundingClientRect(); if(r.right>p.right+2&&!e.closest('.dre-tw')) o.push(e.tagName+'.'+e.className)}); return o.slice(0,8)})()")
        print('overflow',ov)
        await pg.evaluate('toggleTheme()'); await pg.wait_for_timeout(600)
        await pg.set_viewport_size({'width':390,'height':900}); await pg.wait_for_timeout(500)
        await pg.screenshot(path=f'{OUT}/i2_mobile.png',full_page=True)
        hs=await pg.evaluate("document.documentElement.scrollWidth>window.innerWidth+1"); print('mobile hscroll',hs)
        print('errors',errs[:10])
        await b.close()
asyncio.run(main())
