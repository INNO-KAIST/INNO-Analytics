import asyncio, json
from playwright.async_api import async_playwright
D='/home/claude/dre/test/dy033/'
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await (await b.new_context(viewport={'width':1440,'height':1000})).new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e))); pg.on('dialog',lambda d: asyncio.ensure_future(d.accept()))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(700); await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(700)
        await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(300)
        # 1) 연구원 JSON 불러오기 → 점검 표시
        await pg.set_input_files('#dreImp',D+'model.json'); await pg.wait_for_timeout(800)
        await pg.click('#dreTabs button[data-t="data"]'); await pg.wait_for_timeout(300)
        iss=await pg.evaluate("Array.from(document.querySelectorAll('tr.dre-iss .note')).map(n=>n.textContent.trim().slice(0,90))"); print('ISSUES',json.dumps(iss,ensure_ascii=False,indent=0))
        await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(200)
        print('FITNOTE', (await pg.evaluate("(document.querySelector('#dreBody .note.warn')||{}).textContent||''"))[:200])
        await pg.click('#fRun'); await pg.wait_for_timeout(400); print('RUN blocked?', await pg.evaluate("!!document.querySelector('#fStop:not([disabled])')")==False)
        # 고치기: rise 4개 삭제, split, bg
        await pg.click('#dreTabs button[data-t="data"]'); await pg.wait_for_timeout(200)
        for _ in range(4):
            await pg.click('button[data-fix$="|del"]'); await pg.wait_for_timeout(200)
        await pg.click('button[data-fix$="|split"]'); await pg.wait_for_timeout(300)
        for _ in range(2):
            bb=await pg.query_selector('button[data-fix$="|bg"]')
            if bb: await bb.click(); await pg.wait_for_timeout(300)
        st=await pg.evaluate("DRE_UI._state().samples[0].data.map(d=>({type:d.type,n:d.I.length,dir:d.dir,sbg:d.sbg,bg:d.bgSub, y0:d.y.slice(0,3), iss:(d.type)}))"); print('AFTER FIX',json.dumps(st,ensure_ascii=False))
        iss=await pg.evaluate("Array.from(document.querySelectorAll('tr.dre-iss .note')).map(n=>n.textContent.trim().slice(0,80))"); print('ISSUES2',json.dumps(iss,ensure_ascii=False))
        # 2) 새 시료에 CSV 넣기 (자동 분리·σ·Interp·배경)
        await pg.click('#dreAddSmp'); await pg.wait_for_timeout(200); await pg.select_option('#dreDSmp','1')
        await pg.select_option('#dreDType','scurve'); await pg.set_input_files('#dreFile',D+'power_scan_DY033_oil_step0.5_rep3_ND2.csv'); await pg.wait_for_timeout(500)
        mp=await pg.inner_text('#dreMap'); print('MAP',mp[:700])
        await pg.select_option('#mPU','I|kW/cm²'); await pg.click('#mOK'); await pg.wait_for_timeout(400)
        st=await pg.evaluate("DRE_UI._state().samples[1].data.map(d=>({n:d.I.length,dir:d.dir,sbg:+(d.sbg||0).toFixed(2),syNote:d.syNote,nInterp:d.nInterp,bg:d.bgSub,I0:d.I[0],sy0:d.sy&&d.sy.slice(0,2)}))"); print('CSV',json.dumps(st,ensure_ascii=False))
        # 3) pulse CSV
        await pg.select_option('#dreDSmp','1'); await pg.set_input_files('#dreFile',D+'pulse_8uW.csv'); await pg.wait_for_timeout(600)
        print('PULSE',(await pg.inner_text('#dreMap'))[:600])
        await pg.click('#pOK'); await pg.wait_for_timeout(300)
        st=await pg.evaluate("(()=>{const d=DRE_UI._state().samples[1].data.slice(-1)[0]; return {type:d.type,n:d.t.length,I:d.I,t1:d.t[d.t.length-1],stats:d.stats}})()"); print('TRACE',json.dumps(st))
        print('errors',errs[:5]); await b.close()
asyncio.run(main())
