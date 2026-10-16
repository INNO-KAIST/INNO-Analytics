import asyncio, json
from playwright.async_api import async_playwright
U='/root/.claude/uploads/0c3ebbc7-7b2b-513c-a132-6d2362578318/'
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await (await b.new_context(viewport={'width':1440,'height':1000})).new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e))); pg.on('console',lambda m: errs.append(m.text) if m.type=='error' and 'ERR_TUNNEL' not in m.text else None)
        pg.on('dialog',lambda d: asyncio.ensure_future(d.accept()))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(700); await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(700)
        await pg.click('#nav button[data-v="dre"]'); await pg.click('#dreTabs button[data-t="data"]')
        await pg.set_input_files('#dreFile',U+'d1870d20-260702_power_scan_WJ-0128c_14_Oil_NIR.h5'); await pg.wait_for_timeout(2000); await pg.click('#hOK'); await pg.wait_for_timeout(400)
        st=lambda: pg.evaluate("DRE_UI._state().samples[0].data.map(d=>({calOK:d.calOK,I0:d.I[0],P0:d.raw&&d.raw.P[0],cal:d.cal&&(d.cal.name||d.cal.kPerMW),src:d.cal&&d.cal.src}))")
        print('after h5', json.dumps(await st(),ensure_ascii=False))
        await pg.click('#dreTabs button[data-t="data"]'); await pg.click('button[data-cal="0|0"]'); await pg.wait_for_timeout(300)
        await pg.select_option('#cPre','1'); await pg.wait_for_timeout(100); print('preview:', await pg.inner_text('#cPrev'))
        await pg.click('#cOne'); await pg.wait_for_timeout(300)
        s1=await st(); print('blue on up: ratio', s1[0]['I0']/s1[1]['I0']*(s1[1]['P0']/s1[0]['P0']) if s1[0]['P0'] else None, s1[0]['cal'])
        # 휠 재구성
        await pg.click('#dreTabs button[data-t="data"]'); await pg.click('button[data-cal="0|1"]'); await pg.wait_for_timeout(300)
        await pg.check('#cWh'); await pg.click('#cAuto'); await pg.wait_for_timeout(200); print('wheel fit:', await pg.inner_text('#cFit'))
        await pg.click('#cOne'); await pg.wait_for_timeout(300)
        r=await pg.evaluate("(()=>{const d=DRE_UI._state().samples[0].data[1]; const k=268.13e3*1e3; return d.I.map((v,i)=>v/(d.raw.P[i]*1e-3*2.6813e8)).filter((_,i)=>d.raw.P[i]>0.05*Math.max(...d.raw.P))})()")
        print('wheel vs PM ratio range', min(r), max(r))
        # CSV mW 가져오기 → 미확인 경고
        await pg.click('#dreTabs button[data-t="data"]')
        await pg.fill('#drePaste','P_mW,sig\n1,5\n2,10\n5,200\n10,5000\n20,20000\n'); await pg.click('#dreParse'); await pg.click('#mOK'); await pg.wait_for_timeout(300)
        await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(200)
        print('fit warning shown:', await pg.evaluate("/보정을 확인하지 않은 데이터 1개/.test(document.querySelector('#drePanel').innerText)"))
        await pg.click('#dreTabs button[data-t="data"]'); await pg.click('button[data-cal="0|2"]'); await pg.wait_for_timeout(200); await pg.click('#cAll'); await pg.wait_for_timeout(300)
        await pg.click('#dreTabs button[data-t="fit"]'); await pg.wait_for_timeout(200)
        print('warning after apply-all:', await pg.evaluate("/보정을 확인하지 않은/.test(document.querySelector('#drePanel').innerText)"))
        print('final', json.dumps(await st(),ensure_ascii=False))
        # 프로필 저장
        await pg.click('#dreTabs button[data-t="data"]'); await pg.click('button[data-cal="0|0"]'); await pg.wait_for_timeout(200)
        await pg.fill('#cC','3.1'); await pg.fill('#cA','1.3'); await pg.fill('#cNm','테스트 프로필'); await pg.click('#cSave'); await pg.wait_for_timeout(300)
        print('profile saved', await pg.evaluate("JSON.parse(localStorage.getItem('innoanalytics_dre_cal_v1')||'[]').map(p=>p.name)"))
        await pg.click('#dreTabs button[data-t="data"]'); await pg.click('button[data-cal="0|0"]'); await pg.wait_for_timeout(200)
        el=await pg.query_selector('#dreMap'); await el.screenshot(path='/home/claude/dre/shots/cal.png')
        print('errors',errs[:5]); await b.close()
asyncio.run(main())
