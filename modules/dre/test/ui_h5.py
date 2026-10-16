import asyncio, json
from playwright.async_api import async_playwright
U='/root/.claude/uploads/0c3ebbc7-7b2b-513c-a132-6d2362578318/'
F=[U+'d1870d20-260702_power_scan_WJ-0128c_14_Oil_NIR.h5',U+'c013fe13-260702_rise_time_WJ-0128c_800_nm_at_Threshold.h5']
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await (await b.new_context(viewport={'width':1440,'height':950})).new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e))); pg.on('console',lambda m: errs.append(m.text) if m.type=='error' and 'ERR_TUNNEL' not in m.text else None)
        pg.on('dialog',lambda d: asyncio.ensure_future(d.accept()))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(800)
        await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(800)
        await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(300)
        await pg.click('#dreTabs button[data-t="data"]')
        await pg.set_input_files('#dreFile',F[0]); await pg.wait_for_timeout(2500)
        print('scan UI:', (await pg.inner_text('#dreMap'))[:400].replace('\n',' | '))
        await pg.click('#hOK'); await pg.wait_for_timeout(500)
        await pg.click('#dreTabs button[data-t="data"]')
        await pg.set_input_files('#dreFile',F[1]); await pg.wait_for_timeout(2500)
        print('trace UI:', (await pg.inner_text('#dreMap'))[:400].replace('\n',' | '))
        await pg.click('#hOK'); await pg.wait_for_timeout(500)
        d=await pg.evaluate("DRE_UI._state().samples[0].data.map(d=>({type:d.type,dir:d.dir,n:(d.I.length>1?d.I.length:d.t.length),I0:d.I[0],I1:d.I[d.I.length-1],y0:(d.y||[])[0],label:d.label}))")
        print(json.dumps(d,ensure_ascii=False))
        await pg.click('#dreTabs button[data-t="data"]'); await pg.wait_for_timeout(300)
        await pg.screenshot(path='/home/claude/dre/shots/h5.png',full_page=True)
        # ED 루프 중복·삭제
        await pg.click('#nav button[data-v="energy"]'); await pg.wait_for_timeout(800)
        names=await pg.evaluate("ED_UI.PRESETS.map(p=>p.name)"); i=[k for k,n in enumerate(names) if 'Tm' in n and '1064' in n][0]; print('preset',names[i])
        await pg.evaluate(f"ED_UI._loadPreset({i})"); await pg.evaluate("ED_UI._state().procs=[]; ED_UI.rerender()"); await pg.wait_for_timeout(500)
        btn=await pg.query_selector('button[data-pa]')
        print('pa button', bool(btn))
        if btn:
            await btn.click(); await pg.wait_for_timeout(400)
            n1=await pg.evaluate("ED_UI._state?ED_UI._state().procs.length:null")
            rm=await pg.query_selector('button[data-parm]'); print('procs after add',n1,'remove button',bool(rm))
            # 같은 루프를 다른 경로로 다시 추가 시도: data-pa 버튼이 사라졌는지
            again=await pg.query_selector('button[data-pa]'); print('pa button still present (other loops)',bool(again))
            # 같은 과정을 수동으로 하나 남긴 뒤 다시 추가 → 중복 건너뜀
            await rm.click(); await pg.wait_for_timeout(300); print('procs after remove',await pg.evaluate("ED_UI._state().procs.length"))
            await (await pg.query_selector('button[data-pa]')).click(); await pg.wait_for_timeout(300)
            await pg.evaluate("(()=>{const S=ED_UI._state(); const p=Object.assign({},S.procs[0]); delete p.grp; p.uid='x1'; S.procs=[p]; ED_UI.rerender();})()"); await pg.wait_for_timeout(300)
            await (await pg.query_selector('button[data-pa]')).click(); await pg.wait_for_timeout(300)
            print('after re-add with 1 pre-existing:',await pg.evaluate("ED_UI._state().procs.length"),'(expect',n1,')')
            await pg.screenshot(path='/home/claude/dre/shots/ed_procs.png',full_page=False)
            await pg.click('#edProcClr'); await pg.wait_for_timeout(300); print('after clear',await pg.evaluate("ED_UI._state().procs.length"))
        print('errors',errs[:6]); await b.close()
asyncio.run(main())
