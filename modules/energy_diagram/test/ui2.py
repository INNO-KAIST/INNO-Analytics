import asyncio, json
from playwright.async_api import async_playwright
OUT='/home/claude/ed/shots'
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch()
        pg=await b.new_page(viewport={'width':1440,'height':950})
        errs=[]
        pg.on('pageerror',lambda e:errs.append('pageerror '+str(e)))
        pg.on('console',lambda m: errs.append('console '+m.text) if m.type=='error' and 'ERR_TUNNEL' not in m.text else None)
        async def _ce(msg):
            if msg.type=='error' and 'ERR_TUNNEL' not in msg.text:
                try:
                    vals=[await a.evaluate('x=>x&&x.stack||String(x)') for a in msg.args]; errs.append('STACK '+' | '.join(vals))
                except Exception as e: pass
        pg.on('console',lambda m: asyncio.ensure_future(_ce(m)))
        await pg.goto('file:///home/claude/ed/out.html'); await pg.wait_for_timeout(600)
        await pg.click('#nav button[data-v="energy"]'); await pg.wait_for_timeout(200)
        # 1) Yb/Er preset + suggestions
        await pg.select_option('#edPreset','0'); await pg.wait_for_timeout(300)
        await pg.click('#edSugRun'); await pg.wait_for_timeout(200)
        await pg.screenshot(path=f'{OUT}/10_yber.png',full_page=True)
        # 2) 편집: 에너지 전달 모드로 Er CR 4S3/2->4I9/2 + 4I15/2->4I13/2 추가 (Er 열 복제 후)
        await pg.fill('#edF','NaYF4:2%Er'); await pg.click('#edApply'); await pg.wait_for_timeout(200)
        await pg.click('.ed-ion[data-i="0"] button[data-a="d"]'); await pg.wait_for_timeout(200)
        await pg.click('#edMode button[data-m="et"]')
        cols=await pg.evaluate("ED_UI._state().cols.map(c=>c.uid)")
        async def pick(c,id):
            await pg.click(f'.lvhit[data-c="{c}"][data-id="{id}"]', force=True); await pg.wait_for_timeout(80)
        await pick(cols[0],'4S3/2'); await pick(cols[0],'4I9/2'); await pick(cols[1],'4I15/2'); 
        await pg.screenshot(path=f'{OUT}/11_pick_step4.png')
        await pick(cols[1],'4I13/2')
        m=await pg.evaluate("(()=>{const M=ED_UI._model();return M.procs.filter(p=>p.c&&p.kind==='ET').map(p=>({t:p.type,dE:Math.round(p.c.dE),g:p.c.gradeName,d:Math.round(p.c.delta)}))})()")
        print('CR test',m)
        await pg.click('#edMode button[data-m="view"]'); await pg.wait_for_timeout(100)
        await pg.hover('.proc[data-u]', force=True); await pg.wait_for_timeout(150)
        await pg.click('#edMode button[data-m="et"]')
        await pg.screenshot(path=f'{OUT}/12_er_cr_tip.png')
        # 3) 광흡수 모드: 에러 없이 추가
        await pg.click('#edMode button[data-m="ph"]')
        await pick(cols[0],'4I15/2'); await pick(cols[0],'4I11/2')
        # 4) 방출 모드
        await pg.click('#edMode button[data-m="em"]')
        await pick(cols[0],'4F7/2'); await pick(cols[0],'2H11/2(1)')
        m=await pg.evaluate("ED_UI._model().procs.map(p=>p.type+':'+(p.c?p.c.gradeName:''))")
        print('procs',m)
        await pg.click('#edMode button[data-m="view"]')
        await pg.screenshot(path=f'{OUT}/13_er_edit.png',full_page=True)
        # 5) 다양한 조성
        for f,name in [('K2SiF6:Mn4+','mn4'),('Al2O3:Cr','cr'),('CaF2:Eu','eu2'),('SrF2:Sm2+','sm2'),('Y2O3:5%Eu','eu3'),('LaF3:Ce','ce'),('NaGdF4:49%Yb,1%Tm@NaGdF4:15%Tb','gdtb'),('ZnF2:Mn','mn2'),('NaYF4:Pr,Nd,Sm,Dy,Ho,Tb','many')]:
            await pg.fill('#edF',f); await pg.click('#edApply'); await pg.wait_for_timeout(250)
            await pg.screenshot(path=f'{OUT}/2_{name}.png',full_page=True)
        # 6) 다크 · 모바일
        await pg.select_option('#edPreset','2'); await pg.wait_for_timeout(200)
        await pg.click('#chipTheme'); await pg.wait_for_timeout(300)
        await pg.screenshot(path=f'{OUT}/30_pa_dark.png',full_page=True)
        await pg.set_viewport_size({'width':390,'height':844}); await pg.wait_for_timeout(400)
        await pg.screenshot(path=f'{OUT}/31_mobile.png',full_page=True)
        ov=await pg.evaluate("document.documentElement.scrollWidth>window.innerWidth+1")
        print('mobile overflow',ov)
        # 7) 내보내기 문자열 생성
        await pg.set_viewport_size({'width':1440,'height':950}); await pg.wait_for_timeout(300)

        print('\n'.join(errs) or 'no errors')
        await b.close()
asyncio.run(main())
