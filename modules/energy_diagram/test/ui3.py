import asyncio, json, os
from playwright.async_api import async_playwright
OUT='/home/claude/ed/shots'
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); ctx=await b.new_context(viewport={'width':1440,'height':950},accept_downloads=True); pg=await ctx.new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e)))
        await pg.goto('file:///home/claude/ed/out.html'); await pg.wait_for_timeout(500)
        await pg.click('#nav button[data-v="energy"]')
        await pg.select_option('#edPreset','1'); await pg.wait_for_timeout(200)
        for kind in ['svg','png','csv','json']:
            async with pg.expect_download() as dl:
                await pg.select_option('#edExport',kind)
            d=await dl.value; path=f'{OUT}/export.{kind}'; await d.save_as(path); print(kind, d.suggested_filename, os.path.getsize(path))
        # JSON 불러오기 왕복
        st1=await pg.evaluate("JSON.stringify(ED_UI._state().procs.length)")
        await pg.fill('#edF','Y2O3:Eu'); await pg.click('#edApply'); await pg.wait_for_timeout(100)
        await pg.set_input_files('#edImport',f'{OUT}/export.json'); await pg.wait_for_timeout(300)
        st2=await pg.evaluate("JSON.stringify(ED_UI._state().procs.length)")
        print('roundtrip procs',st1,st2, await pg.evaluate("ED_UI._state().formula"))
        # 키보드: 첫 준위에 초점 → Enter → 상세
        await pg.focus('.lvhit[tabindex="0"]'); await pg.keyboard.press('Enter'); await pg.wait_for_timeout(150)
        print('detail hidden?', await pg.evaluate("document.getElementById('edDetCard').hidden"))
        # d 이온 매개변수 변경
        await pg.fill('#edF','Al2O3:Cr'); await pg.click('#edApply'); await pg.wait_for_timeout(150)
        await pg.click('.lvhit[data-id="2E"]',force=True); await pg.wait_for_timeout(150)
        await pg.fill('#edDetCard input[data-p="Dq"]','1500'); await pg.press('#edDetCard input[data-p="Dq"]','Tab'); await pg.wait_for_timeout(200)
        lv=await pg.evaluate("ED_UI._model().cols[0].lv.map(l=>l.id+':'+Math.round(l.E)+':'+l.regime).join(' ')")
        print('Cr Dq1500',lv)
        em=await pg.evaluate("ED_UI._model().emis.map(e=>e.up.id+'>'+e.lo.id+' '+Math.round(e.nm))")
        print('Cr emission',em)
        # Eu2+ 밴드 입력
        await pg.fill('#edF','SrF2:Eu2+'); await pg.click('#edApply'); await pg.wait_for_timeout(150)
        await pg.click('.lvhit[data-id="5d"]',force=True); await pg.wait_for_timeout(150)
        await pg.fill('#edDetCard input[data-b="E"]','24000'); await pg.press('#edDetCard input[data-b="E"]','Tab'); await pg.wait_for_timeout(200)
        print('Eu2+',await pg.evaluate("ED_UI._model().emis.map(e=>e.up.id+'>'+e.lo.id+' '+Math.round(e.nm))"))
        # 경로 제안 → 추가 → 삭제 (Delete 키)
        await pg.select_option('#edPreset','2'); await pg.wait_for_timeout(150)
        await pg.click('#edSugRun'); await pg.wait_for_timeout(150)
        n0=await pg.evaluate("ED_UI._state().procs.length")
        await pg.click('#edSugList button[data-s]'); await pg.wait_for_timeout(150)
        n1=await pg.evaluate("ED_UI._state().procs.length")
        await pg.focus('.proc.sel'); await pg.keyboard.press('Delete'); await pg.wait_for_timeout(150)
        n2=await pg.evaluate("ED_UI._state().procs.length")
        print('sug add/del',n0,n1,n2)
        # PA 루프 추가 버튼
        await pg.evaluate("(()=>{const S=ED_UI._state(); S.procs.length=0;})()"); await pg.select_option('#edPreset','2'); 
        await pg.evaluate("ED_UI._state().procs.splice(0)"); await pg.fill('#edF','NaYF4:8%Tm'); await pg.click('#edApply'); await pg.wait_for_timeout(150)
        await pg.click('button[data-pa]'); await pg.wait_for_timeout(150)
        print('PA add',await pg.evaluate("ED_UI._model().procs.map(p=>p.type+':'+(p.F?p.F.id+'>'+p.T.id:p.Df.id+'>'+p.Dt.id+'+'+p.Af.id+'>'+p.At.id))"))
        # 새로고침 후 상태 유지
        await pg.reload(); await pg.wait_for_timeout(500); await pg.click('#nav button[data-v="energy"]'); await pg.wait_for_timeout(150)
        print('persist',await pg.evaluate("ED_UI._state().formula+' '+ED_UI._state().procs.length"))
        # 파일 탭으로 돌아가면 파일 열 복귀
        await pg.click('#nav button[data-v="interf"]'); await pg.wait_for_timeout(150)
        print('side visible', await pg.evaluate("getComputedStyle(document.getElementById('side')).display"))
        print('errors',errs or 'none')
        await b.close()
asyncio.run(main())
