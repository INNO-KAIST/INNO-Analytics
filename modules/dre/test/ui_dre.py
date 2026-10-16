"""DRE 탭 UI 시험 (Playwright). 사용: python3 ui_dre.py <빌드한 html> <스크린샷 폴더>
  예: python3 ui_dre.py /tmp/out.html /tmp/shots   (modules/inject_all.py 로 만든 html)
  인자를 생략하면 아래 기본 경로를 쓴다. 화면 글은 영어 골격(탭·버튼·라벨)이므로 선택자는 id / data-* 만 쓴다."""
import asyncio, json, sys
from playwright.async_api import async_playwright
HTML=sys.argv[1] if len(sys.argv)>1 else '/home/claude/dre/out.html'
OUT=sys.argv[2] if len(sys.argv)>2 else '/home/claude/dre/shots'
import os; os.makedirs(OUT,exist_ok=True)
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); ctx=await b.new_context(viewport={'width':1440,'height':950},accept_downloads=True); pg=await ctx.new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append('pageerror '+str(e)))
        pg.on('console',lambda m: errs.append('console '+m.text) if m.type=='error' and 'ERR_TUNNEL' not in m.text else None)
        pg.on('dialog',lambda d: asyncio.ensure_future(d.accept()))
        await pg.goto('file://'+HTML); await pg.wait_for_timeout(500)
        await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(400)
        await pg.screenshot(path=f'{OUT}/d0_lee.png',full_page=True)
        # 합성 데이터
        await pg.click('#dreTabs button[data-t="data"]')
        await pg.click('#dreSynth'); await pg.click('#sOK'); await pg.wait_for_timeout(300)
        # 붙여넣기 데이터(두 번째 시료) 
        await pg.click('#dreAddSmp')
        await pg.fill('#drePaste','P_kW,sig\n1,5\n2,10\n5,200\n10,5000\n20,20000\n')
        await pg.select_option('#dreDSmp','1'); await pg.click('#dreParse'); await pg.click('#mOK'); await pg.wait_for_timeout(200)
        await pg.screenshot(path=f'{OUT}/d1_data.png',full_page=True)
        # 두 번째 시료 제거 후 피팅
        await pg.click('button[data-delsmp="1"]'); await pg.wait_for_timeout(200)
        await pg.click('#dreTabs button[data-t="fit"]')
        await pg.fill('#fS','8'); await pg.fill('#fN','7'); await pg.press('#fN','Tab')
        await pg.click('#fRun')
        for i in range(240):
            await pg.wait_for_timeout(1000)
            done=await pg.evaluate("!!(DRE_UI._state().result)")
            if done: break
        print('fit done after',i,'s')
        await pg.screenshot(path=f'{OUT}/d2_fit.png',full_page=True)
        r=await pg.evaluate("(()=>{const R=DRE_UI._state().result; return {chi2:R.best.chi2,dof:R.info.dof,succ:R.ms.filter(q=>q.chi2-R.best.chi2<=0.5).length+'/'+R.ms.length,prof:R.prof.map((p,j)=>R.free[j].key+':'+(p.hi-p.lo).toFixed(2)+(p.openLo?'↓':'')+(p.openHi?'↑':''))}})()")
        print(json.dumps(r,ensure_ascii=False))
        # ED 가져오기
        await pg.click('#nav button[data-v="energy"]'); await pg.wait_for_timeout(200)
        await pg.select_option('#edPreset','0'); await pg.wait_for_timeout(300)
        await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(200)
        await pg.select_option('#dreTpl','ed'); await pg.wait_for_timeout(600)
        await pg.click('#dreTabs button[data-t="model"]')
        await pg.screenshot(path=f'{OUT}/d3_ed.png',full_page=True)
        m=await pg.evaluate("(()=>{const S=DRE_UI._state(); return {sp:S.spec.species.map(s=>s.id+':'+s.levels.join('/')), procs:S.spec.procs.length, ch:S.spec.channels.map(c=>c.id), err:DRE_UI._model().err, simErr:DRE_UI._sim().err}})()")
        print(json.dumps(m,ensure_ascii=False))
        # 다크·모바일
        await pg.click('#chipTheme'); await pg.wait_for_timeout(300)
        await pg.screenshot(path=f'{OUT}/d4_dark.png',full_page=True)
        await pg.set_viewport_size({'width':390,'height':844}); await pg.wait_for_timeout(400)
        await pg.screenshot(path=f'{OUT}/d5_mobile.png',full_page=True)
        print('mobile overflow',await pg.evaluate("document.documentElement.scrollWidth>window.innerWidth+1"))
        print('errors',errs or 'none')
        await b.close()
asyncio.run(main())
