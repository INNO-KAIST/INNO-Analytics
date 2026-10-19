"""Energy Diagram·DRE 결함 수정 UI 시험 (10-19 감사 항목 1–10). 사용: python3 ui_fix.py <빌드한 html> [스크린샷 폴더]"""
import asyncio, sys, os, time, json
from playwright.async_api import async_playwright
HTML=sys.argv[1] if len(sys.argv)>1 else '/tmp/edfix/app.html'
OUT=sys.argv[2] if len(sys.argv)>2 else '/tmp/edfix/shots'; os.makedirs(OUT,exist_ok=True)
F=[]
def ok(c,m):
    print(('ok   ' if c else 'FAIL ')+m)
    if not c: F.append(m)
async def toasts(pg): return await pg.evaluate("[...document.querySelectorAll('#toasts .toast')].map(t=>t.textContent).join(' | ')")
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); ctx=await b.new_context(viewport={'width':1440,'height':950},accept_downloads=True); pg=await ctx.new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e))); pg.on('dialog',lambda d: asyncio.ensure_future(d.accept()))
        await pg.goto('file://'+HTML); await pg.wait_for_timeout(500); await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(500)
        # ---- DRE ----
        await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(400)
        msg=await pg.inner_text('#dreSimMsg'); print(' ',msg)
        ok('32.' in msg or '31.' in msg, '1 Lee default s_max ≈ 32'); ok('6.87' in msg,'1 Lee threshold 6.87 kW/cm²')
        ok(await pg.query_selector('#dreConvBadge') is None,'1 no convergence badge at defaults')
        t=time.time(); await pg.fill('input[data-n="sG"][data-f="v"]','1e-10'); await pg.press('input[data-n="sG"][data-f="v"]','Tab'); await pg.wait_for_timeout(50)
        await pg.evaluate("1"); dt=time.time()-t
        ok(dt<5,f'1 sG=1e-10 recompute {dt:.2f} s (< 5 s)')
        bd=await pg.query_selector('#dreConvBadge'); ok(bd is not None,'1 convergence badge shown: '+(await bd.inner_text() if bd else ''))
        ok(await pg.query_selector('#dreFarBadge') is not None,'1 far-outside-range badge')
        await pg.fill('input[data-n="sG"][data-f="v"]','6e-25'); await pg.press('input[data-n="sG"][data-f="v"]','Tab'); await pg.wait_for_timeout(300)
        ok(await pg.query_selector('#dreConvBadge') is None,'1 badge cleared after fixing value')
        # λ
        for bad in ['-5','0']:
            await pg.fill('#dreNm',bad); await pg.press('#dreNm','Tab'); await pg.wait_for_timeout(150)
            ok(await pg.input_value('#dreNm')=='1064' and await pg.evaluate("DRE_UI._state().nm")==1064,f'3/4 λ {bad} rejected and field restored')
        await pg.fill('#dreNm','5000'); await pg.press('#dreNm','Tab'); await pg.wait_for_timeout(300)
        ok(await pg.evaluate("DRE_UI._state().nm")==3000 and await pg.input_value('#dreNm')=='3000','3 λ 5000 clamped to 3000')
        nn=await pg.query_selector('#dreNmNote'); ok(nn is not None and '1064' in (await nn.inner_text()),'3 cross-section wavelength note')
        await pg.fill('#dreNm','1064'); await pg.press('#dreNm','Tab'); await pg.wait_for_timeout(300)
        ok(await pg.query_selector('#dreNmNote') is None,'3 note gone at template λ')
        # Range
        r0=await pg.input_value('#dreR0')
        for bad in ['-1','500','']:
            await pg.fill('#dreR0',bad); await pg.press('#dreR0','Tab'); await pg.wait_for_timeout(120)
            ok(await pg.input_value('#dreR0')==r0,f'4 range min "{bad}" rejected and restored')
        ok('Range' in await toasts(pg),'4 range toast')
        # 7 parameter table range not clipped
        for w in [768,1024,1280,1440,1920]:
            await pg.set_viewport_size({'width':w,'height':950}); await pg.wait_for_timeout(300)
            bad=await pg.evaluate("""(()=>{ const tw=document.querySelector('#drePanel .dre-tw'); const R=tw.getBoundingClientRect(); const out=[];
              tw.querySelectorAll('input[type=text]').forEach(i=>{ const r=i.getBoundingClientRect(); if(r.right>R.right+1||i.scrollWidth>i.clientWidth+1) out.push(i.getAttribute('aria-label')); });
              return {out, hs: tw.scrollWidth>tw.clientWidth+1}; })()""")
            ok(not bad['out'] and not bad['hs'],f'7 param table at {w}px: clipped {bad}')
        await pg.set_viewport_size({'width':1440,'height':950}); await pg.wait_for_timeout(300)
        # 2 import after reload without opening Energy tab
        await pg.click('#nav button[data-v="energy"]'); await pg.wait_for_timeout(200)
        await pg.fill('#edF','β-NaYF4:Yb,Er'); await pg.click('#edApply'); await pg.wait_for_timeout(300)
        await pg.select_option('#edPreset','0'); await pg.wait_for_timeout(300)
        await pg.fill('#edF','β-NaYF4:Yb,Er'); await pg.click('#edApply'); await pg.wait_for_timeout(300)
        np=await pg.evaluate("ED_UI._state().procs.length"); print('  ED procs',np)
        await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(300)
        await pg.reload(); await pg.wait_for_timeout(500)
        if await pg.evaluate("!document.getElementById('dreTpl')"): await pg.click('#nav button[data-v="dre"]'); await pg.wait_for_timeout(300)
        ok(await pg.evaluate("ED_UI._model()===null"),'2 Energy tab not mounted since reload (bug path)')
        print('  ED model before import (null = Energy tab not mounted since reload):', await pg.evaluate("ED_UI._model()===null"))
        await pg.select_option('#dreTpl','ed'); await pg.wait_for_timeout(600)
        st=await pg.evaluate("(()=>{const S=DRE_UI._state(); return {t:S.template,n:S.spec.procs.length,vars:S.samples[0].vars}})()")
        ok(st['t']=='ed' and st['n']>0,f'2 import after reload works ({st})')
        cd=await pg.query_selector('#dreCDef'); ok(cd is not None,'5 "c defaulted to 1 %" note: '+(await cd.inner_text() if cd else ''))
        ok('defaulted' in await toasts(pg),'5 c default toast')
        # 6/10 mobile
        await pg.set_viewport_size({'width':390,'height':844}); await pg.wait_for_timeout(500)
        g=await pg.evaluate("""(()=>{ const s=document.querySelector('#dreChart svg'); const vb=s.viewBox.baseVal.width; const w=s.getBoundingClientRect().width;
          const t=document.getElementById('dreTabs'); const tr=t.getBoundingClientRect();
          return {scale:w/vb, vb, w, tabsOver: tr.right>window.innerWidth+1 || t.scrollWidth>t.clientWidth+1, hs: document.documentElement.scrollWidth>window.innerWidth+1}; })()""")
        ok(g['scale']>0.88,f'6 chart text scale at 390 px {g["scale"]:.2f} (viewBox {g["vb"]}, shown {g["w"]:.0f})')
        ok(not g['tabsOver'],'10 DRE sub-tab bar fits at 390 px'); ok(not g['hs'],'mobile no horizontal scroll')
        await pg.screenshot(path=f'{OUT}/fix_dre_mobile.png',full_page=True)
        await pg.set_viewport_size({'width':1440,'height':950}); await pg.wait_for_timeout(400)
        g2=await pg.evaluate("(()=>{ const s=document.querySelector('#dreChart svg'); return s.getBoundingClientRect().width/s.viewBox.baseVal.width; })()")
        ok(0.9<g2<1.1,f'6 chart redrawn on resize (scale {g2:.2f})')
        # ---- ED ----
        await pg.click('#nav button[data-v="energy"]'); await pg.wait_for_timeout(300)
        await pg.fill('#edF','NaYF4:20%Yb,2%Er,1%Tm,1%Ho,1%Nd'); await pg.click('#edApply'); await pg.wait_for_timeout(800)
        i8=await pg.evaluate("(()=>{ const c=document.getElementById('edIntCard'); return {n:c.querySelectorAll('.ed-int-dir').length, sel:!!c.querySelector('#edIntPair'), opts:c.querySelectorAll('#edIntPair option').length, h:c.getBoundingClientRect().height}; })()")
        ok(i8['sel'] and i8['n']<=2 and i8['h']<6000,f'8 inter-ion card collapsed {i8}')
        o2=await pg.evaluate("document.querySelectorAll('#edIntPair option')[1].value"); await pg.select_option('#edIntPair',o2); await pg.wait_for_timeout(300)
        ok(await pg.evaluate("document.getElementById('edIntPair').value")==o2,'8 pair selector switches')
        w5=await pg.inner_text('#edWarn'); await pg.fill('#edF','NaYF4:Yb1.5,Er0.02'); await pg.click('#edApply'); await pg.wait_for_timeout(300)
        ok('혼용' in await pg.inner_text('#edWarn'),'5 mixed fraction/percent warning shown')
        # 9 CSV columns
        await pg.select_option('#edPreset','0'); await pg.wait_for_timeout(300)
        async with pg.expect_download() as dl: await pg.select_option('#edExport','csv')
        txt=open(await (await dl.value).path(),encoding='utf-8-sig').read().split('\n')
        hdr=txt[1].split(','); ok(all(c in hdr for c in ['laser_cm1','transition_cm1','gap_cm1']),'9 CSV new columns')
        rows=[l.split(',') for l in txt[2:] if l]; ph=[r for r in rows if r[1]=='GSA']; mp=[r for r in rows if r[1]=='MPR']
        ok(ph and ph[0][hdr.index('laser_cm1')]!='' and ph[0][hdr.index('donor_E_cm1')]=='','9 absorption row uses laser_cm1')
        ok(mp and mp[0][hdr.index('gap_cm1')]!='' and mp[0][hdr.index('donor_E_cm1')]=='','9 MPR row uses gap_cm1')
        await pg.fill('#edF',''); await pg.click('#edApply'); await pg.wait_for_timeout(300)
        await pg.evaluate("document.getElementById('toasts').innerHTML=''")
        for k in ['svg','png','csv','json']:
            await pg.select_option('#edExport',k); await pg.wait_for_timeout(100)
        ok((await toasts(pg)).count('Nothing to export')>=4,'9 empty exports toast')
        print('  errors',errs or 'none'); ok(not errs,'no page errors')
        await ctx.close()
        # ---- 10 touch ----
        ctx=await b.new_context(viewport={'width':390,'height':844},has_touch=True,is_mobile=True,device_scale_factor=2); pg=await ctx.new_page()
        pg.on('pageerror',lambda e:errs.append(str(e)))
        await pg.goto('file://'+HTML); await pg.wait_for_timeout(500); await pg.evaluate("localStorage.clear()"); await pg.reload(); await pg.wait_for_timeout(500)
        print('  coarse', await pg.evaluate("matchMedia('(pointer: coarse)').matches"))
        await pg.evaluate("document.querySelector('#nav button[data-v=\"energy\"]').click()"); await pg.wait_for_timeout(300)
        await pg.evaluate("ED_UI._loadPreset(2)"); await pg.wait_for_timeout(300)
        hh=await pg.evaluate("[...document.querySelectorAll('.lvhit')].map(r=>+r.getAttribute('height'))")
        ok(max(hh)>=20,f'10 touch hit areas enlarged (max {max(hh)}, min {min(hh)})')
        await pg.evaluate("document.querySelector('#edMode button[data-m=\"ph\"]').click()"); await pg.wait_for_timeout(200)
        pt=await pg.inner_text('#edPrompt'); ok('Esc' not in pt,'10 no "Esc" hint on touch: '+pt)
        await pg.evaluate("document.querySelector('#edMode button[data-m=\"view\"]').click()")
        await ctx.close(); await b.close()
        print('RESULT', 'all pass' if not F else f'{len(F)} FAIL: {F}')
        if F or errs: sys.exit(1)
asyncio.run(main())
