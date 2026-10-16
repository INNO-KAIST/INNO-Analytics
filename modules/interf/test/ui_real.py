import asyncio, json, sys
from playwright.async_api import async_playwright
F='/home/claude/modules_if/real/ex1.emd'
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await (await b.new_context(viewport={'width':1440,'height':1100})).new_page()
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e)))
        await pg.goto('file:///home/claude/dre/out.html'); await pg.wait_for_timeout(800)
        await pg.set_input_files('#fin',F); await pg.wait_for_timeout(4000)
        await pg.click('#nav button[data-v="spectrum"]'); await pg.wait_for_timeout(2500)
        info=await pg.evaluate("(()=>{const S=loadSpectrum(); const ai=autoIdentify(S,{}); return {spec:state.curSpec, nspec:curRec().spectra.length, det:curRec().spectra.map(s=>s.det), total:S.y.reduce((a,b)=>a+b,0), fw:state.spec.fwhmMn, els:state.spec.els, auto:ai.map(x=>[x.el,x.fam,+(x.xMdl||0).toFixed(1),x.tier||''])}})()")
        print('SPEC',json.dumps(info,ensure_ascii=False))
        await pg.screenshot(path='/home/claude/modules_if/real/spec.png')
        els=sys.argv[1] if len(sys.argv)>1 else ", ".join(e[0] for e in info['auto'])
        await pg.click('#nav button[data-v="interf"]'); await pg.wait_for_timeout(400)
        await pg.fill('#ifEls',els); await pg.click('#ifRun'); await pg.wait_for_timeout(2500)
        r=await pg.evaluate("""(()=>{const R=IF_SUG._last(); return {mode:R.mode, fw:R.fw, fit:R.fit.cols.map(c=>[c.el+' '+c.fam, Math.round(c.A)]),
          out:R.out.map(o=>({el:o.el, best:o.best?[o.best.line.nm,+o.best.line.e.toFixed(3),+o.best.err.toFixed(3),+o.best.purity.toFixed(3),Math.round(o.best.S),o.best.top&&o.best.top.line?o.best.top.line.el+' '+o.best.top.line.nm+' '+(100*o.best.top.frac).toFixed(0)+'%':null, !!o.best.unexpl]:null,
          alt:o.alt?[o.alt.line.nm,+o.alt.line.e.toFixed(3),+o.alt.err.toFixed(3),+o.alt.purity.toFixed(3)]:null}))}})()""")
        print('IF',json.dumps(r,ensure_ascii=False))
        el=await pg.query_selector('#ifSug'); await el.screenshot(path='/home/claude/modules_if/real/sug.png')
        q=await pg.evaluate("""(()=>{const R=IF_SUG._last(), S=R.S; const bands=[[0.2,0.45],[0.45,0.8],[0.8,1.2],[1.2,1.9],[1.9,3],[6,7],[7,8.2],[8.2,9.2],[9.2,11]];
          return bands.map(([a,b])=>{ let y=0,m=0,bg=0,sn=0,c2=0,n=0; for(let i=0;i<S.nch;i++) if(S.E[i]>=a&&S.E[i]<b){ y+=S.y[i]; m+=R.fit.model[i]; bg+=R.fit.bg[i]; sn+=R.fit.snip[i]; c2+=Math.pow(S.y[i]-R.fit.model[i]-R.fit.bg[i],2)/Math.max(1,S.y[i]); n++; } return [a,b,Math.round(y),Math.round(m),Math.round(bg),Math.round(sn),+(c2/n).toFixed(2)]; }) })()""")
        print('BANDS [a,b,data,model,bg(snip+hat),snip,chi2/ch]'); [print(x) for x in q]
        print('errors',errs[:5]); await b.close()
asyncio.run(main())
