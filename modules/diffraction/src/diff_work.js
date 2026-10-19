/* ============================================================
   INNO Analytics — Diffraction: 무거운 계산을 Web Worker 로 (2048² 에서 수 초).
   Worker 원본 = 주입된 <script id="diff-src-…"> 들의 글자 + 아래 처리기. Worker 를 못 만들면 같은 함수를 주 스레드에서 (알림).
   명령: analyze · refit · rankSpots · matchRings · zoneOf
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
/* ---- 처리기 (Worker 안·주 스레드 공용) ---- */
function handler(G){
  const DC=G.DIFF_CORE, DP=G.DIFF_PATTERN, DI=G.DIFF_INDEX; const PH=new Map(); let PROG=null; // 기록 id|version → phase · PROG = 진행 알림 (Worker 면 postMessage)
  const phaseOf=r=>{ const k=(r.id||r.name)+"|"+(r.version||0)+"|"+(r.updatedAt||""); let p=PH.get(k); if(!p){ p=DC.phaseFromRecord(r); p.id=r.id; p.layer=r.layer; PH.set(k,p); } return p; };
  const slimFit=L=>L&&{center:L.center,R1:L.R1,R2:L.R2,R3:L.R3,len:L.len,sig:L.sig,ang12:L.ang12,ang13:L.ang13,ang23:L.ang23,sigAng:L.sigAng,sigCenter:L.sigCenter,radial:L.radial,radialPct:L.radialPct,radialUsed:L.radialUsed,
    spots:L.spots.map(s=>({x:s.x,y:s.y,I:s.I,snr:s.snr,saturated:!!s.saturated,radius:s.radius||0,n:s.n,res:s.res,indexed:s.indexed})),nIndexed:L.nIndexed,fracI:L.fracI,occupancy:L.occupancy,sigmaPx:L.sigmaPx,tol:L.tol};
  const slimCand=c=>({phaseId:c.phase.id,phaseName:c.phaseName,formula:c.phase.formula,sg:c.phase.sg,layer:c.phase.layer,zone:c.zone,h1:c.h1,h2:c.h2,h3:c.h3,g:c.g,d:c.d,ang12:c.ang12,ratio21:c.ratio21,fom:c.fom,kappa:c.kappa,scaleImplied:c.scaleImplied,
    forb:c.forb,rho:c.rho,missFrac:c.missFrac,score:c.score,hasF:c.hasF,swap:c.swap,
    verify:c.verify&&{nObs:c.verify.nObs,nStrongPred:c.verify.nStrongPred,completeness:c.verify.completeness,doubleDiff:c.verify.doubleDiff,missing:c.verify.missing.slice(0,12),rows:c.verify.rows.map(r=>({x:r.spot.x,y:r.spot.y,n:r.spot.n,h:r.h,d:r.d,I:r.I,absent:r.absent,weak:r.weak,Id:r.Id,dd:r.dd,gm:r.gm,dz:r.dz}))},kin:c.kin||null,dyn:c.dyn||null,dynError:c.dynError||null});
  return function(cmd,a,onProgress){ PROG=onProgress||null;
    if(cmd==="analyze"){ const t0=Date.now(); const {w,h}=a, o=a.opt||{}; const SZ=DP.sanitize(a.img); const imgF=SZ.img;
      // 큰 이미지(긴 변 > 1024): f×f 평균으로 줄여(≤ 1024) 찾고, 위치·중심·격자·링·스케일바는 원래 해상도에서 다듬는다 (2048² 2–3 배, 4096² 약 5 배 빠름, 메모리 ↓)
      const f=o.full?1:Math.max(1,Math.ceil(Math.max(w,h)/(o.work||1024)));
      // 스케일바 이름표 읽기 (DIFF_OCR, 원래 해상도) — 막대를 못 찾았으면 OCR 의 막대 찾기가 받는다
      const label=(out)=>{ const OC=G.DIFF_OCR; if(!OC) return; try{ const t=Date.now(); const r=OC.readAuto(imgF,w,h,out.scaleBar,{expect:"recip"});
          if(r){ out.scaleLabel={text:r.text,value:r.value,unit:r.unit,valueInvNm:r.valueInvNm,recip:r.recip,conf:r.conf,unitConf:r.unitConf!=null?r.unitConf:r.conf,unitScore:r.unitScore,box:r.box};
            if(!out.scaleBar&&r.bar){ const b=r.bar; out.scaleBar={x0:b.x0,x1:b.x1,y0:b.y0,y1:b.y1,lenPx:b.lenPx||(b.x1-b.x0+1),outerPx:b.outerPx||(b.x1-b.x0+1),hollow:!!b.hollow,fromOCR:true}; } }
          if(out.timing) out.timing.ocr=Date.now()-t; }catch(e){ out.ocrError=String(e.message||e); } };
      if(f===1){ const out=core(imgF,w,h,o,SZ); label(out); out.ms=Date.now()-t0; out.work=1; return out; }
      const D=DP.downsample(imgF,w,h,f); const sc=v=>v==null?null:Object.assign({},v,{x:(v.x-(f-1)/2)/f,y:(v.y-(f-1)/2)/f});
      const co=core(D.img,D.w,D.h,Object.assign({},o,{center:sc(o.center),ringCenter:sc(o.ringCenter)}),null); const T=co.timing; const t1=Date.now();
      const mask=DP.upMask(co.mask,D.w,D.h,w,h,f); if(SZ.bad) for(let i=0;i<mask.length;i++) if(SZ.bad[i]) mask[i]=1;
      const up=v=>({x:v.x*f+(f-1)/2,y:v.y*f+(f-1)/2});
      const out={w,h,work:f,overlay:{bars:co.overlay.bars.map(b=>Object.assign({},b,{x0:b.x0*f,x1:b.x1*f+f-1,y0:b.y0*f,y1:b.y1*f+f-1,len:b.len*f})),text:co.overlay.text.map(t=>({x0:t.x0*f,x1:t.x1*f+f-1,y0:t.y0*f,y1:t.y1*f+f-1}))},
        scaleBar:DP.refineBar(imgF,w,h,co.scaleBar,f),beamStop:co.beamStop.map(p=>Object.assign({},p,{n:p.n*f*f,cx:p.cx*f+(f-1)/2,cy:p.cy*f+(f-1)/2,x0:p.x0*f,x1:p.x1*f+f-1,y0:p.y0*f,y1:p.y1*f+f-1})),clipped:co.clipped,nonFinite:SZ.n,mask,type:co.type,typeWhy:co.typeWhy};
      if(co.spots){ const ref=DP.refineSpots(imgF,w,h,co.spots.list,f,{mask}); const c0=co.spots.center; let c=null;
        if(o.center) c=Object.assign({},o.center,{source:"user",origin:"user"});
        else if(c0){ const cu=up(c0); if(c0.origin==="spot"){ let b=null,bd=1e9; for(const s of ref){ const d=Math.hypot(s.x-cu.x,s.y-cu.y); if(d<bd){ bd=d; b=s; } } c=Object.assign({},c0,b&&bd<3*f?{x:b.x,y:b.y}:cu); }
          else { const r=DP.findCenter(ref,w,h,{hint:cu}); c=Object.assign({},c0,r&&r.pairs?{x:r.x,y:r.y,pairs:r.pairs}:cu); } }
        const subset=o.topN>0?ref.slice().sort((p,q)=>q.I-p.I).slice(0,o.topN):ref; const M=c?DP.fitLattices(subset,c,{size:Math.max(w,h),radial:o.radial==null?"auto":o.radial,max:o.maxLattices||1}):{lattices:[],unindexed:subset};
        out.spots={list:ref.map(s=>({x:s.x,y:s.y,I:s.I,snr:s.snr,saturated:!!s.saturated,radius:s.radius||0})),noise:co.spots.noise,saturated:co.spots.saturated,center:c,lattices:M.lattices.map(slimFit),nUnindexed:M.unindexed.length}; }
      if(co.rings){ const e=co.rings.ellipse; const init={cx:e.cx*f+(f-1)/2,cy:e.cy*f+(f-1)/2,a:e.a,b:e.b,phi:e.phi}; const R=DP.analyzeRings(imgF,w,h,{mask,init,iters:1}); const pr=R.profile;
        out.rings={ellipse:R.ellipse,ratio:R.ratio,center:R.center,coarse:{x:co.rings.coarse.x*f+(f-1)/2,y:co.rings.coarse.y*f+(f-1)/2},unc:R.unc,rStop:R.rStop,peaks:R.peaks,ringFits:R.ringFits.map(q=>({r:q.r,n:q.n,nIn:q.nIn,nSig:q.nSig,medSig:q.medSig,cx:q.cx,cy:q.cy,a:q.a,b:q.b,phi:q.phi})),profile:{r:pr.r,I:pr.I,n:pr.n},background:Array.from(R.background)}; }
      if(o.pattern==="auto"&&out.rings&&out.spots&&!(co.typeWhy&&co.typeWhy.ringsSkipped)){ const K=classify(out.rings.ringFits,out.rings.peaks,out.spots.lattices[0],Math.max(w,h)); out.type=K.type; out.typeWhy=K.why; }
      T.refine=Date.now()-t1; out.timing=T; label(out); out.ms=Date.now()-t0; return out; }
    return rest(cmd,a); };
  /* 링 / 스폿 판정
     · 믿을 만한 스폿 격자 = 지수 6 개 이상, 기저 ≥ 이미지의 3 %, 지수 붙은 세기 ≥ 35 % (링·결정립 얼룩을 맞춘 작은 가짜 격자가 아님)
     · 링 맞춤의 광선 점 수는 '이미지 안에 보이는 원호'(nIn) 기준 — 중심이 치우쳐 가장자리에 잘린 링도 링
     · 엄격: 점 ≥ 60 % · 광선 잡음 대비 돌출 중앙값 ≥ 10 (Fe-Cr 실측의 둘레 번짐 4–7, Au 링 17–127)
     · 느슨(믿을 만한 격자가 없을 때만): 점 ≥ 45 % · 돌출 ≥ 6 — 약한 링, 8-bit 로그 표시 사본
     · 링 2 개 이상이 같은 중심(3 px 안)이면 링 */
  function classify(F,peaks,L,size){ F=F||[]; const nR=360; const vis=f=>Math.max(f.nIn!=null?f.nIn:nR,1); const arc=f=>f.nIn==null||f.nIn>=0.25*nR;
    const strict=F.filter(f=>arc(f)&&f.n>=0.6*vis(f)&&(f.medSig==null||f.medSig>=10)), loose=F.filter(f=>arc(f)&&f.n>=0.45*vis(f)&&(f.medSig==null||f.medSig>=6));
    const conc=G=>{ if(G.length<2) return false; const mx=G.map(f=>f.cx).sort((a,b)=>a-b)[G.length>>1], my=G.map(f=>f.cy).sort((a,b)=>a-b)[G.length>>1]; return G.filter(f=>Math.hypot(f.cx-mx,f.cy-my)<3).length>=2; };
    const latOk=!!(L&&L.nIndexed>=6&&Math.min(L.len[0],L.len[1])>=0.03*size&&(L.fracI==null||L.fracI>=0.35)); const speckle=!!(L&&Math.min(L.len[0],L.len[1])<0.03*size);
    const ringsStrict=conc(strict), ringsLoose=!latOk&&conc(loose); const rings=ringsStrict||ringsLoose;
    let type; if(rings&&!(latOk&&L.fracI>0.8&&strict.length<3)) type="rings"; else type=latOk?"spots":(rings?"rings":"spots");
    return {type,why:{ringFitsGood:strict.length,ringFitsLoose:loose.length,spotsIndexed:L?L.nIndexed:0,latticeOk:latOk,speckle}}; }
  function core(img,w,h,o,SZ){ const t0=Date.now(); const T={};
      const ov=DP.findOverlay(img,w,h); const bs=DP.findBeamStop(img,w,h); const mask=new Uint8Array(w*h); for(let i=0;i<mask.length;i++) mask[i]=bs.mask[i]|ov.mask[i]|(SZ&&SZ.bad?SZ.bad[i]:0); T.mask=Date.now()-t0;
      const out={w,h,overlay:{bars:ov.bars,text:ov.text},scaleBar:DP.scaleBar(ov,w,h),beamStop:bs.parts,clipped:!!bs.clipped,nonFinite:SZ?SZ.n:0,mask};
      const doSpots=o.pattern!=="rings", doRings=o.pattern!=="spots";
      if(doSpots){ const t=Date.now(); const S=DP.findSpots(img,w,h,{mask}); let c=o.center?Object.assign({},o.center,{source:"user",origin:"user"}):DP.findCenter(S.spots,w,h,{beamStop:bs.mask});
        const subset=o.topN>0?S.spots.slice(0,o.topN):S.spots; const M=c?DP.fitLattices(subset,c,{size:Math.max(w,h),radial:o.radial==null?"auto":o.radial,max:o.maxLattices||1}):{lattices:[],unindexed:subset};
        out.spots={list:S.spots.map(s=>({x:s.x,y:s.y,I:s.I,snr:s.snr,saturated:!!s.saturated,radius:s.radius||0})),noise:S.noise,saturated:S.saturated,center:c,lattices:M.lattices.map(slimFit),nUnindexed:M.unindexed.length}; T.spots=Date.now()-t; }
      // 자동: 스폿 격자가 뚜렷하면(지수 12 개 이상, 설명 85 % 이상, 기저가 이미지의 3 % 이상 — 링의 얼룩 격자가 아님) 링 분석을 건너뛴다 (휴대폰 시간 절약)
      let skipRings=false; if(doSpots&&doRings&&o.pattern==="auto"){ const L=out.spots.lattices[0]; skipRings=!!(L&&L.nIndexed>=12&&L.fracI>=0.85&&Math.min(L.len[0],L.len[1])>=0.03*Math.max(w,h)); }
      if(doRings&&!skipRings){ const t=Date.now(); const R=DP.analyzeRings(img,w,h,{mask,center:o.ringCenter||null}); const pr=R.profile;
        out.rings={ellipse:R.ellipse,ratio:R.ratio,center:R.center,coarse:R.coarse,unc:R.unc,rStop:R.rStop,peaks:R.peaks,ringFits:R.ringFits.map(f=>({r:f.r,n:f.n,nIn:f.nIn,nSig:f.nSig,medSig:f.medSig,cx:f.cx,cy:f.cy,a:f.a,b:f.b,phi:f.phi})),profile:{r:pr.r,I:pr.I,n:pr.n},background:Array.from(R.background)}; T.rings=Date.now()-t; }
      // 자동 판정 (classify): 링이냐 스폿이냐 — 큰 이미지는 원래 해상도 결과로 다시 판정한다 (아래 analyze)
      if(skipRings){ out.type="spots"; out.typeWhy={ringsSkipped:true,spotsIndexed:out.spots.lattices[0].nIndexed}; }
      else if(doSpots&&doRings){ const K=classify(out.rings.ringFits,out.rings.peaks,out.spots.lattices[0],Math.max(w,h)); out.type=K.type; out.typeWhy=K.why; }
      else out.type=doRings?"rings":"spots";
      out.ms=Date.now()-t0; out.timing=T; return out; }
  function rest(cmd,a){ // analyze 밖의 명령
    if(cmd==="refit"){ if(a.hint){ const c=DP.findCenter(a.spots,a.w||a.size,a.h||a.size,{hint:a.center,beamStop:null}); if(c) a.center={x:c.x,y:c.y,source:"user",origin:c.origin}; } const M=DP.fitLattices(a.spots,a.center,{size:a.size,radial:a.radial==null?"auto":a.radial,max:a.max||1,weightPow:a.weightPow}); return {lattices:M.lattices.map(slimFit),nUnindexed:M.unindexed.length,center:a.center}; }
    if(cmd==="rankSpots"){ const phases=[]; const bad=[]; for(const r of a.records){ try{ phases.push(phaseOf(r)); }catch(e){ bad.push({id:r.id,error:String(e.message||e)}); } }
      const R=DI.rankSpots(a.fit,phases,Object.assign({},a.opt||{},{onProgress:PROG})); return {meas:R.meas,cands:R.cands.map(slimCand),bad,refused:R.refused||null,truncated:!!R.truncated,priorEmpty:!!R.priorEmpty,ms:R.ms}; }
    if(cmd==="matchRings"){ const phases=[]; const bad=[]; for(const r of a.records){ try{ phases.push(phaseOf(r)); }catch(e){ bad.push({id:r.id,error:String(e.message||e)}); } }
      const M=DI.matchRings(a.peaks,phases,a.opt||{}); return {cands:M.map(m=>({phaseId:m.phase.id,phaseName:m.phaseName,formula:m.phase.formula,sg:m.phase.sg,layer:m.phase.layer,scale:m.scale,scaleErr:m.scaleErr,kappa:m.kappa,explained:m.explained,missing:m.missing,dev:m.dev,rho:m.rho,score:m.score,pairs:m.pairs,
        rings:m.rings.filter(r=>r.g<=a.gShow).map(r=>({g:r.g,d:r.d,m:r.m,hkl:r.hkls[0],I:r.I}))})),bad,priorEmpty:!!M.priorEmpty}; }
    if(cmd==="zoneOf"){ const p=phaseOf(a.record); return DC.zoneBasis(p,a.zone,{gmax:a.gmax||15}); }
    throw new Error("알 수 없는 명령 "+cmd); };
}
const DW={}; let worker=null, seq=0; const pend=new Map(); let local=null; DW.mode="none";
function srcText(){ const ids=["diff-src-sg","diff-src-scat","diff-src-core","diff-src-dyn","diff-src-ocr","diff-src-pattern","diff-src-index"]; let s=""; const opt=new Set(["diff-src-dyn","diff-src-ocr"]); /* 선택 모듈은 없어도 됨 */ for(const id of ids){ const el=document.getElementById(id); if(!el){ if(opt.has(id)) continue; return null; } s+=el.textContent+"\n"; } return s; }
function ensure(){ if(worker||local) return; try{ const src=srcText(); if(!src||typeof Worker==="undefined") throw new Error("no src");
    const code=src+"\nconst H=("+handler.toString()+")(self);\nself.onmessage=e=>{ const {id,cmd,a}=e.data; try{ const r=H(cmd,a,p=>self.postMessage({id,progress:p})); const tr=[]; if(r&&r.mask&&r.mask.buffer) tr.push(r.mask.buffer); self.postMessage({id,ok:true,r},tr); }catch(err){ self.postMessage({id,ok:false,e:String(err&&err.message||err)}); } };";
    const url=URL.createObjectURL(new Blob([code],{type:"text/javascript"})); worker=new Worker(url); DW.mode="worker";
    worker.onmessage=e=>{ const p=pend.get(e.data.id); if(!p) return; if(e.data.progress!=null){ if(p.onProgress) try{ p.onProgress(e.data.progress); }catch(err){} return; } pend.delete(e.data.id); if(e.data.ok) p.res(e.data.r); else p.rej(new Error(e.data.e)); };
    worker.onerror=e=>{ for(const p of pend.values()) p.rej(new Error("Worker 오류: "+(e.message||""))); pend.clear(); worker=null; local=handler(G); DW.mode="main"; };
  }catch(e){ local=handler(G); DW.mode="main"; } }
/* run(cmd, args, transfer) → Promise. analyze 의 img 는 복사본을 넘긴다(전송 후 원본이 비지 않게 호출 쪽이 결정) */
/* opt.onProgress(frac) — Worker 의 진행 알림 (rankSpots). DW.cancel(): Worker 를 끝내 계산을 멈추고 기다리던 약속은 "cancelled" 로 거절 (다음 run 에서 새로 만듦) */
DW.run=function(cmd,a,transfer,opt){ opt=opt||{}; ensure(); if(local) return new Promise((res,rej)=>setTimeout(()=>{ try{ res(local(cmd,a,opt.onProgress)); }catch(e){ rej(e); } },0));
  return new Promise((res,rej)=>{ const id=++seq; pend.set(id,{res,rej,onProgress:opt.onProgress}); worker.postMessage({id,cmd,a},transfer||[]); }); };
DW.cancel=function(){ if(!worker) return false; try{ worker.terminate(); }catch(e){} worker=null; for(const p of pend.values()) p.rej(Object.assign(new Error("cancelled"),{cancelled:true})); pend.clear(); return true; };
DW.handler=handler;
if(typeof module!=="undefined"&&module.exports) module.exports=DW;
G.DIFF_WORK=DW;
})();
