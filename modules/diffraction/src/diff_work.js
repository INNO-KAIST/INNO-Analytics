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
  const DC=G.DIFF_CORE, DP=G.DIFF_PATTERN, DI=G.DIFF_INDEX; const PH=new Map(); // 기록 id|version → phase
  const phaseOf=r=>{ const k=(r.id||r.name)+"|"+(r.version||0)+"|"+(r.updatedAt||""); let p=PH.get(k); if(!p){ p=DC.phaseFromRecord(r); p.id=r.id; p.layer=r.layer; PH.set(k,p); } return p; };
  const slimFit=L=>L&&{center:L.center,R1:L.R1,R2:L.R2,R3:L.R3,len:L.len,sig:L.sig,ang12:L.ang12,ang13:L.ang13,ang23:L.ang23,sigAng:L.sigAng,sigCenter:L.sigCenter,radial:L.radial,radialPct:L.radialPct,radialUsed:L.radialUsed,
    spots:L.spots.map(s=>({x:s.x,y:s.y,I:s.I,snr:s.snr,saturated:!!s.saturated,radius:s.radius||0,n:s.n,res:s.res,indexed:s.indexed})),nIndexed:L.nIndexed,fracI:L.fracI,occupancy:L.occupancy,sigmaPx:L.sigmaPx,tol:L.tol};
  const slimCand=c=>({phaseId:c.phase.id,phaseName:c.phaseName,formula:c.phase.formula,sg:c.phase.sg,layer:c.phase.layer,zone:c.zone,h1:c.h1,h2:c.h2,h3:c.h3,g:c.g,d:c.d,ang12:c.ang12,ratio21:c.ratio21,fom:c.fom,kappa:c.kappa,scaleImplied:c.scaleImplied,
    forb:c.forb,rho:c.rho,missFrac:c.missFrac,score:c.score,hasF:c.hasF,swap:c.swap,
    verify:c.verify&&{nObs:c.verify.nObs,nStrongPred:c.verify.nStrongPred,completeness:c.verify.completeness,doubleDiff:c.verify.doubleDiff,missing:c.verify.missing.slice(0,12),rows:c.verify.rows.map(r=>({x:r.spot.x,y:r.spot.y,n:r.spot.n,h:r.h,d:r.d,I:r.I,absent:r.absent,weak:r.weak}))}});
  return function(cmd,a){
    if(cmd==="analyze"){ const {w,h}=a, o=a.opt||{}; const t0=Date.now(); const T={}; const SZ=DP.sanitize(a.img); const img=SZ.img;
      const ov=DP.findOverlay(img,w,h); const bs=DP.findBeamStop(img,w,h); const mask=new Uint8Array(w*h); for(let i=0;i<mask.length;i++) mask[i]=bs.mask[i]|ov.mask[i]|(SZ.bad?SZ.bad[i]:0); T.mask=Date.now()-t0;
      const out={w,h,overlay:{bars:ov.bars,text:ov.text},scaleBar:DP.scaleBar(ov,w,h),beamStop:bs.parts,clipped:!!bs.clipped,nonFinite:SZ.n,mask};
      const doSpots=o.pattern!=="rings", doRings=o.pattern!=="spots";
      if(doSpots){ const t=Date.now(); const S=DP.findSpots(img,w,h,{mask}); let c=o.center?Object.assign({},o.center,{source:"user",origin:"user"}):DP.findCenter(S.spots,w,h,{beamStop:bs.mask});
        const subset=o.topN>0?S.spots.slice(0,o.topN):S.spots; const M=c?DP.fitLattices(subset,c,{size:Math.max(w,h),radial:o.radial==null?"auto":o.radial,max:o.maxLattices||1}):{lattices:[],unindexed:subset};
        out.spots={list:S.spots.map(s=>({x:s.x,y:s.y,I:s.I,snr:s.snr,saturated:!!s.saturated,radius:s.radius||0})),noise:S.noise,saturated:S.saturated,center:c,lattices:M.lattices.map(slimFit),nUnindexed:M.unindexed.length}; T.spots=Date.now()-t; }
      // 자동: 스폿 격자가 뚜렷하면(지수 12 개 이상, 설명 85 % 이상, 기저가 영상의 3 % 이상 — 링의 얼룩 격자가 아님) 링 분석을 건너뛴다 (휴대폰 시간 절약)
      let skipRings=false; if(doSpots&&doRings&&o.pattern==="auto"){ const L=out.spots.lattices[0]; skipRings=!!(L&&L.nIndexed>=12&&L.fracI>=0.85&&Math.min(L.len[0],L.len[1])>=0.03*Math.max(w,h)); }
      if(doRings&&!skipRings){ const t=Date.now(); const R=DP.analyzeRings(img,w,h,{mask,center:o.ringCenter||null}); const pr=R.profile;
        out.rings={ellipse:R.ellipse,ratio:R.ratio,center:R.center,coarse:R.coarse,unc:R.unc,rStop:R.rStop,peaks:R.peaks,ringFits:R.ringFits.map(f=>({r:f.r,n:f.n,nSig:f.nSig,medSig:f.medSig,cx:f.cx,cy:f.cy,a:f.a,b:f.b,phi:f.phi})),profile:{r:pr.r,I:pr.I,n:pr.n},background:Array.from(R.background)}; T.rings=Date.now()-t; }
      // 자동 판정: 링 맞춤이 둘 이상이고 광선 점이 넉넉(≥ 60 %)하며 서로 중심이 맞으면 링, 아니면 스폿
      if(skipRings){ out.type="spots"; out.typeWhy={ringsSkipped:true,spotsIndexed:out.spots.lattices[0].nIndexed}; }
      else if(doSpots&&doRings){ const F=out.rings.ringFits; const good=F.filter(f=>f.n>=0.6*360&&(f.medSig==null||f.medSig>=10)); /* 광선마다 뚜렷한 봉우리(돌출 ÷ 광선 잡음 중앙값 ≥ 10)여야 링 — 둘레 번짐의 매끈한 꼭대기는 4–7 (실측 Fe-Cr), Au 링은 17–68 */ let rings=false; if(good.length>=2){ const cx=good.map(f=>f.cx), cy=good.map(f=>f.cy); const sp=Math.hypot(Math.max(...cx)-Math.min(...cx),Math.max(...cy)-Math.min(...cy)); rings=sp<3; }
        const L=out.spots.lattices[0]; const spotsOk=L&&L.nIndexed>=6; out.type=rings&&!(spotsOk&&L.fracI>0.8&&good.length<3)?"rings":(spotsOk?"spots":(rings?"rings":"spots")); out.typeWhy={ringFitsGood:good.length,spotsIndexed:L?L.nIndexed:0}; }
      else out.type=doRings?"rings":"spots";
      out.ms=Date.now()-t0; out.timing=T; return out; }
    if(cmd==="refit"){ if(a.hint){ const c=DP.findCenter(a.spots,a.w||a.size,a.h||a.size,{hint:a.center,beamStop:null}); if(c) a.center={x:c.x,y:c.y,source:"user",origin:c.origin}; } const M=DP.fitLattices(a.spots,a.center,{size:a.size,radial:a.radial==null?"auto":a.radial,max:a.max||1,weightPow:a.weightPow}); return {lattices:M.lattices.map(slimFit),nUnindexed:M.unindexed.length,center:a.center}; }
    if(cmd==="rankSpots"){ const phases=[]; const bad=[]; for(const r of a.records){ try{ phases.push(phaseOf(r)); }catch(e){ bad.push({id:r.id,error:String(e.message||e)}); } }
      const R=DI.rankSpots(a.fit,phases,a.opt||{}); return {meas:R.meas,cands:R.cands.map(slimCand),bad}; }
    if(cmd==="matchRings"){ const phases=[]; const bad=[]; for(const r of a.records){ try{ phases.push(phaseOf(r)); }catch(e){ bad.push({id:r.id,error:String(e.message||e)}); } }
      const M=DI.matchRings(a.peaks,phases,a.opt||{}); return {cands:M.map(m=>({phaseId:m.phase.id,phaseName:m.phaseName,formula:m.phase.formula,sg:m.phase.sg,layer:m.phase.layer,scale:m.scale,scaleErr:m.scaleErr,kappa:m.kappa,explained:m.explained,missing:m.missing,dev:m.dev,rho:m.rho,score:m.score,pairs:m.pairs,
        rings:m.rings.filter(r=>r.g<=a.gShow).map(r=>({g:r.g,d:r.d,m:r.m,hkl:r.hkls[0],I:r.I}))})),bad}; }
    if(cmd==="zoneOf"){ const p=phaseOf(a.record); return DC.zoneBasis(p,a.zone,{gmax:a.gmax||15}); }
    throw new Error("알 수 없는 명령 "+cmd); };
}
const DW={}; let worker=null, seq=0; const pend=new Map(); let local=null; DW.mode="none";
function srcText(){ const ids=["diff-src-sg","diff-src-scat","diff-src-core","diff-src-pattern","diff-src-index"]; let s=""; for(const id of ids){ const el=document.getElementById(id); if(!el) return null; s+=el.textContent+"\n"; } return s; }
function ensure(){ if(worker||local) return; try{ const src=srcText(); if(!src||typeof Worker==="undefined") throw new Error("no src");
    const code=src+"\nconst H=("+handler.toString()+")(self);\nself.onmessage=e=>{ const {id,cmd,a}=e.data; try{ const r=H(cmd,a); const tr=[]; if(r&&r.mask&&r.mask.buffer) tr.push(r.mask.buffer); self.postMessage({id,ok:true,r},tr); }catch(err){ self.postMessage({id,ok:false,e:String(err&&err.message||err)}); } };";
    const url=URL.createObjectURL(new Blob([code],{type:"text/javascript"})); worker=new Worker(url); DW.mode="worker";
    worker.onmessage=e=>{ const p=pend.get(e.data.id); if(!p) return; pend.delete(e.data.id); if(e.data.ok) p.res(e.data.r); else p.rej(new Error(e.data.e)); };
    worker.onerror=e=>{ for(const p of pend.values()) p.rej(new Error("Worker 오류: "+(e.message||""))); pend.clear(); worker=null; local=handler(G); DW.mode="main"; };
  }catch(e){ local=handler(G); DW.mode="main"; } }
/* run(cmd, args, transfer) → Promise. analyze 의 img 는 복사본을 넘긴다(전송 후 원본이 비지 않게 호출 쪽이 결정) */
DW.run=function(cmd,a,transfer){ ensure(); if(local) return new Promise((res,rej)=>setTimeout(()=>{ try{ res(local(cmd,a)); }catch(e){ rej(e); } },0));
  return new Promise((res,rej)=>{ const id=++seq; pend.set(id,{res,rej}); worker.postMessage({id,cmd,a},transfer||[]); }); };
DW.handler=handler;
if(typeof module!=="undefined"&&module.exports) module.exports=DW;
G.DIFF_WORK=DW;
})();
