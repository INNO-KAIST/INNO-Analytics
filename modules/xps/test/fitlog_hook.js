/* 피팅 기록 고리 (SPEC §3.3b 검증용): node -r ./test/fitlog_hook.js test/<시험>.js
   XPS_FITLOG=<파일>  fitOnce 가 끝날 때마다 한 줄(JSON): 자유 매개변수 이름·값·σ, χ², 반복 수, 정지 이유.
   XPS_LMLOG=<파일>   LM 한 번마다 한 줄: 반복별 χ²·받아들인 걸음의 λ·시도 수 (멈춤 규칙을 기록에서 다시 돌려 볼 때).
   XPS_NOSTALL=1      멈춤 규칙 끔(옛 동작과 비트 단위로 같음). trig = 그 피팅에서 푼 시도 갈래 수(0 이면 옛 경로 그대로).
   XPS_STALLCFG=JSON  멈춤 규칙 상수 바꿔 보기(예: {"W":10}).
   같은 시험을 고치기 전·후로 돌려 줄 번호(seq)로 맞대어 비교한다(test/fitlog_compare.js). */
"use strict";
const fs=require("fs"), path=require("path");
const out=process.env.XPS_FITLOG, lmOut=process.env.XPS_LMLOG;
if(out||lmOut||process.env.XPS_NOSTALL||process.env.XPS_STALLCFG){
  const C=require(path.join(__dirname,"../src/xps_core.js"));
  const tag=path.basename(process.argv[1]||"?"); let seq=0, sec="";
  const log0=console.log; console.log=(...a)=>{ const t=String(a[0]==null?"":a[0]); if(/^== /.test(t)) sec=t.slice(3,60); return log0.apply(console,a); };   // 지금 시험 절(== …) 기록
  if(process.env.XPS_NOSTALL) C._hook.stall=false;
  if(process.env.XPS_STALLCFG) Object.assign(C.LM_STALL,JSON.parse(process.env.XPS_STALLCFG));
  if(lmOut){ fs.writeFileSync(lmOut,""); C._hook.lm=o=>{ fs.appendFileSync(lmOut,JSON.stringify(Object.assign({tag,fit:seq},o))+"\n"); }; }
  if(out) C._hook.fit=(prob,r,opt)=>{
    const map=r.map, free=map.free; let sig=null;
    try{ const cv=C.covariance(prob,r); sig=cv.sig; }catch(e){ sig=null; }
    const p=free.map(f=>f.ci<0?prob.bg[f.k]:(f.k==="dE"||f.k==="r"||f.k==="q")?prob.comps[f.ci].doublet[f.k]:prob.comps[f.ci][f.k]);
    let ys=0; for(const t of prob.y) ys+=t;
    const rec={tag,seq,sec,n:prob.E.length,np:free.length,ys:+ys.toPrecision(10),names:free.map(f=>f.id+"|"+f.k),lo:free.map(f=>f.lo),hi:free.map(f=>f.hi),
      p,sig,chi2:r.chi2,chi2p1:r.chi2p1,iter:r.iter,iter1:r.iter1,conv:r.conv,stop:r.stop||null,stop1:r.stop1||null,kicks:r.kicks||0,trig:r.trig||0,pick:r.pick||null,branch:r.branch?{kind:r.branch.kind,reset:r.branch.reset.map(q=>q.id+"|"+q.k),cut:r.branch.cut}:null,stall:r.stall||null,opt:{maxit:(opt||{}).maxit||null,optimizer:r.optimizer||"lm"}};
    fs.appendFileSync(out,JSON.stringify(rec)+"\n");
  };
  if(!out) C._hook.fit=()=>{ seq++; }; else { const f0=C._hook.fit; C._hook.fit=(a,b,c)=>{ f0(a,b,c); seq++; }; }
}
