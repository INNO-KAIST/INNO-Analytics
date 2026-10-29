/* LM 멈춤 규칙(SPEC §3.3b) 무작위 비교 — 같은 합성 스펙트럼을 두 엔진으로 피팅해 χ²·반복·시간을 줄마다 적는다.
   node test/stall_stress.js <엔진.js> <출력.jsonl> [첫 seed=1] [개수=320]   (XPS_NOSTALL=1 이면 멈춤 규칙을 끈 옛 동작)
   비교: node test/stall_stress_cmp.js <옛.jsonl> <새.jsonl>
   스펙트럼: 성분 1–6 (모든 모양 계열 + Cookbook 2026 새 모양 · ST/SB 수식어 · doublet), 배경 linear/shirley/smart/tougaard/active,
   가중 poisson/casa/const, 세기 10²–10⁵ (Poisson 잡음). 자유 매개변수: 'tight' = BE·넓이·주 폭, 'loose'(약 25 %) = 모양 매개변수까지 모두
   (퇴화·의미 없는 매개변수가 섞인 ill-posed 피팅). 생성은 항상 이 파일의 난수로(엔진과 무관). */
"use strict";
const fs=require("fs"), path=require("path");
const C=require(path.resolve(process.argv[2])); if(process.env.XPS_NOSTALL) C._hook.stall=false;
const out=process.argv[3]; const s0=+(process.argv[4]||1), N=+(process.argv[5]||320);
const G=require(path.join(__dirname,"../src/xps_core.js")).rng;   // 생성 난수(같은 seed → 같은 스펙트럼)
const SH=[["GLS",{m:0.3}],["GLP",{m:0.3}],["V",{}],["DS",{alpha:0.12}],["LA",{alpha:1.2,beta:2.0}],["GLT",{m:0.2,T:0.15,tau:0.8}],["AGT",{m:0.3,ET:0.05,TM:0.2,CT:0.05}],
  ["TLA",{alpha:1.3,mu:2}],["GLN",{m:0.3,tail:"T",tk:2}],["SGLN",{m:0.3,tail:"PHI",tj:1,tk:2}],["GLS",{m:0.2,mod:"ST",stMu:1.5,stG:0.3}],["LA",{alpha:1,beta:1,mod:"SB"}],["V",{mod:"ST",stMu:2,stG:0}],["GLN",{m:0.2}]];
const BG=["linear","shirley","smart","tougaard","active"];
fs.writeFileSync(out,"");
for(let s=s0;s<s0+N;s++){
  const R=G(1000+s), U=(a,b)=>a+(b-a)*R.u(), I=k=>Math.floor(R.u()*k);
  const n=150+I(250), st=U(0.04,0.1), E=Float64Array.from({length:n},(_,i)=>280+st*i), ctx={emin:E[0],emax:E[n-1]};
  const K=1+I(6), loose=R.u()<0.25, comps=[], cons={};
  const span=E[n-1]-E[0];
  for(let k=0;k<K;k++){ const [sh,o]=SH[I(SH.length)]; const F=U(0.7,1.8);
    const c=C.newComp(sh,E[0]+span*U(0.25,0.75),Object.assign({id:"c"+k,F},o)); c.area=U(0.2,1)*Math.pow(10,U(2,4.5))*F;
    if(sh==="V"){ c.wG=F*0.8; c.wL=F*0.3; } if(sh==="DS"){ c.wL=F*0.25; c.wG=F*0.6; } if(sh==="LA"||sh==="TLA"){ c.wG=F*0.4; }
    if(sh==="GLN"||sh==="SGLN"){ c.wG=F*0.3; }
    if(R.u()<0.3&&sh!=="TP") c.doublet={dE:U(0.6,5),r:0.5,q:1};
    comps.push(c);
    const main=sh==="V"||sh==="DS"?"wL":sh==="MEAS"?"wG":"F"; cons[c.id]={};
    C.paramNames(c).forEach(q=>{ if(q==="be"||q==="area"||q===main) return;
      if(q==="dE"||q==="r"||q==="q"){ if(!loose||q==="q") cons[c.id][q]={mode:"fixed"}; return; }
      if(!loose||c[q]==null||q==="stG") cons[c.id][q]={mode:"fixed"}; else cons[c.id][q]={mode:"free"}; }); }   // ST 의 γ(stG) 는 풀지 않음: γ → 0 이면 격자가 40 만 점까지 잘게 — 비용 문제(이 시험의 대상 아님)
  if(R.u()<0.12){ const c=C.newComp("TP",E[0]+span*U(0.4,0.6),{id:"tp",F:span*U(0.3,0.6),ta:1,tb:3,wL:0.5}); c.area=comps.reduce((a,c)=>a+c.area,0)*U(0.3,1); comps.push(c); cons.tp={ta:{mode:"fixed"},tb:{mode:"fixed"},wL:loose?{mode:"free"}:{mode:"fixed"}}; }
  if(R.u()<0.15&&comps.length>1&&comps[1].F!=null&&comps[0].F!=null){ cons[comps[1].id].F={mode:"link",ref:comps[0].id,k:1}; }
  const mu=new Float64Array(n); comps.forEach(c=>{ const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) mu[i]+=v[i]; });
  let mx=0; for(const t of mu) mx=Math.max(mx,t); const b0=mx*U(0.05,0.6), slope=U(-0.2,0.2)*b0/span, kk=U(0,0.3);
  const Q=new Float64Array(n); for(let i=1;i<n;i++) Q[i]=Q[i-1]+0.5*(mu[i]+mu[i-1])*st; const Qt=Q[n-1]||1;
  const cf=Math.pow(10,U(-1,1.5));   // 세기 척도(잡음 수준)
  const y=Float64Array.from(mu,(v,i)=>{ const lam=(v+b0+slope*(E[i]-E[0])+kk*b0*Q[i]/Qt)*cf; const g=R.g(); return Math.max(0,lam+Math.sqrt(Math.max(lam,0))*g)/cf; });
  const bgT=BG[I(BG.length)]; const bg={type:bgT,n:5}; if(bgT==="tougaard"){ bg.B=2866; bg.C=1643; }
  const wr=R.u(), weight=wr<0.1?"casa":wr<0.2?"const":"poisson";
  const start=C.clone(comps).map(c=>{ c.be+=0.15*R.g(); ["F","wL","wG"].forEach(k=>{ if(c[k]>0) c[k]*=U(0.85,1.15); }); return c; });
  let P; try{ P=C.makeProblem({be:E,y,c:cf},{lo:E[0],hi:E[n-1],bg,comps:start,cons,weight}); }catch(e){ continue; }
  const t0=process.hrtime.bigint(); let r, err=null; try{ r=C.fitOnce(P,{}); }catch(e){ err=e.message; }
  const ms=Number(process.hrtime.bigint()-t0)/1e6;
  const rec={seed:s,K:comps.length,shapes:comps.map(c=>c.shape+(c.mod?"+"+c.mod:"")+(c.doublet?"(d)":"")).join(" "),bg:bgT,weight,loose,n,ms};
  if(err) rec.err=err; else { const S=C.fitStats(P,r); const map=r.map; Object.assign(rec,{np:map.free.length,chi2:r.chi2,chi2p1:r.chi2p1,chi2d:S.chi2d,iter:r.iter,iter1:r.iter1,stop:r.stop||null,stop1:r.stop1||null,kicks:r.kicks||0,
    flags:r.stall?r.stall.flags.map(q=>q.id+"|"+q.k+":"+q.why.join("+")):null,names:map.free.map(f=>f.id+"|"+f.k),p:map.free.map(f=>f.ci<0?P.bg[f.k]:(f.k==="dE"||f.k==="r"||f.k==="q")?P.comps[f.ci].doublet[f.k]:P.comps[f.ci][f.k])});
    try{ rec.sig=C.covariance(P,r).sig; }catch(e){} }
  fs.appendFileSync(out,JSON.stringify(rec)+"\n");
}
