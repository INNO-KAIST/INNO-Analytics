// Tm 5/6준위: (1) 3준위 극한에서 Lee 템플릿과 일치 (2) 점유 보존·음수 없음 (3) scipy 독립 구현 대조용 출력
const D=require('../src/dre_core.js'); const fs=require('fs');
let pass=0,fail=0; const ok=(c,m)=>{ if(c) pass++; else { fail++; console.log('✗ '+m); } };
const envOf=(T,extra)=>Object.assign(Object.fromEntries(Object.entries(T.params).map(([k,p])=>[k,p.v])),T.samples[0].vars,extra||{});
/* (1) 3준위 극한 */
{ const L=D.templateLee(), T5=D.templateTm(5); const c=8;
  const eL=envOf(L,{c}); const s31=eL.k31*c*c, q22=eL.kQ22*c**3/(c*c+4.3**2), q23=eL.kQ23*c**3/(c*c+4.3**2), N=1.42e22*c/100;
  const e5=envOf(T5,{c,sG:eL.sG,sE1:eL.sE,sE2:0,w2r:eL.w2r,w2nr:eL.w2nr,w35:1e10,b51:0,w3r:eL.w3r,b32:eL.b32,w3nr:eL.w3nr,cCR:s31/N,cUC23:q23/N,cUC22:q22/N});
  const ML=D.compile(L.spec), M5=D.compile(T5.spec); const I=D.logGrid(1e3,2e5,20);
  const uL=D.sweep(ML,eL,I,1,1064).n, u5=D.sweep(M5,e5,I,1,1064).n;
  let mx=0; I.forEach((_,i)=>{ mx=Math.max(mx,Math.abs(u5[i][1]/uL[i][1]-1),Math.abs(u5[i][3]/uL[i][2]-1)); });
  console.log('3준위 극한 정상상태 최대 상대차',mx.toExponential(2)); ok(mx<1e-4,'3준위 극한 정상상태');
  const tL=D.riseTime(ML,eL,2e4,1064,0,0.9).t, t5=D.riseTime(M5,e5,2e4,1064,0,0.9).t; console.log('상승 t90',tL,t5); ok(Math.abs(t5/tL-1)<1e-3,'3준위 극한 상승시간');
  const hL=D.sweep(ML,eL,I,-1,1064).n, h5=D.sweep(M5,e5,I,-1,1064).n; let md=0; I.forEach((_,i)=>md=Math.max(md,Math.abs(h5[i][3]/hL[i][2]-1))); ok(md<1e-4,'3준위 극한 내림 가지'); }
/* (2) 보존·양수, (3) 기준값 출력 */
const ref={};
for(const L of [5,6]){ const T=D.templateTm(L); const M=D.compile(T.spec); const e=envOf(T,{c:30});
  const I=D.logGrid(1e3,1e7,10); const u=D.sweep(M,e,I,1,1064).n, d=D.sweep(M,e,I,-1,1064).n;
  const bad=u.concat(d).some(n=>Math.abs(n.reduce((a,b)=>a+b,0)-1)>1e-10||n.some(v=>v< -1e-14)); ok(!bad,`${L}준위 보존·양수`);
  const tOut=[1e-5,1e-4,1e-3,1e-2]; const r=D.integrate(M,e,1e5,1064,1e-2,null,{tOut});
  ref[L]={I,up:u,dn:d,kin:{I:1e5,t:r.outT,n:r.outN},params:e,procs:T.spec.procs,levels:T.spec.species[0].levels}; }
fs.writeFileSync(__dirname+'/tm5_ref_js.json',JSON.stringify(ref));
console.log(`Tm5/6 시험 ${pass} 통과 / ${fail} 실패`); process.exit(fail?1:0);
