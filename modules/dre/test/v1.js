/* V1: Lee 2021 모델 재현 — JS 엔진 vs Python(scipy) 독립 구현 */
const D=require('../src/dre_core.js'); const ref=require('./lee_ref.json');
let pass=0,fail=0; const ok=(c,m)=>{ if(c) pass++; else { fail++; console.log('FAIL',m); } };
const T=D.templateLee(); const M=D.compile(T.spec);
let worstSS=0, worstKin=0;
for(const [k,r] of Object.entries(ref)){
  const env={}; Object.entries(T.params).forEach(([n,p])=>env[n]=p.v); env.c=r.c; env.w2nr=r.w2nr; env.w3nr=r.w3nr;
  const up=D.sweep(M,env,r.I,1,1064), dn=D.sweep(M,env,r.I,-1,1064);
  ok(up.fails===0&&dn.fails===0,`${k} converge`);
  r.steady.forEach((roots,i)=>{ const [n2,n3]=roots[0]; const a=up.n[i];
    const e=Math.max(Math.abs(a[1]/n2-1),Math.abs(a[2]/n3-1)); worstSS=Math.max(worstSS,e);
    ok(e<1e-6,`${k} I=${r.I[i]} ss n2 ${a[1]} vs ${n2}, n3 ${a[2]} vs ${n3}`);
    ok(Math.abs(dn.n[i][2]/n3-1)<1e-6,`${k} down-sweep same (no hysteresis)`); });
  for(const [I,kin] of Object.entries(r.kin)){
    const res=D.integrate(M,env,+I,1064,2.0,null,{rtol:1e-7,tOut:kin.t});
    kin.t.forEach((t,j)=>{ const n3=res.outN[j][2]; const rr=kin.n3[j]; if(rr>1e-12){ const e=Math.abs(n3/rr-1); worstKin=Math.max(worstKin,e); ok(e<2e-3,`${k} kin I=${I} t=${t} ${n3} vs ${rr}`); } });
  }
}
/* 스케일 대칭: 모든 계수 ×k → 정상상태 동일 */
{ const env={}; Object.entries(T.params).forEach(([n,p])=>env[n]=p.v); env.c=8; const e2=Object.assign({},env);
  ['sG','sE','w2r','w2nr','w3r','w3nr','k31','kQ22','kQ23'].forEach(n=>e2[n]*=7.3);
  const I=[2e3,7e3,3e4]; const a=D.sweep(M,env,I,1,1064), b=D.sweep(M,e2,I,1,1064);
  ok(I.every((_,i)=>Math.abs(a.n[i][2]/b.n[i][2]-1)<1e-8),'scale symmetry'); }
/* IFT 감도 vs 유한차분 */
{ const env={}; Object.entries(T.params).forEach(([n,p])=>env[n]=p.v); env.c=8; const I=7e3, phi=I/D.photonEnergyJ(1064);
  const names=['sG','sE','w2nr','w3nr','k31','kQ22','kQ23']; const n=D.sweep(M,env,[I],1,1064).n[0];
  const S=D.sensitivity(M,env,n,phi,names);
  names.forEach((nm,j)=>{ const h=1e-4; const e1=Object.assign({},env), e2=Object.assign({},env); e1[nm]*=Math.pow(10,h); e2[nm]*=Math.pow(10,-h);
    const a=D.sweep(M,e1,[I],1,1064).n[0][2], b=D.sweep(M,e2,[I],1,1064).n[0][2]; const fd=(a-b)/(2*h);
    ok(Math.abs(S[2*names.length+j]-fd)<=1e-4*Math.abs(fd)+1e-16,`IFT ${nm} ${S[2*names.length+j]} vs FD ${fd}`); }); }
/* 임계 부근 상승 시간(0.5 s 적분 근사 오류 재현) */
{ const env={}; Object.entries(T.params).forEach(([n,p])=>env[n]=p.v); env.c=8; M.channels;
  const r=D.riseTime(M,env,7e3,1064,0,0.95); ok(r.t>0.5&&r.t<0.8,`t95 at 7 kW/cm² = ${r.t}`); console.log('t95(7 kW/cm2) =',r.t.toFixed(3),'s'); }
console.log(`worst steady rel err ${worstSS.toExponential(2)}, worst kinetic rel err ${worstKin.toExponential(2)}`);
console.log(`pass ${pass}, fail ${fail}`); process.exit(fail?1:0);
