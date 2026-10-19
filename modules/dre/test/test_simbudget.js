/* 화면 시뮬레이션 상한 시험 — node test/test_simbudget.js
   observeMulti(스윕 한 번 → 모든 채널) = observe(채널마다) · 엉뚱한 값(sG 1e-10)에서 시간 상한·실패 집계 · Lee 2021 기본값 결과 유지 */
const D=require('../src/dre_core.js');
let pass=0,fail=0; const ok=(c,m)=>{ if(c) pass++; else { fail++; console.log('FAIL',m); } };
const envOf=T=>{ const e={}; Object.entries(T.params).forEach(([n,p])=>e[n]=p.v); return Object.assign(e,T.samples[0].vars); };
/* 1. 같은 값: 여러 채널 한 번 = 채널마다 */
for(const [T,I] of [[D.templateLee(),D.logGrid(100,2e5,40)],[D.templateTm(6),D.logGrid(1e3,1e7,40)]]) for(const beam of ['peak','gauss']){
  const M=D.compile(T.spec), env=envOf(T), ch=M.channels.map((_,k)=>k);
  const mu=D.observeMulti(M,env,I,1,T.nm,ch,{beam});
  let worst=0; ch.forEach(k=>{ const o=D.observe(M,env,I,1,T.nm,k,{beam}); o.y.forEach((v,i)=>{ worst=Math.max(worst,Math.abs(v-mu.Y[k][i])/Math.max(1e-300,Math.abs(v))); }); });
  ok(worst===0,`${T.name.slice(0,20)} ${beam} multi==single (${worst})`); ok(mu.fails===0&&!mu.aborted,'default converges');
}
/* 2. Lee 기본값: threshold 6.87 kW/cm², 최대 기울기 ≈ 32 (화면 격자 0.1–200 kW/cm², 40/decade) */
{ const T=D.templateLee(), M=D.compile(T.spec), I=D.logGrid(100,2e5,40);
  const o={beam:'peak',maxFails:8,deadline:Date.now()+1500,stOpt:{maxMs:60}};
  const u=D.observeMulti(M,envOf(T),I,1,T.nm,[0],o); const s=D.slopes(I,u.Y[0]); let i=0; s.forEach((v,j)=>{ if(v>s[i]) i=j; });
  ok(Math.abs(I[i]/1e3-6.87)<0.05,`threshold ${(I[i]/1e3).toFixed(3)} kW/cm²`); ok(Math.abs(s[i]-32)<1.5,`s_max ${s[i].toFixed(2)}`); }
/* 3. 오타 값: 시간 상한과 실패 보고 */
{ const T=D.templateLee(), M=D.compile(T.spec), I=D.logGrid(100,2e5,40);
  for(const v of [1e-10,1e10]){ const env=Object.assign(envOf(T),{sG:v}); const t0=Date.now();
    const o=()=>({beam:'peak',maxFails:8,deadline:Date.now()+1500,stOpt:{maxMs:60}});
    const u=D.observeMulti(M,env,I,1,T.nm,[0],o()), d=D.observeMulti(M,env,I,-1,T.nm,[0],o()); const dt=Date.now()-t0;
    ok(dt<3500,`sG=${v} took ${dt} ms`); ok(u.fails>0&&u.aborted,`sG=${v} failures reported (${u.fails}, aborted ${u.aborted})`);
    ok(u.Y[0].some(x=>!isFinite(x))||u.done<u.total,'uncomputed points are NaN'); }
  /* 상한 없는 기본 sweep 은 예전과 같다(피팅용): 중단하지 않는다 */
  const sw=D.sweep(M,Object.assign(envOf(T),{sG:1e-10}),I.slice(0,20),1,T.nm); ok(!sw.aborted&&sw.n.every(Boolean),'unbounded sweep unchanged'); }
console.log(`sim budget tests ${pass} pass / ${fail} fail`); if(fail) process.exit(1);
