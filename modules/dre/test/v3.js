/* V3: 쌍안정 — up 스윕 = 낮은 안정 가지, down 스윕 = 높은 안정 가지 (Python 전근 + 안정성 판정과 대조) */
const D=require('../src/dre_core.js'); const ref=require('./v3_ref.json');
const T=D.templateLee(); const M=D.compile(T.spec);
const env={c:ref.c}; Object.entries(T.params).forEach(([n,p])=>env[n]=p.v); Object.assign(env,ref.e);
const I=[...Array(60).keys()].map(i=>2e4*Math.pow(30,i/59)).concat(ref.res.map(r=>r.I));
I.sort((a,b)=>a-b);
const up=D.sweep(M,env,I,1,1064), dn=D.sweep(M,env,I,-1,1064);
let pass=0,fail=0; const ok=(c,m)=>{ if(c) pass++; else { fail++; console.log('FAIL',m); } };
ok(up.fails===0&&dn.fails===0,'converged');
ref.res.forEach(r=>{
  const k=I.indexOf(r.I); const st=r.roots.map((x,i)=>({n3:x[1],s:r.stable[i]})).filter(q=>q.s).map(q=>q.n3).sort((a,b)=>a-b);
  const lo=st[0], hi=st[st.length-1];
  ok(Math.abs(up.n[k][2]/lo-1)<1e-6,`up @${r.I.toFixed(0)}: ${up.n[k][2]} vs low ${lo}`);
  ok(Math.abs(dn.n[k][2]/hi-1)<1e-6,`down @${r.I.toFixed(0)}: ${dn.n[k][2]} vs high ${hi}`);
});
console.log(`V3 bistable: pass ${pass}, fail ${fail}`); process.exit(fail?1:0);
