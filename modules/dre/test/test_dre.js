/* DRE 엔진 회귀 시험 — node test/test_dre.js  (참조 데이터: lee_ref.json, v3_ref.json 은 scipy 독립 구현으로 생성) */
const D=require('../src/dre_core.js'); const ref=require('./lee_ref.json'); const bi=require('./v3_ref.json');
let pass=0,fail=0; const ok=(c,m)=>{ if(c) pass++; else { fail++; console.log('FAIL',m); } };
const T=D.templateLee(); const M=D.compile(T.spec);
const base=()=>{ const e={}; Object.entries(T.params).forEach(([n,p])=>e[n]=p.v); return e; };
/* 1. 식 해석기 */
const ex=D.compileExpr("k31*c^2+ (1-b32)*w3r/2 - -1"); ok(Math.abs(ex({k31:2,c:3,b32:0.5,w3r:4})-(18+1+1))<1e-12,'expr');
let thrown=false; try{ D.compileExpr("a+*b"); }catch(e){ thrown=true; } ok(thrown,'expr error');
/* 2. V1 — Lee 6개 시료 정상상태·동역학 */
let wS=0,wK=0;
for(const [k,r] of Object.entries(ref)){
  const env=Object.assign(base(),{c:r.c,w2nr:r.w2nr,w3nr:r.w3nr});
  const up=D.sweep(M,env,r.I,1,1064), dn=D.sweep(M,env,r.I,-1,1064); ok(!up.fails&&!dn.fails,k+' conv');
  r.steady.forEach((q,i)=>{ const [n2,n3]=q[0]; const e=Math.max(Math.abs(up.n[i][1]/n2-1),Math.abs(up.n[i][2]/n3-1),Math.abs(dn.n[i][2]/n3-1)); wS=Math.max(wS,e); });
  for(const [I,kin] of Object.entries(r.kin)){ const res=D.integrate(M,env,+I,1064,2.0,null,{rtol:1e-7,tOut:kin.t}); kin.t.forEach((t,j)=>{ if(kin.n3[j]>1e-12) wK=Math.max(wK,Math.abs(res.outN[j][2]/kin.n3[j]-1)); }); }
}
ok(wS<1e-8,`V1 steady worst ${wS}`); ok(wK<1e-5,`V1 kinetic worst ${wK}`);
/* 3. V3 — 쌍안정 가지 선택 */
{ const env=Object.assign(base(),{c:bi.c},bi.e); const I=bi.res.map(r=>r.I); const up=D.sweep(M,env,I,1,1064), dn=D.sweep(M,env,I,-1,1064); let nb=0;
  bi.res.forEach((r,i)=>{ const st=r.roots.map((x,j)=>({n3:x[1],s:r.stable[j]})).filter(q=>q.s).map(q=>q.n3).sort((a,b)=>a-b);
    if(st.length>1) nb++; ok(Math.abs(up.n[i][2]/st[0]-1)<1e-6&&Math.abs(dn.n[i][2]/st[st.length-1]-1)<1e-6,`V3 @${r.I}`); }); ok(nb>=10,'V3 bistable points'); }
/* 4. 스케일 대칭 (정확) */
{ const e=Object.assign(base(),{c:8}), e2=Object.assign({},e); ['sG','sE','w2r','w2nr','w3r','w3nr','k31','kQ22','kQ23'].forEach(n=>e2[n]*=7.3);
  const a=D.sweep(M,e,[2e3,7e3,3e4],1,1064), b=D.sweep(M,e2,[2e3,7e3,3e4],1,1064); ok(a.n.every((n,i)=>Math.abs(n[2]/b.n[i][2]-1)<1e-8),'scale symmetry');
  const ta=D.riseTime(M,e,7e3,1064,0,0.95).t, tb=D.riseTime(M,e2,7e3,1064,0,0.95).t; ok(Math.abs(ta/tb/7.3-1)<1e-3,`rise time scales 1/k (${ta/tb})`); }
/* 5. IFT 감도 = 유한차분 */
{ const env=Object.assign(base(),{c:8}), I=7e3, phi=I/D.photonEnergyJ(1064), names=['sG','sE','w2nr','w3nr','k31','kQ22','kQ23'];
  const n=D.sweep(M,env,[I],1,1064).n[0], S=D.sensitivity(M,env,n,phi,names);
  names.forEach((nm,j)=>{ const h=1e-4, e1=Object.assign({},env), e2=Object.assign({},env); e1[nm]*=Math.pow(10,h); e2[nm]*=Math.pow(10,-h);
    const fd=(D.sweep(M,e1,[I],1,1064).n[0][2]-D.sweep(M,e2,[I],1,1064).n[0][2])/(2*h); ok(Math.abs(S[2*names.length+j]-fd)<=1e-4*Math.abs(fd)+1e-18,'IFT '+nm); }); }
/* 6. 잔차 Jacobian(계층 배율·절단 포함) = 유한차분 ; 가우시안 빔 포함 */
for(const beam of ["peak","gauss"]){
  const env=Object.assign(base(),{c:8}); const I=[...Array(25).keys()].map(i=>300*Math.pow(200,i/24));
  const y=D.observe(M,env,I,1,1064,0,{beam}).y.map((v,i)=>v*1e6*(1+0.03*Math.sin(i)));
  const FREE=['sG','sE','w2nr','w3nr','k31','kQ22','kQ23']; const params={}; Object.entries(T.params).forEach(([n,p])=>params[n]={v:p.v*(FREE.includes(n)?1.3:1),lo:p.v/100,hi:p.v*100,fixed:!FREE.includes(n)});
  const B=D.buildProblem({spec:T.spec,nm:1064,beam,params,samples:[{id:"S",vars:{c:8},data:[{type:"scurve",ch:T.spec.channels[0].id,dir:1,I,y,sbg:0.05,k:3,relErr:0.05}]}]});
  const th=D.thetaFromParams(B), R=D.residuals(B,th,true); let w=0;
  for(let j=0;j<FREE.length;j++){ const h=1e-5, a=Float64Array.from(th), b=Float64Array.from(th); a[j]+=h; b[j]-=h; const ra=D.residuals(B,a,false).r, rb=D.residuals(B,b,false).r;
    for(let i=0;i<ra.length;i++){ const fd=(ra[i]-rb[i])/(2*h); w=Math.max(w,Math.abs(fd-R.J[i][j])/(Math.abs(fd)+1e-2)); } }
  ok(w<1e-4,`residual J vs FD (${beam}) ${w}`);
  /* 잡음 없는 데이터에서 LM 이 χ² ≈ 0 까지 */
  const yc=D.observe(M,env,I,1,1064,0,{beam}).y.map(v=>v*1e6); B.problem.samples[0].data[0].y=yc;
  const r=D.lm(B,th,{maxit:300}); ok(r.chi2<1e-6,`LM noise-free (${beam}) χ²=${r.chi2}`);
}
/* 7. 가우시안 빔 적분 수렴 (격자 40 vs 160 점/decade) */
{ const env=Object.assign(base(),{c:8}), I=[1e3,5e3,7e3,2e4,8e4];
  const a=D.observe(M,env,I,1,1064,0,{beam:"gauss"}).y, b=D.observe(M,env,I,1,1064,0,{beam:"gauss",perDec:1920}).y;
  ok(a.every((v,i)=>Math.abs(v/b[i]-1)<(i===2?1e-2:1e-4)),'gauss beam grid convergence '+a.map((v,i)=>(v/b[i]-1).toExponential(1)).join(',')); }
/* 8. 2준위 감광 + 3준위 활성 UC: 저파워 기울기 1(감광)·2(2광자 준위) */
{ const spec={species:[{id:"S",levels:["g","e"],cvar:"cs"},{id:"A",levels:["g","m","u"],cvar:"ca"}],
   procs:[{type:"abs",sp:"S",from:"g",to:"e",sigma:"sS"},{type:"dec",sp:"S",from:"e",to:"g",k:"kS"},
          {type:"et",d:["S","e","g"],a:["A","g","m"],W:"W1",scale:"x"},{type:"et",d:["S","e","g"],a:["A","m","u"],W:"W2",scale:"x"},
          {type:"dec",sp:"A",from:"m",to:"g",k:"km"},{type:"dec",sp:"A",from:"u",to:"g",k:"ku"}],
   channels:[{id:"S",terms:[["S","e",1]]},{id:"u",terms:[["A","u",1]]}]};
  const M2=D.compile(spec); const env={sS:1e-20,kS:1e3,W1:1e5,W2:1e5,km:1e3,ku:1e4,cs:20,ca:2};
  const I=[0.1,0.2,1e5,2e5]; const y=[0,1].map(k=>D.observe(M2,env,I,1,980,k).y);
  const s=(yy,a,b)=>Math.log(yy[b]/yy[a])/Math.log(I[b]/I[a]);
  ok(Math.abs(s(y[0],0,1)-1)<0.01&&Math.abs(s(y[1],0,1)-2)<0.01,`UC low-power slopes ${s(y[0],0,1).toFixed(3)}, ${s(y[1],0,1).toFixed(3)}`);
  ok(s(y[1],2,3)<0.9,'UC saturation slope < 1');
  /* 에너지(여기) 보존: 흡수 = 감쇠 합 (정상상태) */
  const c=D.coeffs(M2,env), x=D.molFrac(M2,env), phi=1e3/D.photonEnergyJ(980); const n=D.steady(M2,c,x,phi,null).n;
  const Nabs=0.2*env.sS*phi*n[0]; ok(n.every(v=>v>=0),'populations ≥ 0'); ok(Math.abs(n[0]+n[1]-1)<1e-12&&Math.abs(n[2]+n[3]+n[4]-1)<1e-12,'conservation per species');
}
/* 9. 상승 시간 (임계 부근) */
{ const env=Object.assign(base(),{c:8}); const t=D.riseTime(M,env,7e3,1064,0,0.95).t; ok(t>0.6&&t<0.68,`t95(7 kW/cm²)=${t}`);
  const res=D.integrate(M,env,7e3,1064,0.5,null,{tOut:[0.5]}); const ss=D.sweep(M,env,[7e3],1,1064).n[0][2]; const q=res.outN[0][2]/ss; ok(q<0.9,`0.5 s 적분/정상 = ${q.toFixed(3)} (Lee 코드 근사의 오차)`); console.log('n3(0.5 s)/n3(∞) @7 kW/cm² =',q.toFixed(3)); }
console.log(`pass ${pass}, fail ${fail}`); process.exit(fail?1:0);
