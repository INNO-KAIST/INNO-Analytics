/* V2 파일럿 (수정 엔진) — 합성 S-curve(참값 앎): 전역해 도달률·식별성·수명 데이터의 효과 */
const D=require('../src/dre_core.js'); const fs=require('fs');
const T=D.templateLee(); const M=D.compile(T.spec);
const rnd=D.mulberry(7); const gauss=()=>{ let u=0,v=0; while(!u) u=rnd(); v=rnd(); return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v); };
const truthEnv={c:8}; Object.entries(T.params).forEach(([n,p])=>truthEnv[n]=p.v);
const I=[...Array(45).keys()].map(i=>300*Math.pow(60000/300,i/44));
const n3=D.sweep(M,truthEnv,I,1,1064).n.map(n=>n[2]);
const S=3e6,bg=40; const ynet=n3.map(v=>{ const y0=S*v; return Math.max(y0+bg+gauss()*Math.sqrt(y0+bg),1)*Math.exp(0.05*gauss())-bg; });
const tau800=D.decayTime(M,truthEnv,1000,1064,0).t*Math.exp(0.03*gauss());
const out={};
function run(label,FREE,extra,opt){
  const params={}; Object.entries(T.params).forEach(([n,p])=>{ params[n]={v:p.v,fixed:!FREE.includes(n),lo:p.v/100,hi:p.v*100}; });
  const data=[{type:"scurve",ch:T.spec.channels[0].id,dir:1,I,y:ynet,sbg:Math.sqrt(bg),k:3,relErr:0.05}].concat(extra||[]);
  const B=D.buildProblem({spec:T.spec,nm:1064,beam:"peak",params,samples:[{id:"S1",vars:{c:8},data}]});
  const th=D.thetaFromParams(B);
  let t0=Date.now(); const ms=D.multistart(B,{starts:opt.starts||40,seed:3,maxit:200});
  const best=ms[0], info=D.dofInfo(B,best), thr=3.84*info.s2;
  const r={label,time:(Date.now()-t0)/1000,best:best.chi2,dof:info.dof,succ:ms.filter(q=>q.chi2-best.chi2<=0.5).length/ms.length,
    waterfall:ms.map(q=>q.chi2),err:{},prof:{},fim:null};
  FREE.forEach((n,j)=>{ r.err[n]=best.theta[j]-th[j]; const p=D.profile(B,best,j,{span:2,npts:17,thr,maxit:120}); r.prof[n]={lo:p.lo-th[j],hi:p.hi-th[j],w:p.hi-p.lo,openLo:p.openLo,openHi:p.openHi}; });
  const ev=D.fim(B,best.theta); r.fim=ev.map(e=>Math.log10(Math.max(Math.abs(e.val),1e-30))); r.sloppy=Object.fromEntries(FREE.map((n,j)=>[n,+ev[0].vec[j].toFixed(3)]));
  console.log(`\n[${label}] ${r.time.toFixed(1)}s  χ²=${best.chi2.toFixed(1)} (dof ${info.dof})  전역해 도달 ${Math.round(r.succ*100)}%`);
  FREE.forEach(n=>console.log(`  ${n.padEnd(5)} 오차 ${r.err[n].toFixed(2).padStart(6)}  95% [${r.prof[n].lo.toFixed(2)}, ${r.prof[n].hi.toFixed(2)}] ${r.prof[n].openLo?'↓열림 ':''}${r.prof[n].openHi?'↑열림':''}`));
  console.log('  FIM log10 eig:',r.fim.map(v=>v.toFixed(1)).join(' '),' | 가장 sloppy:',JSON.stringify(r.sloppy));
  out[label]=r; return r;
}
run("A. S-curve 1개, 7개 자유",['sG','sE','w2nr','w3nr','k31','kQ22','kQ23'],null,{});
run("B. + 800 nm 감쇠 1/e 1점(1 kW/cm²)",['sG','sE','w2nr','w3nr','k31','kQ22','kQ23'],[{type:"decay",ch:T.spec.channels[0].id,I:[1000],t:[tau800],relErr:0.03}],{starts:24});
run("C. S-curve 1개, ³F₄ 비복사 W₂,nr 고정(수명 앎)",['sG','sE','w3nr','k31','kQ22','kQ23'],null,{});
fs.writeFileSync('v2_pilot.json',JSON.stringify({tau800,out},null,1));
