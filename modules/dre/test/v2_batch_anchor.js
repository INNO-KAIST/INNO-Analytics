/* V2 일괄: 참값 30세트(Lee 값 ±1 decade 무작위) × S-curve 1개 → 전역해 도달, 매개변수 오차 분포 */
const D=require('../src/dre_core.js'); const fs=require('fs');
const T=D.templateLee(); const M=D.compile(T.spec);
const FREE=['sG','sE','w3nr','k31','kQ22','kQ23'];   /* W2,nr 고정(³F₄ 수명 앎) */
const ALL=['sG','sE','w2nr','w3nr','k31','kQ22','kQ23'];
const rnd=D.mulberry(2026); const gauss=()=>{ let u=0,v=0; while(!u) u=rnd(); v=rnd(); return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v); };
const I=[...Array(45).keys()].map(i=>100*Math.pow(1000,i/44));
const rows=[];
for(let s=0;s<15;s++){
  const truth={c:8}; Object.entries(T.params).forEach(([n,p])=>truth[n]=p.v);
  ALL.forEach(n=>truth[n]*=Math.pow(10,(rnd()*2-1)));
  const sw=D.sweep(M,truth,I,1,1064); const n3=sw.n.map(n=>n[2]);
  const Sc=3e6/Math.max(...n3)*0.3, bg=40;
  const y=n3.map(v=>{ const y0=Sc*v; return Math.max(y0+bg+gauss()*Math.sqrt(y0+bg),1)*Math.exp(0.05*gauss())-bg; });
  const params={}; Object.entries(T.params).forEach(([n,p])=>{ params[n]={v:(n==='w2nr'?truth.w2nr:p.v),fixed:!FREE.includes(n),lo:T.params[n].v/1000,hi:T.params[n].v*1000}; });
  const tau=D.decayTime(M,truth,I[0],1064,0).t*Math.exp(0.03*gauss());
  const B=D.buildProblem({spec:T.spec,nm:1064,beam:"peak",params,samples:[{id:"S",vars:{c:8},data:[{type:"scurve",ch:T.spec.channels[0].id,dir:1,I,y,sbg:Math.sqrt(bg),k:3,relErr:0.05},{type:"decay",ch:T.spec.channels[0].id,I:[I[0]],t:[tau],relErr:0.03}]}]});
  const thT=Float64Array.from(FREE.map(n=>Math.log10(truth[n])));
  const chiT=D.lm(B,thT,{maxit:200}).chi2;          /* 참값에서 출발한 국소 최적 = 전역해 기준 */
  const t0=Date.now(); const ms=D.multistart(B,{starts:12,seed:s+1,maxit:200});
  const best=ms[0];
  const reach=ms.filter(r=>r.chi2<=Math.min(chiT,best.chi2)+0.5).length/ms.length;
  const slopeMax=Math.max(...D.slopes(I,n3));
  rows.push({s,slopeMax,chiT,best:best.chi2,reach,time:(Date.now()-t0)/1000,err:FREE.map((n,j)=>best.theta[j]-thT[j])});
  console.log(s,'slope',slopeMax.toFixed(1),'χ² truth-local',chiT.toFixed(1),'best',best.chi2.toFixed(1),'reach',Math.round(reach*100)+'%',((Date.now()-t0)/1000).toFixed(1)+'s','err',rows[rows.length-1].err.map(v=>v.toFixed(2)).join(' '));
}
fs.writeFileSync('v2_batch_anchor.json',JSON.stringify(rows,null,1));
const med=a=>{ const b=[...a].sort((x,y)=>x-y); return b[Math.floor(b.length/2)]; };
console.log('median reach',med(rows.map(r=>r.reach)),' found-global(best≤truth-local+0.5)',rows.filter(r=>r.best<=r.chiT+0.5).length+'/'+rows.length);
FREE.forEach((n,j)=>console.log(n,'|err| median',med(rows.map(r=>Math.abs(r.err[j]))).toFixed(2),' 90%',[...rows.map(r=>Math.abs(r.err[j]))].sort((a,b)=>a-b)[13].toFixed(2)));
