// 적합성 검사 시험: (T1) 모델 자신이 만든 데이터 → 통과, (T3) 구조가 다른 데이터 → 실패, (T2) WJ-0128c 실측 → 실패
const D=require('../src/dre_core.js'); const fs=require('fs');
const K=2.6813e5, V4=process.env.V4DIR||require('path').join(__dirname,'fixtures_v4/');
const spec={species:[{id:"Tm",levels:["3H6","3F4","3H4"]}],procs:[
 {type:"abs",sp:"Tm",from:"3H6",to:"3F4",sigma:"sG"},{type:"abs",sp:"Tm",from:"3F4",to:"3H4",sigma:"sE"},
 {type:"dec",sp:"Tm",from:"3F4",to:"3H6",k:"W2"},{type:"dec",sp:"Tm",from:"3H4",to:"3H6",k:"W3*(1-b)"},{type:"dec",sp:"Tm",from:"3H4",to:"3F4",k:"W3*b"},
 {type:"et",d:["Tm","3H4","3F4"],a:["Tm","3H6","3F4"],W:"s31",scale:"none"},{type:"et",d:["Tm","3F4","3H6"],a:["Tm","3F4","3F4"],W:"q22",scale:"none"},{type:"et",d:["Tm","3F4","3H6"],a:["Tm","3F4","3H4"],W:"q23",scale:"none"}],
 channels:[{id:"800",terms:[["Tm","3H4",1]]}]};
const names=["sG","sE","W3","b","s31","q22","q23"];
const mkParams=th=>{ const p={W2:{v:1000,lo:1,hi:1e6,fixed:true}}; names.forEach((n,j)=>p[n]={v:Math.pow(10,th[j]),lo:Math.pow(10,th[j]-3),hi:Math.min(n==="b"?0.999:1e30,Math.pow(10,th[j]+3))}); return p; };
const fmtI=v=>(v/K).toFixed(2)+" mW";
function report(tag,A){ console.log(`\n[${tag}] ${A.worst.toUpperCase()} — ${A.verdict}`); A.items.forEach(q=>console.log(`  ${q.key} ${q.level.padEnd(4)} ${q.title}: ${q.summary}`)); }
const rnd=D.mulberry(42); const g=()=>{ let u=0; while(!u) u=rnd(); return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*rnd()); };
let pass=0,fail=0; const expect=(c,m)=>{ if(c){pass++;} else {fail++; console.log('  ✗ 기대 불일치: '+m);} };

/* ---- T1: 모델 자신이 만든 데이터 ---- */
const bi=JSON.parse(fs.readFileSync(V4+'bi_h12.json')); const th0=bi.th.slice(); th0[0]+=Math.log10(bi.f.Pdn/6.96); th0[1]+=Math.log10(bi.f.Pdn/6.96);
function synth(th,mod,rel,seed){
  const M=D.compile(spec); const env=Object.assign({W2:1000},Object.fromEntries(names.map((n,j)=>[n,Math.pow(10,th[j])])));
  const Pu=D.logGrid(4.4,13.4,210).slice(0,100), Pd=Pu.slice().reverse();
  const yu=D.sweep(M,env,Pu.map(p=>p*K),1,1064).n.map(n=>n[2]), yd=D.sweep(M,env,Pd.map(p=>p*K),-1,1064).n.map(n=>n[2]);
  const S=1.3e6/Math.max(...yd); const mk=(P,y)=>y.map((v,i)=>S*v*(mod?mod(P[i]):1)*(1+rel*g())+400*g());
  return {M,env,data:[{type:"scurve",ch:"800",dir:1,I:Pu.map(p=>p*K),y:mk(Pu,yu),sbg:400,k:3,relErr:rel},{type:"scurve",ch:"800",dir:-1,I:Pd.map(p=>p*K),y:mk(Pd,yd),sbg:400,k:3,relErr:rel}]};
}
function traceOf(M,env,P,tEnd,npt,noise){ const T=D.logGrid(tEnd/2000,tEnd,npt); const r=D.integrate(M,env,P*K,1064,tEnd,null,{tOut:T,maxSteps:300000}); const y=r.outN.map(n=>n[2]); const f=y[y.length-1]; return {type:"trace",ch:"800",I:[P*K],t:r.outT,y:y.map(v=>1e3*v/f*(1+noise*g()))}; }
function runCase(tag,data,th,fit){
  const prob={spec,nm:1064,beam:"peak",params:mkParams(th),samples:[{id:"S",vars:{},data}]}; const B=D.buildProblem(prob);
  let theta=D.thetaFromParams(B); if(fit){ const r=D.lm(B,theta,{maxit:60}); theta=r.theta; }
  const A=D.adequacy(B,theta,{fmtI}); report(tag,A); return A; }
{ const s=synth(th0,null,0.03); const P0=8.35*1.05; s.data.push(traceOf(s.M,s.env,P0,0.6,600,0.02));
  const A=runCase("T1 참값 = 데이터 생성 모델",s.data,th0,true); expect(A.worst!=="fail","T1 은 fail 이 없어야"); }
/* 오탐률: 잡음 시드 8개 (트레이스 없이) */
{ let nf=0,nw=0; for(let k=0;k<8;k++){ const s=synth(th0,null,0.05); const prob={spec,nm:1064,beam:"peak",params:mkParams(th0),samples:[{id:"S",vars:{},data:s.data}]}; const B=D.buildProblem(prob); const r=D.lm(B,D.thetaFromParams(B),{maxit:40}); const A=D.adequacy(B,r.theta,{fmtI}); if(A.worst==="fail") nf++; if(A.worst==="warn") nw++; console.log('   rep',k,A.items.filter(q=>q.level!=='ok'&&q.level!=='na').map(q=>q.key+':'+q.level+' '+q.summary).join(' | ')); }
  console.log(`\n[T1-반복] 잡음 5 %, 8회: fail ${nf}, warn ${nw}`); expect(nf===0,"오탐 fail 0"); }
/* ---- T3: 구조가 다른 데이터 (위 가지 기울기 +0.35) ---- */
{ const s=synth(th0,P=>P>6.9?Math.pow(P/6.9,0.35):1,0.03); const A=runCase("T3 위 가지 모양이 다른 데이터",s.data,th0,true);
  expect(A.items.find(q=>q.key==="A3").level==="fail"||A.items.find(q=>q.key==="A2").level==="fail","T3 은 A2/A3 fail"); }
/* ---- T3b: 과도 상승이 있는 트레이스 (느린 음되먹임 모사) ---- */
{ const s=synth(th0,null,0.03); const tr=traceOf(s.M,s.env,8.35*1.05,0.6,600,0.01); tr.y=tr.y.map((v,i)=>v*(1+0.12*Math.exp(-tr.t[i]/0.08)*(1-Math.exp(-tr.t[i]/0.01))));
  const f=tr.y.slice(-120).reduce((a,b)=>a+b,0)/120; s.data.push(tr); const A=runCase("T3b 과도 상승 트레이스",s.data,th0,false); expect(A.items.find(q=>q.key==="A6").level==="fail","T3b 는 A6 fail"); }
/* ---- T2: WJ-0128c 실측 ---- */
{ const S=JSON.parse(fs.readFileSync(V4+'scans.json')).nir; const bg=p=>1516.0*p-521.9;
  const up={type:"scurve",ch:"800",dir:1,I:[],y:[],sbg:450,k:3,relErr:0.05}, dn={type:"scurve",ch:"800",dir:-1,I:[],y:[],sbg:450,k:3,relErr:0.05};
  S.P.forEach((p,i)=>{ const o=i<=S.turn?up:dn; o.I.push(p*K); o.y.push(S.R[i]-bg(p)); });
  const rt=JSON.parse(fs.readFileSync(V4+'rise_trace.json'));
  const best=JSON.parse(fs.readFileSync(V4+'refine.json'))[0].th;
  const A=runCase("T2 WJ-0128c NIR + 상승 트레이스 (3준위 최적)",[up,dn,{type:"trace",ch:"800",I:[5.699*K],t:rt.t,y:rt.y}],best,false);
  expect(A.worst==="fail","T2 는 fail"); A.items.forEach(q=>console.log(`   ${q.key}: ${q.detail||""}`)); }
console.log(`\n적합성 시험 ${pass} 통과 / ${fail} 실패`); process.exit(fail?1:0);
