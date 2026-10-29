const fs=require('fs'), path=require('path');
const D='/home/claude/w/src/modules/xps';
global.window=global; 
const C=require(D+'/src/xps_core.js'); const X=C.XPS_CORE||C||global.XPS_CORE;
const tm=fs.readFileSync(D+'/test/data/casaval/Shape_GLP30.vms','utf8').split('\n'); const k=tm.indexOf('762');
const head=tm.slice(0,k+5), tail=['end of experiment',''];
const E=Array.from({length:381},(_,i)=>279.08+0.05*i);
const cases={ K_GL30n100_F1_2:'GL(30,100)', K_SGL30n100_F1_2:'SGL(30,100)', K_TLA12_3_0_F1:'TLA(1.2,3,0)', K_TLA12_3_100_F1:'TLA(1.2,3,100)',
 K_DS01_100_F1:'DS(0.1,100)', K_ST17_LA11200_F1:'ST(1.7)LA(1,1,200)', K_LA2_3_0_F1:'LA(2,3,0)', K_LA1_1_250_F1:'LA(1,1,250)'};
const fw={K_GL30n100_F1_2:1.2,K_SGL30n100_F1_2:1.2};
const truth={};
for(const [name,s] of Object.entries(cases)){
  const r=X.casaToComp(s,{be:285,area:15000,fwhm:fw[name]||1.0,strict:true}); const c=r.comp||r;
  const y=X.compEval(c,Float64Array.from(E)); 
  const v=Array.from(y,t=>1000+t);
  const body=[]; for(const t of v.slice().reverse()) body.push(t.toFixed(6),'1');
  fs.writeFileSync('/tmp/k3/'+name+'.vms',head.concat(body,tail).join('\n'));
  let ar=0; for(let i=1;i<E.length;i++) ar+=0.5*(y[i]+y[i-1])*0.05;
  truth[name]={casa:s,comp:JSON.parse(JSON.stringify(c)),areaInWindow:ar,max:Math.max(...y),note:r.note};
}
fs.writeFileSync('/tmp/k3/K_truth.json',JSON.stringify(truth,null,1)); console.log(Object.entries(truth).map(([k,v])=>k+' '+v.areaInWindow.toFixed(1)+' max '+v.max.toFixed(1)).join('\n'));
