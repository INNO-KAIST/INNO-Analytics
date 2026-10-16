const D=require('../src/dre_core.js');
// Lee 템플릿 합성 데이터로 insight 산출 점검
const T=D.templateLee(); const M=D.compile(T.spec); const env=Object.assign({},Object.fromEntries(Object.entries(T.params).map(([k,p])=>[k,p.v])),{c:8});
const I=D.logGrid(1e3,2e5,40); const y=D.observe(M,env,I,1,1064,0,{}).y.map(v=>v*1e6);
const rnd=D.mulberry(3); const g=()=>{ let u=0; while(!u) u=rnd(); return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*rnd()); };
const prob={spec:T.spec,nm:1064,beam:"peak",params:JSON.parse(JSON.stringify(T.params)),samples:[{id:"S1",vars:{c:8},data:[{type:"scurve",ch:"800 nm (³H₄)",dir:1,I,y:y.map(v=>v*(1+0.04*g())+30*g()),sbg:30,k:3,relErr:0.04},{type:"decay",ch:"800 nm (³H₄)",I:[1e3],t:[D.decayTime(M,env,1e3,1064,0).t],relErr:0.05}]}]};
prob.params.w2nr.fixed=true; prob.params.w3nr.fixed=true; prob.params.w2nr.per=undefined; prob.params.w3nr.per=undefined;
const B=D.buildProblem(prob); const ms=D.multistart(B,{starts:12,seed:2,maxit:60}); const best=ms[0];
const t0=Date.now(); const X=D.insight(B,best,ms); console.log('insight ms',Date.now()-t0);
console.log('free',B.free.map(f=>f.key).join(' '),'\nse(log10)',X.se.map(v=>v.toFixed(3)).join(' '));
console.log('corr row0',X.corr[0].map(v=>v.toFixed(2)).join(' '));
console.log('clusters',X.clusters.map(c=>`${c.chi2.toFixed(1)}×${c.count}`).join(', '));
console.log('flux rows',X.flux[0].rows.length,'procs',X.procs.map(p=>p.name).join(' | '));
X.flux[0].rows.forEach(r=>console.log((r.I).toExponential(2),r.branch,r.share.map(v=>(100*v).toFixed(1)).join(' ')));
const A=D.adequacy(B,best.theta); console.log('adequacy',A.worst,A.items.map(q=>q.key+':'+q.level).join(' '));
const t1=Date.now(); const sl=D.slice2D(B,best.theta,0,1,{n:9}); console.log('slice ms',Date.now()-t1, 'min',Math.min(...sl.z.flat()).toFixed(2),'best',best.chi2.toFixed(2));
