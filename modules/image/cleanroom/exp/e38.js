const H=require("./hyp.js");
const root=(T,E)=>{ const i=T.vCl.findIndex(c=>c<0); return i<0?0:i; };
const slabs=(e)=>e.slabs; const withV=(e,T)=>T.vVox[e.v1].concat(e.slabs,e.v2!==e.v1?T.vVox[e.v2]:[]); const withV2=(e,T)=>T.vVox[e.v1].concat(e.slabs,T.vVox[e.v2]);
const b={dfs:"stackE",order:"fwd",le:false,allRoots:true,root,vox:slabs,vox2:slabs,nb:"cal",oob:"zero",cycleStartAll:true,metric:"slabs"};
const N=+(process.env.N||1000), seed=+(process.env.SEED||77);
console.log("P1",H.evaluate(1,{a:Object.assign({},b,{method:1}),noAll:Object.assign({},b,{method:1,allRoots:false})},N,seed,0).res);
console.log("P2",H.evaluate(2,{a:Object.assign({},b,{method:2}),noAll:Object.assign({},b,{method:2,allRoots:false})},N,seed,0).res);
console.log("P3",H.evaluate(3,{a:Object.assign({},b,{method:3,vox:withV}),twice:Object.assign({},b,{method:3,vox:withV2})},N,seed,0).res);
