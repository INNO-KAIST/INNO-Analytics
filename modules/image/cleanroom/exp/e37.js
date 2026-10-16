const H=require("./hyp.js");
const root=(T,E)=>{ const i=T.vCl.findIndex(c=>c<0); return i<0?0:i; };
const slabs=(e)=>e.slabs; const withV=(e,T)=>T.vVox[e.v1].concat(e.slabs,e.v2!==e.v1?T.vVox[e.v2]:[]); const withF=(e,T)=>[T.vFirst[e.v1]].concat(e.slabs,[T.vFirst[e.v2]]);
const b={dfs:"stackE",order:"fwd",le:false,allRoots:true,root,method:3,vox:slabs,vox2:slabs,nb:"cal",oob:"zero",cycleStartAll:true};
const V={win:b, raw:Object.assign({},b,{avgRaw:true}), winV:Object.assign({},b,{vox:withV}), rawV:Object.assign({},b,{vox:withV,avgRaw:true}), winF:Object.assign({},b,{vox:withF}), rawF:Object.assign({},b,{vox:withF,avgRaw:true}), sum:Object.assign({},b,{sumOnly3:true})};
console.log(H.evaluate(3,V,+(process.env.N||300),+(process.env.SEED||21),0).res);
