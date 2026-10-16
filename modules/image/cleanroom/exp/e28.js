const H=require("./hyp.js");
const root=(T,E)=>{ const i=T.vCl.findIndex(c=>c<0); return i<0?0:i; };
const V={}; for(const oob of ["zero","skip"]) for(const order of ["fwd","rev"]) for(const le of [false,true]) V[`${oob}-${order}-${le?"le":"lt"}`]={dfs:"stackE",order,le,allRoots:true,root,method:2,vox:(e)=>e.slabs,nb:"3x3x3",oob};
console.log(H.evaluate(2,V,+(process.env.N||300),+(process.env.SEED||21),0).res);
