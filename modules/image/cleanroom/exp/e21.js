const H=require("./hyp.js");
const root=(T,E)=>{ const i=T.vCl.findIndex(c=>c<0); return i<0?0:i; };
const base={dfs:"stackE",order:"fwd",le:false,allRoots:true,root};
const slabs=(e)=>e.slabs;
const V2={}, V3={};
for(const nb of ["3x3","3x3x3"]) for(const oob of ["skip","zero"]){ V2[nb+"-"+oob]=Object.assign({},base,{method:2,vox:slabs,nb,oob}); V3[nb+"-"+oob]=Object.assign({},base,{method:3,vox:slabs,vox2:slabs,nb,oob}); }
const N=+(process.env.N||600), seed=+(process.env.SEED||21);
console.log("prune2",H.evaluate(2,V2,N,seed,+(process.env.SHOW||0)).res);
console.log("prune3",H.evaluate(3,V3,N,seed,+(process.env.SHOW||0)).res);
