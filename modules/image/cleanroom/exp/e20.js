const H=require("./hyp.js");
const base={dfs:"stackE",method:1,metric:"slabs",order:"fwd",le:false,allRoots:true,cycleStart:true};
const V={first:Object.assign({},base,{root:"first"}),
  firstEnd:Object.assign({},base,{root:(T,E)=>{ const i=T.vCl.findIndex(c=>c<0); return i<0?0:i; }}),
  firstEdgeV1:Object.assign({},base,{root:(T,E)=>E.length?E[0].v1:0}),
  len:Object.assign({},base,{root:"first",metric:"len"}),
};
for(const seed of (process.env.SEEDS||"11,12").split(",")) console.log(seed,H.evaluate(1,V,+(process.env.N||1000),+seed,+(process.env.SHOW||0)).res);
