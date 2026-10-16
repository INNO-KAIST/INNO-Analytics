const H=require("./hyp.js");
const V={}; for(const dfs of ["rec","stackE"]) for(const metric of ["slabs"]) for(const order of ["fwd","rev"]) for(const le of [false,true]) V[`${dfs}-${metric}-${order}-${le?"le":"lt"}`]={dfs,method:1,metric,order,le,root:"first",allRoots:true,cycleStart:true};
console.log(H.evaluate(1,V,+(process.env.N||400),+(process.env.SEED||3),+(process.env.SHOW||0)).res);
