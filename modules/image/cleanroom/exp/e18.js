const H=require("./hyp.js");
const V={"slabs-rev-lt":{method:1,metric:"slabs",order:"rev",le:false,root:"first",allRoots:true,cycleStart:true}};
const r=H.evaluate(1,V,+(process.env.N||400),+(process.env.SEED||3),20); console.log(r.res);
