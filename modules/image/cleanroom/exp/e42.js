// find volumes where the slab counter differs from the number of slab voxels, then compare both end-pruning variants
const gen=require("./gen.js"), C=require("./cmp.js"), CR=require("../askel_cr.js");
const found=[]; for(let seed=400;seed<460&&found.length<60;seed++){ for(const c of gen(seed,200)){ const r=CR.analyzeSkeleton(c,{pw:c.cal[0],ph:c.cal[1],pd:c.cal[2]},{});
    let nslab=new Map(); // slab voxels per tree: count tag 127 per tree via edges is hard; use total
    const tot=r.trees.reduce((s,t)=>s+t.slabs,0); let cnt=0; for(const v of r.tag) if(v===127) cnt++; if(tot!==cnt) found.push(Object.assign({},c,{pruneEnds:true,prune:0,name:c.name+" s"+seed+" counter "+tot+" voxels "+cnt})); } }
console.log("found",found.length); found.slice(0,5).forEach(c=>console.log(c.name));
console.log("current:",C.cmp(found,false).nbad);
