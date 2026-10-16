const gen=require("./gen.js"), CR=require("../askel_cr.js");
const [seed,N,idx,...pts]=process.argv.slice(2); const c=gen(+seed,+N)[+idx];
const g=CR._t.debugGraph(c,{pw:c.cal[0],ph:c.cal[1],pd:c.cal[2]});
const want=new Set(pts); const trees=new Set(); g.edges.forEach(e=>{ if(e.slabs.some(p=>want.has(p.join()))) trees.add(e.tree); });
for(const t of trees){ console.log("tree",t,"vertices",JSON.stringify(g.trees[t].vVox)); g.edges.forEach((e,k)=>{ if(e.tree===t) console.log(" edge",k,"v",e.v1,e.v2,"slabs",e.slabs.map(p=>p.join()).join(" ")); }); }
