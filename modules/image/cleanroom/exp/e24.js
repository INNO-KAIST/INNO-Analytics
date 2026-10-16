const O=require("./oracle.js"), CR=require("../askel_cr.js");
const rows=[
"............",
"...######...",
"..#......#..",
"..#......#..",
"...######...",
"...#........",
"...#........",
"............"];
const v=O.fromRows(rows); const g=CR._t.debugGraph(v,{pw:1,ph:1,pd:1}); const L=g.edges[1].slabs; console.log("vertices",JSON.stringify(g.trees[0].vVox));
const cases=[]; const add=(name,f)=>{ const o=new Uint8Array(v.w*v.h); for(let i=0;i<o.length;i++) o[i]=f(i%v.w,(i/v.w)|0); for(const P of [2,3]) cases.push(Object.assign({},v,{orig:o,prune:P,lsp:false,name:name+" P"+P})); };
const idxOf=(x,y)=>L.findIndex(q=>q[0]===x&&q[1]===y);
add("decreasing along list",(x,y)=>{ const k=idxOf(x,y); return k>=0?200-k:250; });
add("decreasing, bg 0",(x,y)=>{ const k=idxOf(x,y); return k>=0?200-k:0; });
add("increasing along list",(x,y)=>{ const k=idxOf(x,y); return k>=0?100+k:250; });
add("all 100, junction voxels 5",(x,y)=>{ const k=idxOf(x,y); if(k>=0) return 100; return 5; });
cases.push(Object.assign({},v,{prune:0,lsp:false,name:"base"}));
const os=O.run(cases); const base=os[os.length-1];
cases.slice(0,-1).forEach((c,i)=>{ const rem=[]; for(const [p] of base.tag) if(!os[i].tag.has(p)) rem.push(p%c.w+","+((p/c.w)|0)); console.log(c.name.padEnd(30),"removed",rem.join(" "),"idx",rem.map(s=>L.findIndex(q=>q[0]+","+q[1]===s)).join(" "),"nslab",L.length); });
