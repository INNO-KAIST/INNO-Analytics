// lowest intensity voxel on a single loop with a controlled intensity image
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
const v=O.fromRows(rows); const cases=[];
const mk=(f)=>{ const o=new Uint8Array(v.w*v.h); for(let i=0;i<o.length;i++) o[i]=f(i%v.w,(i/v.w)|0); return o; };
const variants={
  "dark at (7,1)":mk((x,y)=>x===7&&y===1?10:200),
  "dark at (2,2) [slab]":mk((x,y)=>x===2&&y===2?10:200),
  "dark at (3,4) [junction]":mk((x,y)=>x===3&&y===4?10:200),
  "dark at (3,6) [tail end]":mk((x,y)=>x===3&&y===6?10:200),
  "dark background only":mk((x,y)=>v.data[x+v.w*y]?200:5),
  "ramp x":mk((x,y)=>50+10*x),
  "ramp -x":mk((x,y)=>200-10*x),
};
for(const [k,o] of Object.entries(variants)) for(const P of [2,3]) cases.push(Object.assign({},v,{orig:o,prune:P,lsp:false,name:k+" P"+P}));
cases.push(Object.assign({},v,{prune:0,lsp:false,name:"base"}));
const os=O.run(cases); const base=os[os.length-1];
cases.forEach((c,i)=>{ if(c.name==="base") return; const rem=[]; for(const [p] of base.tag) if(!os[i].tag.has(p)) rem.push([p%c.w,(p/c.w)|0]); console.log(c.name.padEnd(32),"removed",JSON.stringify(rem)); });
const g=CR._t.debugGraph(v,{pw:1,ph:1,pd:1}); g.edges.forEach((e,k)=>console.log("edge",k,e.v1,e.v2,e.slabs.map(p=>p[0]+","+p[1]).join(" ")));
