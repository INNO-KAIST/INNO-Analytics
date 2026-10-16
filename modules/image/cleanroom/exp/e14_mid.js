// which voxel is removed by "shortest branch" for self-loops of various slab counts, rings, parallel edges
const O=require("./oracle.js"), C=require("./cmp.js");
const cases=[];
function circle(w,h,cx,cy,r){ const v={w,h,d:1,data:new Uint8Array(w*h)}; for(let a=0;a<360;a+=0.5){ const x=Math.round(cx+r*Math.cos(a*Math.PI/180)), y=Math.round(cy+r*Math.sin(a*Math.PI/180)); v.data[x+w*y]=255; } return v; }
// rings (cycle components)
for(const [cx,cy,r] of [[10,10,4],[10,10,5],[11,9,6],[12,12,7],[10,11,3]]){ const v=circle(25,25,cx,cy,r); cases.push(Object.assign(v,{prune:1,lsp:true,name:`ring ${cx},${cy},${r}`})); }
// rectangle loop with a tail (self-loop on a cluster), varying sizes
for(const [lw,lh] of [[3,2],[4,2],[4,3],[5,3],[6,4]]){ const w=lw+12,h=lh+8, v={w,h,d:1,data:new Uint8Array(w*h)}; const x0=3,y0=3,x1=x0+lw+1,y1=y0+lh+1;
  for(let x=x0+1;x<x1;x++){ v.data[x+w*y0]=255; v.data[x+w*y1]=255; } for(let y=y0+1;y<y1;y++){ v.data[x0+w*y]=255; v.data[x1+w*y]=255; }
  const ym=y0+1; for(let x=x1+1;x<x1+5;x++) v.data[x+w*ym]=255; cases.push(Object.assign(v,{prune:1,lsp:true,name:`rect ${lw}x${lh}+tail`})); }
const base=cases.map(c=>Object.assign({},c,{prune:0}));
const os=O.run(base.concat(cases)); const nb=base.length;
cases.forEach((c,i)=>{ const o0=os[i], o1=os[nb+i]; const removed=[]; for(const [p] of o0.tag) if(!o1.tag.has(p)) removed.push([p%c.w,(p/c.w)|0]);
  const r=C.mine(base[i]); // my un-pruned analysis: list edges with slab lists
  console.log(c.name,"removed",JSON.stringify(removed),"edges",o0.edges.map(e=>e.slice(1).join(" ")).join(" | ")); });
const CR=require("../askel_cr.js");
cases.forEach((c,i)=>{ const o0=os[i], o1=os[nb+i]; const removed=[]; for(const [p] of o0.tag) if(!o1.tag.has(p)) removed.push([p%c.w,(p/c.w)|0,0].join());
  const g=CR._t.debugGraph(base[i],{pw:1,ph:1,pd:1});
  g.edges.forEach((e,k)=>{ const idx=e.slabs.map((p,j)=>removed.includes(p.join())?j:-1).filter(j=>j>=0); console.log("  ",c.name,"edge",k,"len",e.len.toFixed(3),"v",e.v1,e.v2,"nslab",e.slabs.length,"removed idx",JSON.stringify(idx)); }); });
