// kernel probing: rectangle loop (self-loop on a junction with a tail), uniform 200, one dark voxel at an offset from slab (10,4)
const O=require("./oracle.js"), CR=require("../askel_cr.js");
const w=24,h=16, v={w,h,d:1,data:new Uint8Array(w*h)}; const put=(x,y)=>v.data[x+w*y]=255;
for(let x=4;x<=19;x++){ put(x,4); put(x,11); } for(let y=5;y<=10;y++){ put(3,y); put(20,y); } put(3,4); // corner makes junction? keep simple
for(let y=12;y<=14;y++) put(3,y); // tail from corner (3,11)
const g=CR._t.debugGraph(v,{pw:1,ph:1,pd:1}); const L=g.edges.find(e=>e.slabs.length>10).slabs; const idx=(s)=>L.findIndex(q=>q[0]+","+q[1]===s);
const cases=[]; const offs=[[0,0],[0,1],[0,-1],[1,1],[0,2],[0,-2],[1,2],[2,2],[0,3]];
for(const [dx,dy] of offs){ for(const val of [0,150]){ const o=new Uint8Array(w*h).fill(200); o[(10+dx)+w*(4+dy)]=val; cases.push(Object.assign({},v,{orig:o,prune:2,lsp:false,name:`dark ${val} at (10+${dx},4+${dy})`})); } }
cases.push(Object.assign({},v,{prune:0,lsp:false,name:"base"}));
const os=O.run(cases); const base=os[os.length-1];
cases.slice(0,-1).forEach((c,i)=>{ const rem=[]; for(const [p] of base.tag) if(!os[i].tag.has(p)) rem.push(p%w+","+((p/w)|0)); console.log(c.name.padEnd(26),"removed",rem.join(" "),"idx",rem.map(idx).join(" ")); });
console.log("slab list",L.map(p=>p[0]+","+p[1]).join(" "));
