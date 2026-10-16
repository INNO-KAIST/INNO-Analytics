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
const v=O.fromRows(rows); const g=CR._t.debugGraph(v,{pw:1,ph:1,pd:1}); const L=g.edges[1].slabs;
const cases=[]; const lvl=+(process.env.LVL||10);
for(let k=0;k<L.length;k++){ const o=new Uint8Array(v.w*v.h).fill(200); o[L[k][0]+v.w*L[k][1]]=lvl; cases.push(Object.assign({},v,{orig:o,prune:+(process.env.P||2),lsp:false,k})); }
cases.push(Object.assign({},v,{prune:0,lsp:false}));
const os=O.run(cases); const base=os[os.length-1];
cases.slice(0,-1).forEach((c,i)=>{ const rem=[]; for(const [p] of base.tag) if(!os[i].tag.has(p)) rem.push(p%c.w+","+((p/c.w)|0)); const ri=rem.map(s=>L.findIndex(q=>q[0]+","+q[1]===s)); console.log("dark idx",c.k,"("+L[c.k]+")","removed",rem.join(" "),"idx",ri.join(" ")); });
