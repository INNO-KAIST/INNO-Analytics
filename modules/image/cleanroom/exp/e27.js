const O=require("./oracle.js"), CR=require("../askel_cr.js");
const rows=[
"..............",
"...########...",
"..#........#..",
"..#........#..",
"...########...",
"...#..........",
"...#..........",
".............."];
const v=O.fromRows(rows); const g=CR._t.debugGraph(v,{pw:1,ph:1,pd:1}); const L=g.edges.find(e=>e.slabs.length>5).slabs; const at=(k)=>L[k][0]+v.w*L[k][1];
const cases=[]; const add=(name,f)=>{ const o=new Uint8Array(v.w*v.h).fill(200); f(o); cases.push(Object.assign({},v,{orig:o,prune:2,lsp:false,name})); };
add("s5=10, s4=250",o=>{ o[at(5)]=10; o[at(4)]=250; });
add("s5=10, s6=250",o=>{ o[at(5)]=10; o[at(6)]=250; });
add("s5=10, s4=250, s6=250",o=>{ o[at(5)]=10; o[at(4)]=250; o[at(6)]=250; });
add("s5=10 s6=10",o=>{ o[at(5)]=10; o[at(6)]=10; });
add("bg above s5 = 0",o=>{ o[at(5)-v.w]=0; });
add("bg above s5 = 0, s4 = 255",o=>{ o[at(5)-v.w]=0; o[at(4)]=255; });
cases.push(Object.assign({},v,{prune:0,lsp:false,name:"base"}));
const os=O.run(cases); const base=os[os.length-1];
cases.slice(0,-1).forEach((c,i)=>{ const rem=[]; for(const [p] of base.tag) if(!os[i].tag.has(p)) rem.push(p%c.w+","+((p/c.w)|0)); console.log(c.name.padEnd(28),"removed",rem.join(" "),"idx",rem.map(s=>L.findIndex(q=>q[0]+","+q[1]===s)).join(" ")); });
console.log(L.map(p=>p[0]+","+p[1]).join(" "));
