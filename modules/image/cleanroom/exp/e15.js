const O=require("./oracle.js"), CR=require("../askel_cr.js");
module.exports=function removedIdx(shapes,opts){ const cases=[], base=[]; for(const [k,rows] of Object.entries(shapes)){ const v=O.fromRows(rows); base.push(Object.assign({},v,{lsp:true,name:k})); cases.push(Object.assign({},v,{lsp:true,name:k},opts)); }
  const os=O.run(base.concat(cases)); const nb=base.length;
  cases.forEach((c,i)=>{ const o0=os[i], o1=os[nb+i]; const removed=[]; for(const [p] of o0.tag) if(!o1.tag.has(p)) removed.push([p%c.w,(p/c.w)|0,0].join());
    console.log("==",c.name,JSON.stringify(opts),"removed",removed.join(" ")); console.log(O.ascii(c.w,c.h,c.data,new Map(removed.map(s=>{const [x,y]=s.split(",").map(Number); return [x+c.w*y,"X"];}))));
    const g=CR._t.debugGraph(base[i],{pw:1,ph:1,pd:1});
    g.edges.forEach((e,k)=>{ const idx=e.slabs.map((p,j)=>removed.includes(p.join())?j:-1).filter(j=>j>=0); console.log("  edge",k,"len",e.len.toFixed(3),"v",e.v1,e.v2,"("+e.p1+")-("+e.p2+")","nslab",e.slabs.length,"removed idx",JSON.stringify(idx)); });
    console.log(o1.raw.split("\n").filter(l=>l.startsWith("T ")).join("\n")); });
  return os; };
if(require.main===module) module.exports({
lenVsSlab:[
"............",
"...#####....",
"..#.....#...",
".#.......#..",
"..#.....#...",
"...#...#....",
"....#.#.....",
".....#......",
"............"],
lenVsSlab2:[ // two clusters joined by a diagonal path (4 slabs, long) and straight path (5 slabs)
"..............",
"..#.......#...",
"...#######....",
"...#.....#....",
"....#...#.....",
".....#.#......",
"......#.......",
".............."],
},{prune:1});
