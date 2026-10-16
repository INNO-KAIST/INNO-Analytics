// rare structures: end vertex with a self-loop; end-point adjacent to a visited junction (zero-length-walk edge); compare with prune ends
const gen=require("./gen.js"), C=require("./cmp.js"), CR=require("../askel_cr.js");
const selfEnd=[], zeroWalk=[];
for(let seed=500;seed<540;seed++) for(const c of gen(seed,200)){ const g=CR._t.debugGraph(c,{pw:c.cal[0],ph:c.cal[1],pd:c.cal[2]});
  if(g.edges.some(e=>e.v1===e.v2&&g.trees[e.tree].vCl[e.v1]<0&&!g.trees[e.tree].cycle)) selfEnd.push(c);
  { const m=new Map(); for(const e of g.edges){ if(e.slabs.length===0&&g.trees[e.tree].vCl[e.v1]<0&&g.trees[e.tree].vCl[e.v2]>=0){ const k=e.tree+':'+e.v2; m.set(k,(m.get(k)||0)+1); } } if([...m.values()].some(x=>x>1)) zeroWalk.push(c); } }
console.log("end self-loops",selfEnd.length,"zero-walk edges",zeroWalk.length);
const cases=[]; for(const c of selfEnd.slice(0,150).concat(zeroWalk.slice(0,150))) for(const [P,pe] of [[0,true],[1,true],[2,false],[1,false]]) cases.push(Object.assign({},c,{prune:P,pruneEnds:pe}));
console.log("mismatching",C.cmp(cases).nbad,"of",cases.length);
