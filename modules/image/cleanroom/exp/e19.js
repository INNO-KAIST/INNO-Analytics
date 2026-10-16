const gen=require("./gen.js"), CR=require("../askel_cr.js"), O=require("./oracle.js");
const [seed,N,idx,x0,x1,y0,y1]=process.argv.slice(2).map(Number); const c=gen(seed,N)[idx];
const g=CR._t.debugGraph(c,{pw:c.cal[0],ph:c.cal[1],pd:c.cal[2]});
const o=O.run([Object.assign({},c,{prune:0}),Object.assign({},c,{prune:+(process.env.P||1)})]);
const rem=new Map(); for(const [p] of o[0].tag) if(!o[1].tag.has(p)) rem.set(p,"X");
const tagc=new Map([...o[0].tag].map(([p,v])=>[p,rem.has(p)?"X":v===30?"E":v===70?"J":"s"]));
let s=""; for(let y=y0;y<=y1;y++){ let l=String(y).padStart(3)+" "; for(let x=x0;x<=x1;x++){ const i=x+c.w*y; l+=tagc.get(i)||"."; } s+=l+"\n"; } console.log("cols from",x0); console.log(s);
g.edges.forEach((e,k)=>{ const inwin=e.slabs.concat([e.p1,e.p2]).some(p=>p[0]>=x0&&p[0]<=x1&&p[1]>=y0&&p[1]<=y1); if(inwin) console.log("edge",k,"tree",e.tree,"v",e.v1,e.v2,"("+e.p1+")-("+e.p2+")","len",e.len.toFixed(3),"slabs",e.slabs.map(p=>p.join(",")).join(" ")); });
