const gen=require("./gen.js"), CR=require("../askel_cr.js"), O=require("./oracle.js");
const [seed,N,idx,x0,y0,x1,y1,z,tr]=process.argv.slice(2).map(Number); const c=gen(seed,N)[idx]; const P=+(process.env.P||2);
const o=O.run([Object.assign({},c,{prune:0}),Object.assign({},c,{prune:P})]);
const zz=z||0; for(let y=y0;y<=y1;y++){ let l=String(y).padStart(3)+": "; for(let x=x0;x<=x1;x++){ const i=x+c.w*(y+c.h*zz); const t=o[0].tag.get(i); const r=t&&!o[1].tag.has(i); l+=(r?"X":t===30?"E":t===70?"J":t===127?"s":".")+String(c.orig[i]).padStart(3)+" "; } console.log(l); }
const g=CR._t.debugGraph(c,{pw:c.cal[0],ph:c.cal[1],pd:c.cal[2]}); console.log("tree",tr,"vertices",JSON.stringify(g.trees[tr].vVox));
g.edges.forEach((e,k)=>{ if(e.tree===tr) console.log("edge",k,"v",e.v1,e.v2,"len",e.len.toFixed(3),"slabs",e.slabs.map(p=>p.join(",")).join(" ")); });
