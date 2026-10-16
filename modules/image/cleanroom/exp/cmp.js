// compare clean-room output with oracle for given cases; prints diffs
const O=require("./oracle.js"); const CR=require(process.env.CRMOD?require("path").resolve(process.env.CRMOD):"../askel_cr.js");
const PR=["none","shortest branch","lowest intensity voxel","lowest intensity branch"];
const close=(a,b)=>a===b||Math.abs(a-b)<=1e-12*Math.max(1,Math.abs(b));
function mine(c){ return CR.analyzeSkeleton({w:c.w,h:c.h,d:c.d||1,data:c.data},{pw:(c.cal||[1,1,1])[0],ph:(c.cal||[1,1,1])[1],pd:(c.cal||[1,1,1])[2]},{lsp:c.lsp!==false,prune:PR[c.prune||0],pruneEnds:!!c.pruneEnds,orig:c.orig?{data:c.orig}:null}); }
function diff(c,o,r){ const ex=[]; const lsp=c.lsp!==false;
  if(r.trees.length!==o.trees.length) ex.push(`trees ${r.trees.length} vs ${o.trees.length}`);
  else r.trees.forEach((t,k)=>{ const a=[t.branches,t.junctions,t.endPoints,t.junctionVoxels,t.slabs,t.avgLength,t.triples,t.quadruples,t.maxLength].concat(lsp?[t.lsp,t.spx,t.spy,t.spz]:[]); a.forEach((v,j)=>{ if(!close(v,o.trees[k][j])) ex.push(`tree ${k} f${j} ${v} vs ${o.trees[k][j]}`); }); });
  if(r.edges.length!==o.edges.length) ex.push(`edges ${r.edges.length} vs ${o.edges.length}`);
  else r.edges.forEach((e,k)=>{ const a=[e.tree,e.length,...e.v1,...e.v2,e.lengthRa,e.nSlabs]; a.forEach((v,j)=>{ if(!close(v,o.edges[k][j])) ex.push(`edge ${k} f${j} ${v} vs ${o.edges[k][j]}`); }); });
  let bad=0; for(let p=0;p<r.tag.length;p++) if(r.tag[p]!==(o.tag.get(p)||0)) bad++; if(bad) ex.push(`tag ${bad}`);
  if(lsp){ let b2=0; for(let p=0;p<r.tag.length;p++) if(((r.lspImage&&r.lspImage[p]===96)?1:0)!==(o.lsp.has(p)?1:0)) b2++; if(b2) ex.push(`lspimg ${b2}`); }
  return ex; }
function cmp(cases,verbose){ const os=O.run(cases); let nbad=0; cases.forEach((c,i)=>{ const r=mine(c), ex=diff(c,os[i],r); if(ex.length){ nbad++; if(verbose!==false) console.log(`case ${i} ${c.name||""}: `+ex.slice(0,8).join("; ")); } }); return {os,nbad}; }
module.exports={mine,diff,cmp,PR};
