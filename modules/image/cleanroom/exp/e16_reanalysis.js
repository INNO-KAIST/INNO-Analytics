// Is the cycle-pruned output = plain analysis of (input minus removed voxels)?
const O=require("./oracle.js"), C=require("./cmp.js"), gen=require("./gen.js");
const P=+(process.env.P||1), N=+(process.env.N||300);
const cases=gen(+(process.env.SEED||1),N); const base=cases.map(c=>Object.assign({},c,{prune:0,pruneEnds:false})); const pr=cases.map(c=>Object.assign({},c,{prune:P,pruneEnds:false}));
const os=O.run(base.concat(pr)); let bad=0, cutN=0;
pr.forEach((c,i)=>{ const o0=os[i], o1=os[N+i]; const data=c.data.slice(); let nrem=0; for(const [p] of o0.tag) if(!o1.tag.has(p)){ data[p]=0; nrem++; } if(nrem) cutN++;
  // also voxels whose tag vanished but which were foreground: removed. Compare my plain analysis on the cut image with the oracle's pruned output
  const cc=Object.assign({},c,{data,prune:0}); const r=C.mine(cc); const ex=C.diff(cc,o1,r); if(ex.length){ bad++; if(bad<6) console.log(c.name,"removed",nrem,ex.slice(0,4).join("; ")); } });
console.log("prune",P,"cases",N,"with cuts",cutN,"mismatch",bad);
