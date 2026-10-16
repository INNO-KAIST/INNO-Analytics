// wrapper: force the non-matrix LSP path (FW_MAX=1; BB env → block budget) and log warnings, for harness runs
const CR=require("../askel_cr.js"); CR._t.setFwMax(1); if(process.env.BB) CR._t.setBlockBudget(+process.env.BB);
let k=0; const fs=require("fs");
module.exports={PRUNE:CR.PRUNE,analyzeSkeleton:(v,c,o)=>{ const r=CR.analyzeSkeleton(v,c,o); if(r.warn.length) fs.appendFileSync(process.env.WARNLOG||"/dev/null",`${k} ${r.warn.join(" | ")}\n`); k++; return r; }};
