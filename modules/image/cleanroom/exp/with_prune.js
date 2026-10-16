// wrapper exporting PRUNE=true so the harness adds its pruning cases (used before the flag is flipped in askel_cr.js)
const CR=require("../askel_cr.js"); module.exports=Object.assign({},CR,{PRUNE:true});
