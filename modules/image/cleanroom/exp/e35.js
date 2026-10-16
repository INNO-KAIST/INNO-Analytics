const H=require("./hyp.js");
const root=(T,E)=>{ const i=T.vCl.findIndex(c=>c<0); return i<0?0:i; };
const b={dfs:"stackE",order:"fwd",le:false,allRoots:true,root,method:+(process.env.P||2),vox:(e)=>e.slabs,vox2:(e)=>e.slabs,nb:"cal",cycleStartAll:true};
const V={zero:Object.assign({},b,{oob:"zero"}), skip:Object.assign({},b,{oob:"skip"}), zeroR2d:Object.assign({},b,{oob:"zero",R2d:true}), zeroZin2d:Object.assign({},b,{oob:"zero",zIn2d:true})};
console.log(H.evaluate(+(process.env.P||2),V,+(process.env.N||300),+(process.env.SEED||21),0).res);
