// sensitivity: crop the region around a mismatch, perturb each orig pixel, see which removal results
const gen=require("./gen.js"), O=require("./oracle.js");
const [seed,N,idx,x0,y0,x1,y1]=process.argv.slice(2).map(Number); const c=gen(seed,N)[idx];
const w=x1-x0+1, h=y1-y0+1; const data=new Uint8Array(w*h), orig=new Uint8Array(w*h);
for(let y=0;y<h;y++) for(let x=0;x<w;x++){ data[x+w*y]=c.data[(x+x0)+c.w*(y+y0)]; orig[x+w*y]=c.orig[(x+x0)+c.w*(y+y0)]; }
const cases=[{w,h,d:1,data,orig,prune:0,lsp:false},{w,h,d:1,data,orig,prune:2,lsp:false}];
const val=+(process.env.VAL||0);
for(let i=0;i<w*h;i++){ const o=orig.slice(); o[i]=val; cases.push({w,h,d:1,data,orig:o,prune:2,lsp:false,i}); }
const os=O.run(cases); const rem=(o)=>[...os[0].tag.keys()].filter(p=>!o.tag.has(p)).map(p=>(p%w+x0)+","+(((p/w)|0)+y0)).sort().join(" ");
const r0=rem(os[1]); console.log("crop removed:",r0);
for(let k=2;k<cases.length;k++){ const r=rem(os[k]); if(r!==r0){ const i=cases[k].i; console.log("pixel",(i%w+x0)+","+(((i/w)|0)+y0),"=",val,"→ removed",r); } }
