const CR=require("../askel_cr.js");
function T(name,f){ const t=process.hrtime.bigint(); const r=f(); console.log(name.padEnd(48),(Number(process.hrtime.bigint()-t)/1e9).toFixed(3)+" s",JSON.stringify(r)); }
const mk=(w,h,d)=>({w,h,d,data:new Uint8Array(w*h*d)});
const orig=(v)=>{ const o=new Uint8Array(v.w*v.h*v.d); for(let i=0;i<o.length;i++) o[i]=(i*2654435761>>>24)&255; return {data:o}; };
const meshes={};
{ const v=mk(1300,1300,1); for(let y=0;y<v.h;y+=9) for(let x=0;x<v.w;x++) v.data[x+v.w*y]=255; for(let x=0;x<v.w;x+=9) for(let y=0;y<v.h;y++) v.data[x+v.w*y]=255; meshes["grid mesh 1300² (~21k vertices)"]=v; }
{ const n=256; let a=7; const rnd=()=>{ a=(a+0x6D2B79F5)>>>0; let t=a; t=Math.imul(t^(t>>>15),t|1); t^=t+Math.imul(t^(t>>>7),t|61); return ((t^(t>>>14))>>>0)/4294967296; }; const v=mk(n,n,n);
  for(let s=0;s<400;s++){ let x=Math.floor(rnd()*n),y=Math.floor(rnd()*n),z=Math.floor(rnd()*n); for(let k=0;k<600;k++){ v.data[x+n*(y+n*z)]=255; x=Math.max(0,Math.min(n-1,x+Math.floor(rnd()*3)-1)); y=Math.max(0,Math.min(n-1,y+Math.floor(rnd()*3)-1)); z=Math.max(0,Math.min(n-1,z+Math.floor(rnd()*3)-1)); } }
  meshes["256³ random walks"]=v; }
{ const v=mk(1500,1500,1); v.data.fill(255); meshes["1500² full"]=v; }
for(const [name,v] of Object.entries(meshes)){ const o=orig(v);
  for(const [P,pe] of [["shortest branch",false],["lowest intensity voxel",true],["lowest intensity branch",false]]) T(`${name} prune=${P}${pe?" +ends":""} +lsp`,()=>{ const r=CR.analyzeSkeleton(v,{pw:1,ph:1,pd:1},{lsp:true,prune:P,pruneEnds:pe,orig:o}); return [r.trees.length,r.edges.length,r.warn.length]; }); }
