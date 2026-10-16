const CR=require("../askel_cr.js");
function T(name,f){ const t=process.hrtime.bigint(); const r=f(); console.log(name, (Number(process.hrtime.bigint()-t)/1e9).toFixed(3)+" s", r); }
const run=(v,o)=>CR.analyzeSkeleton(v,{pw:1,ph:1,pd:1},o);
{ const w=1500,h=1500; const v={w,h,d:1,data:new Uint8Array(w*h).fill(255)}; T("1500x1500 full",()=>{ const r=run(v,{lsp:false}); return [r.trees.length,r.trees[0].junctions,r.edges.length]; });
  T("1500x1500 full +lsp",()=>{ const r=run(v,{lsp:true}); return [r.trees[0].lsp]; }); }
{ const w=2000,h=2000; const data=new Uint8Array(w*h); for(let y=0;y+4<h;y+=4){ for(let x=1;x<w-1;x++) data[x+w*y]=255; const cx=((y>>2)&1)?0:w-1; for(let k=1;k<4;k++) data[cx+w*(y+k)]=255; }
  const v={w,h,d:1,data}; T("2000x2000 serpentine +lsp",()=>{ const r=run(v,{lsp:true}); return [r.trees.length,r.edges.length,r.trees[0].lsp]; }); }
{ const n=256; let seed=7; const rnd=()=>((seed=(seed*1103515245+12345)&0x7fffffff)/0x7fffffff); const data=new Uint8Array(n*n*n);
  for(let s=0;s<400;s++){ let x=Math.floor(rnd()*n),y=Math.floor(rnd()*n),z=Math.floor(rnd()*n); for(let k=0;k<600;k++){ data[x+n*(y+n*z)]=255; x=Math.max(0,Math.min(n-1,x+Math.floor(rnd()*3)-1)); y=Math.max(0,Math.min(n-1,y+Math.floor(rnd()*3)-1)); z=Math.max(0,Math.min(n-1,z+Math.floor(rnd()*3)-1)); } }
  const v={w:n,h:n,d:n,data}; let c=0; for(const q of data) if(q) c++;
  T(`256^3 sparse random walks (${c} fg) +lsp`,()=>{ const r=run(v,{lsp:true}); return [r.trees.length,r.edges.length,Math.max(...r.trees.map(t=>t.lsp)).toFixed(3),r.warn.length]; }); }
{ // binary-tree-like skeleton with ~20000 vertices in one tree: comb of combs
  const w=4100,h=420; const data=new Uint8Array(w*h); for(let x=0;x<w;x++) data[x+w*10]=255; // spine
  for(let x=5;x<w-5;x+=4){ for(let y=11;y<410;y++) data[x+w*y]=255; for(let y=20;y<405;y+=8) { data[x+1+w*y]=255; } }
  const v={w,h,d:1,data}; T("comb (~20k vertices) +lsp",()=>{ const r=run(v,{lsp:true}); return [r.trees.length,r.edges.length,r.trees[0].lsp,r.warn]; }); }
{ // grid-of-loops skeleton: big biconnected block
  const w=1300,h=1300; const data=new Uint8Array(w*h); for(let y=0;y<h;y+=9) for(let x=0;x<w;x++) data[x+w*y]=255; for(let x=0;x<w;x+=9) for(let y=0;y<h;y++) data[x+w*y]=255;
  const v={w,h,d:1,data}; T("grid mesh +lsp",()=>{ const r=run(v,{lsp:true}); return [r.trees.length,r.edges.length,r.trees[0].lsp,r.warn]; }); }
