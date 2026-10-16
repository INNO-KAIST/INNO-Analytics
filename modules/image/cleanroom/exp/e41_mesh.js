// large meshes (one big biconnected block, n > FW_MAX): bounding-eccentricity path vs the oracle; also blobs with BLOCK_BUDGET=0
const C=require("./cmp.js"), CR=require("../askel_cr.js"); const IC=require("../../test/load_core.js");
let a=+(process.env.SEED||9); const rnd=()=>{ a=(a+0x6D2B79F5)>>>0; let t=a; t=Math.imul(t^(t>>>15),t|1); t^=t+Math.imul(t^(t>>>7),t|61); return ((t^(t>>>14))>>>0)/4294967296; };
const cases=[];
for(let k=0;k<+(process.env.NM||2);k++){ // jittered grid: polyline rows/cols with random offsets, then skeletonized
  const w=330,h=330, p={w,h,bit:8,data:new Float32Array(w*h)}; const step=11; const ox=[], oy=[];
  for(let i=0;i*step<w;i++) ox.push(i*step+4+Math.floor(rnd()*5)); for(let j=0;j*step<h;j++) oy.push(j*step+4+Math.floor(rnd()*5));
  const dot=(x,y)=>{ for(let dy=-1;dy<=1;dy++) for(let dx=-1;dx<=1;dx++){ const X=x+dx,Y=y+dy; if(X>=0&&Y>=0&&X<w&&Y<h) p.data[X+w*Y]=255; } };
  for(const y0 of oy){ let y=y0; for(let x=0;x<w;x++){ dot(x,y); if(rnd()<0.15) y+=rnd()<0.5?-1:1; } }
  for(const x0 of ox){ let x=x0; for(let y=0;y<h;y++){ dot(x,y); if(rnd()<0.15) x+=rnd()<0.5?-1:1; } }
  for(let i=0;i<w*h;i++) if(rnd()<0.002) p.data[i]=0;
  const sk=IC.skeletonize(p); cases.push({w,h,d:1,data:Uint8Array.from(sk.data,x=>x?255:0),lsp:true,cal:[0.731+0.1*rnd(),1.137,1],name:"jittered mesh "+k}); }
const mode=process.env.MODE||"default"; if(mode==="bounding"){ CR._t.setBlockBudget(0); CR._t.setFwMax(1); }
for(const c of cases){ const r=C.mine(c); console.log(c.name,"trees",r.trees.length,"edges",r.edges.length,"warn",r.warn.join(" | ")); }
const t0=Date.now(); const {nbad}=C.cmp(cases); console.log("mode",mode,"oracle s",(Date.now()-t0)/1000,"mismatching",nbad,"of",cases.length);
