// LSP on large graphs (n > FW_MAX) vs the oracle (whose O(n^3) matrix takes a while)
const C=require("./cmp.js"), CR=require("../askel_cr.js"); const IC=require("../../test/load_core.js");
let seed=+(process.env.SEED||4242); const rnd=()=>((seed=(seed*1103515245+12345)&0x7fffffff)/0x7fffffff);
const cases=[];
function skelBlobs(w,h,nb,holes){ const p={w,h,bit:8,data:new Float32Array(w*h)};
  for(let b=0;b<nb;b++){ const cx=rnd()*w, cy=rnd()*h, rx=4+rnd()*30, ry=4+rnd()*30; for(let y=Math.max(0,Math.floor(cy-ry));y<Math.min(h,cy+ry+1);y++) for(let x=Math.max(0,Math.floor(cx-rx));x<Math.min(w,cx+rx+1);x++) if(((x-cx)/rx)**2+((y-cy)/ry)**2<=1) p.data[y*w+x]=255; }
  for(let i=0;i<w*h;i++) if(rnd()<holes) p.data[i]=0;
  const sk=IC.skeletonize(p); return {w,h,d:1,data:Uint8Array.from(sk.data,x=>x?255:0)}; }
const N=+(process.env.NC||4);
for(let k=0;k<N;k++){ const v=skelBlobs(+(process.env.BW||300),+(process.env.BH||260),+(process.env.NB||120),+(process.env.HOLES||0.01)); cases.push(Object.assign(v,{lsp:true,cal:[0.5+rnd(),0.5+rnd(),1],name:"blobs+holes "+k})); }
// a pure tree: random branching walk (no loops) on a big canvas, with irrational calibration
for(let k=0;k<(+process.env.NT||1);k++){ const w=+(process.env.TW||900),h=w, data=new Uint8Array(w*h); const c0=w>>1; const occ=(x,y)=>x<1||y<1||x>=w-1||y>=h-1||data[x+w*y];
  const tips=[[c0,c0,0]]; data[c0+w*c0]=255; const dirs=[[1,0],[0,1],[-1,0],[0,-1]];
  for(let it=0;it<+(process.env.TI||200000)&&tips.length;it++){ const ti=Math.floor(rnd()*tips.length); const [x,y,dd]=tips[ti]; const d=rnd()<0.85?dd:Math.floor(rnd()*4); const [dx,dy]=dirs[d]; const nx=x+2*dx, ny=y+2*dy;
    // free if the 2-step corridor and its sides are empty
    let free=!occ(nx,ny)&&!occ(x+dx,y+dy); for(const [ex,ey] of [[nx+dy,ny+dx],[nx-dy,ny-dx],[nx+dx,ny+dy],[x+dx+dy,y+dy+dx],[x+dx-dy,y+dy-dx]]) if(occ(ex,ey)) free=false;
    if(!free){ if(rnd()<0.3) tips.splice(ti,1); continue; }
    data[x+dx+w*(y+dy)]=255; data[nx+w*ny]=255; tips[ti]=[nx,ny,d]; if(rnd()<0.08) tips.push([nx,ny,(d+1+2*Math.floor(rnd()*2))%4]); }
  cases.push({w,h,d:1,data,lsp:true,cal:[0.7310,1.1370,1],name:"random tree "+k}); }
for(const c of cases){ const r=C.mine(c); c._nv=r.trees.map((t,i)=>new Set(r.edges.filter(e=>e.tree===i).flatMap(e=>[e.v1.join(),e.v2.join()])).size); c._warn=r.warn; }
cases.forEach((c,i)=>console.log(i,c.name,c.w+"x"+c.h,"trees",c._nv.length,"max vertices/tree (distinct first voxels)",Math.max(...c._nv),"warn",c._warn.length));
if(process.env.DRY) process.exit(0);
const t0=Date.now(); const {os,nbad}=C.cmp(cases); console.log("oracle time",(Date.now()-t0)/1000,"s; mismatching cases",nbad,"of",cases.length);
