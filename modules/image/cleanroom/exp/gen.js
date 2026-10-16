// random test volumes (same flavour as the harness generators), deterministic by seed
const IC=require("../../test/load_core.js");
module.exports=function gen(seed0,n){ let a=seed0>>>0; const rnd=()=>{ a=(a+0x6D2B79F5)>>>0; let t=a; t=Math.imul(t^(t>>>15),t|1); t^=t+Math.imul(t^(t>>>7),t|61); return ((t^(t>>>14))>>>0)/4294967296; };
  const blank=(w,h,d)=>({w,h,d,data:new Uint8Array(w*h*d)}); const put=(v,x,y,z)=>{ if(x>=0&&x<v.w&&y>=0&&y<v.h&&z>=0&&z<v.d) v.data[x+v.w*(y+v.h*z)]=255; };
  const out=[];
  for(let k=0;k<n;k++){ const kind=k%5; let c;
    if(kind<=1){ const w=20+Math.floor(rnd()*50), h=15+Math.floor(rnd()*40); const p={w,h,bit:8,data:new Float32Array(w*h)}; const nb=1+Math.floor(rnd()*5);
      for(let b=0;b<nb;b++){ const cx=rnd()*w, cy=rnd()*h, rx=2+rnd()*9, ry=2+rnd()*9; for(let y=0;y<h;y++) for(let x=0;x<w;x++) if(((x-cx)/rx)**2+((y-cy)/ry)**2<=1) p.data[y*w+x]=255; }
      if(kind===1) for(let i=0;i<w*h;i++) if(rnd()<0.03) p.data[i]=p.data[i]?0:255;
      const sk=IC.skeletonize(p); c={w,h,d:1,data:Uint8Array.from(sk.data,x=>x?255:0),cal:[0.3+rnd(),0.3+rnd(),1]}; }
    else if(kind===2){ const w=8+Math.floor(rnd()*30), h=6+Math.floor(rnd()*25); c=blank(w,h,1); const dens=0.05+rnd()*0.4; for(let i=0;i<w*h;i++) if(rnd()<dens) c.data[i]=255; c.cal=[1,1+rnd(),1]; }
    else if(kind===3){ const w=6+Math.floor(rnd()*14),h=6+Math.floor(rnd()*12),d=2+Math.floor(rnd()*6); c=blank(w,h,d); const dens=0.03+rnd()*0.15; for(let i=0;i<w*h*d;i++) if(rnd()<dens) c.data[i]=255; c.cal=[0.5+rnd(),0.5+rnd(),0.5+2*rnd()]; }
    else { const w=16+Math.floor(rnd()*14),h=14+Math.floor(rnd()*12),d=3+Math.floor(rnd()*6); c=blank(w,h,d); for(let s2=0;s2<2+Math.floor(rnd()*4);s2++){ let x=Math.floor(rnd()*w),y=Math.floor(rnd()*h),z=Math.floor(rnd()*d); for(let m=0;m<15+rnd()*40;m++){ put(c,x,y,z); x=Math.max(0,Math.min(w-1,x+Math.floor(rnd()*3)-1)); y=Math.max(0,Math.min(h-1,y+Math.floor(rnd()*3)-1)); z=Math.max(0,Math.min(d-1,z+Math.floor(rnd()*3)-1)); } } c.cal=[1,1,1+rnd()]; }
    const orig=new Uint8Array(c.w*c.h*c.d); const fx=rnd()*0.7, fy=rnd()*0.7, ph=rnd()*6; for(let i=0;i<orig.length;i++){ const x=i%c.w, y=((i/c.w)|0)%c.h, z=(i/(c.w*c.h))|0; orig[i]=Math.max(0,Math.min(255,Math.round(128+60*Math.sin(fx*x+ph)+50*Math.cos(fy*y+z)+20*(rnd()-0.5)))); }
    c.orig=orig; c.lsp=rnd()<0.7; c.name="g"+k+" kind"+kind; out.push(c); }
  return out; };
