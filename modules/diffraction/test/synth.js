// 시험용 합성 회절도 (분석 코드와 공유하지 않는 독립 구현). 정답 위치 = 결정학 코어의 정대축 반사 × 축척 · 회전 · (선택) 아핀 찌그러짐 T · 방사 왜곡 k
const fs=require("fs"), path=require("path");
require("../src/diff_sg.js"); require("../src/diff_scat.js"); const DC=require("../src/diff_core.js");
function rng(seed){ let a=seed>>>0; return ()=>{ a=(a+0x6D2B79F5)>>>0; let t=a; t=Math.imul(t^(t>>>15),t|1); t^=t+Math.imul(t^(t>>>7),t|61); return ((t^(t>>>14))>>>0)/4294967296; }; }
function poisson(R,lam){ if(lam<30){ const L=Math.exp(-lam); let k=0,p=1; do{ k++; p*=R(); }while(p>L); return k-1; } const u=R()||1e-12, v=R(); return Math.max(0,Math.round(lam+Math.sqrt(lam)*Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v))); }
const PH={ Si:DC.parseCIF(fs.readFileSync(path.join(__dirname,"cif/Si_Fd-3m.cif"),"utf8"))[0], Au:DC.parseCIF(fs.readFileSync(path.join(__dirname,"cif/Au_Fm-3m.cif"),"utf8"))[0] };

/* 스폿 회절도 합성: 정대축, 축척 s(1/nm/px), 중심, 회전 θ, 스폿 σ, 세기 ∝ |F|²·exp(−g²/g0²), 배경, 빔스토퍼, 글자 */
function synthSpots(o){ const R=rng(o.seed); const w=o.w, h=o.h; const img=new Float32Array(w*h);
  const zp=DC.zonePattern(o.phase,o.zone,o.gmax,{F:true}); const truth=[]; const c=Math.cos(o.rot), s=Math.sin(o.rot); let Imax=0; for(const p of zp.spots) Imax=Math.max(Imax,p.I);
  for(const p of zp.spots){ if(p.I<1e-6*Imax) continue; let u=(p.x*c-p.y*s)/o.scale, v=(p.x*s+p.y*c)/o.scale; if(o.T){ const uu=o.T[0][0]*u+o.T[0][1]*v, vv=o.T[1][0]*u+o.T[1][1]*v; u=uu; v=vv; } if(o.k){ const f=1+o.k*(u*u+v*v); u*=f; v*=f; } const x=o.cx+u, y=o.cy+v; if(x<5||y<5||x>w-5||y>h-5) continue; truth.push({x,y,h:p.h,g:p.g}); const A=o.peak*Math.sqrt(p.I/Imax)*Math.exp(-(p.g*p.g)/(o.g0*o.g0));
    const r=Math.ceil(4*o.sig); for(let yy=Math.max(0,Math.floor(y)-r);yy<=Math.min(h-1,Math.ceil(y)+r);yy++) for(let xx=Math.max(0,Math.floor(x)-r);xx<=Math.min(w-1,Math.ceil(x)+r);xx++) img[yy*w+xx]+=A*Math.exp(-((xx-x)**2+(yy-y)**2)/(2*o.sig*o.sig)); }
  // 직접빔 + 확산 배경
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const r2=(x-o.cx)**2+(y-o.cy)**2; img[y*w+x]+=o.bg*(1+8*Math.exp(-r2/(2*(w/6)**2)))+(o.direct?o.peak*20*Math.exp(-r2/(2*(o.sig*1.5)**2)):0); }
  for(let i=0;i<w*h;i++) img[i]=poisson(R,img[i]);
  if(o.stop){ // 빔스토퍼: 중심을 덮는 원판 + 가장자리까지 팔
    for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const dx=x-o.cx, dy=y-o.cy; const along=dx*Math.cos(o.stopAng)+dy*Math.sin(o.stopAng), across=-dx*Math.sin(o.stopAng)+dy*Math.cos(o.stopAng);
      if(dx*dx+dy*dy<(o.stopR)**2||(along>0&&Math.abs(across)<o.stopR*0.4)) img[y*w+x]=poisson(R,o.bg*0.05); } }
  if(o.text){ let mx=0; for(const v of img) mx=Math.max(mx,v); const sat=mx*1.2; for(let k=0;k<4;k++){ const x0=20+k*12; for(let y=20;y<34;y++) img[y*w+x0]=sat; for(let x=x0;x<x0+8;x++) img[20*w+x]=sat; } // 글자 비슷한 가는 획
    const by=h-30; for(let y=by;y<by+5;y++) for(let x=30;x<30+o.barLen;x++) img[y*w+x]=sat; }
  return {img,truth}; }


module.exports={rng,poisson,PH,synthSpots,DC};
