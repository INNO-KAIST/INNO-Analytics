// Diffraction 패턴 읽기 대조 — 정답을 아는 합성 회절도 (스폿·링). 합성기는 시험 안에 따로 있다 (분석 코드와 공유 안 함).
const fs=require("fs"), path=require("path");
require("../src/diff_sg.js"); require("../src/diff_scat.js"); const DC=require("../src/diff_core.js"); const DP=require("../src/diff_pattern.js");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } if(!c||process.env.V) console.log((c?"  ✓ ":"  ✗ ")+m+(d?" — "+d:"")); };
const SEED=+(process.env.SEED||20261017);
function rng(seed){ let a=seed>>>0; return ()=>{ a=(a+0x6D2B79F5)>>>0; let t=a; t=Math.imul(t^(t>>>15),t|1); t^=t+Math.imul(t^(t>>>7),t|61); return ((t^(t>>>14))>>>0)/4294967296; }; }
function poisson(R,lam){ if(lam<30){ const L=Math.exp(-lam); let k=0,p=1; do{ k++; p*=R(); }while(p>L); return k-1; } const u=R()||1e-12, v=R(); return Math.max(0,Math.round(lam+Math.sqrt(lam)*Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v))); }
const PH={ Si:DC.parseCIF(fs.readFileSync(path.join(__dirname,"cif/Si_Fd-3m.cif"),"utf8"))[0], Au:DC.parseCIF(fs.readFileSync(path.join(__dirname,"cif/Au_Fm-3m.cif"),"utf8"))[0] };

/* 스폿 회절도 합성: 정대축, 축척 s(1/nm/px), 중심, 회전 θ, 스폿 σ, 세기 ∝ |F|²·exp(−g²/g0²), 배경, 빔스토퍼, 글자 */
function synthSpots(o){ const R=rng(o.seed); const w=o.w, h=o.h; const img=new Float32Array(w*h);
  const zp=DC.zonePattern(o.phase,o.zone,o.gmax,{F:true}); const truth=[]; const c=Math.cos(o.rot), s=Math.sin(o.rot); let Imax=0; for(const p of zp.spots) Imax=Math.max(Imax,p.I);
  for(const p of zp.spots){ if(p.I<1e-6*Imax) continue; const x=o.cx+(p.x*c-p.y*s)/o.scale, y=o.cy+(p.x*s+p.y*c)/o.scale; if(x<5||y<5||x>w-5||y>h-5) continue; truth.push({x,y,h:p.h,g:p.g}); const A=o.peak*Math.sqrt(p.I/Imax)*Math.exp(-(p.g*p.g)/(o.g0*o.g0));
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

console.log("== 스폿: Si [011] (빔스토퍼·글자·잡음, 회전)");
for(const [k,cfg] of [{zone:[0,1,1],rot:0.31,cx:251.37,cy:262.81},{zone:[0,0,1],rot:-0.77,cx:270.2,cy:240.6},{zone:[1,1,1],rot:1.13,cx:255.5,cy:255.5}].entries()){
  const o=Object.assign({seed:SEED+k,w:512,h:512,phase:PH.Si,gmax:14,scale:0.045,sig:1.6,peak:900,g0:12,bg:6,direct:true,stop:true,stopR:14,stopAng:2.2,text:true,barLen:120},cfg);
  const {img,truth}=synthSpots(o); const bs=DP.findBeamStop(img,o.w,o.h); const ov=DP.findOverlay(img,o.w,o.h); const mask=new Uint8Array(o.w*o.h); for(let i=0;i<mask.length;i++) mask[i]=bs.mask[i]|ov.mask[i];
  const sb=DP.scaleBar(ov,o.w,o.h); ok(sb&&Math.abs(sb.lenPx-o.barLen)<=1,`[${o.zone}] 스케일바 길이 ${sb&&sb.lenPx} px (정답 ${o.barLen})`);
  const {spots}=DP.findSpots(img,o.w,o.h,{mask}); const vis=truth.filter(t=>!mask[Math.round(t.y)*o.w+Math.round(t.x)]);
  let hit=0, err=0; for(const t of vis){ let b=1e9; for(const p of spots) b=Math.min(b,Math.hypot(p.x-t.x,p.y-t.y)); if(b<2){ hit++; err=Math.max(err,b); } }
  ok(hit>=0.9*vis.length,`[${o.zone}] 스폿 ${spots.length} 찾음 · 보이는 정답 ${vis.length} 중 ${hit} (최악 ${err.toFixed(2)} px)`);
  const cen=DP.findCenter(spots,o.w,o.h,{beamStop:bs.mask}); ok(cen&&Math.hypot(cen.x-o.cx,cen.y-o.cy)<1.0,`[${o.zone}] 중심 (Friedel) 오차 ${cen&&Math.hypot(cen.x-o.cx,cen.y-o.cy).toFixed(3)} px`);
  const L=DP.fitLattice(spots,cen,{size:o.w}); // 정답 기준: 가장 짧은 두 허용 반사 (Gauss 축약)
  const zp=DC.zonePattern(o.phase,o.zone,o.gmax,{F:true}); let Imax=0; for(const p of zp.spots) Imax=Math.max(Imax,p.I); const strong=zp.spots.filter(p=>p.I>1e-6*Imax); const v=strong.map(p=>[p.x/o.scale,p.y/o.scale]);
  let tb=null; for(let i=0;i<v.length&&!tb;i++) for(let j=i+1;j<v.length;j++){ const d=Math.abs(v[i][0]*v[j][1]-v[i][1]*v[j][0]); if(d>1e-6){ // 기본 격자 = 모든 강한 점을 정수로
        const D=v[i][0]*v[j][1]-v[i][1]*v[j][0]; const all=v.every(q=>{ const n1=(q[0]*v[j][1]-q[1]*v[j][0])/D, n2=(v[i][0]*q[1]-v[i][1]*q[0])/D; return Math.abs(n1-Math.round(n1))<1e-6&&Math.abs(n2-Math.round(n2))<1e-6; }); if(all){ tb=DP.reduce2(v[i],v[j]); break; } } }
  const tl=[Math.hypot(...tb[0]),Math.hypot(...tb[1]),Math.hypot(tb[0][0]+tb[1][0],tb[0][1]+tb[1][1])];
  ok(L&&Math.abs(L.len[0]-tl[0])<0.08&&Math.abs(L.len[1]-tl[1])<0.08&&Math.abs(L.len[2]-tl[2])<0.1,`[${o.zone}] R1·R2·R3 ${L&&L.len.map(x=>x.toFixed(3)).join(" / ")} px (정답 ${tl.map(x=>x.toFixed(3)).join(" / ")})`);
  ok(L&&L.ang12>=90-1e-9&&L.len[2]<=Math.hypot(L.R1[0]-L.R2[0],L.R1[1]-L.R2[1]),`[${o.zone}] REDI 규약: ∠R1·R2 = ${L&&L.ang12.toFixed(2)}° ≥ 90, R3 = 짧은 대각선`);
  ok(L&&Math.hypot(L.center.x-o.cx,L.center.y-o.cy)<0.3,`[${o.zone}] 격자 맞춤 중심 오차 ${L&&Math.hypot(L.center.x-o.cx,L.center.y-o.cy).toFixed(3)} px`);
  // 거짓 스폿 섞기 (문턱 근처 SNR 6–8, 15 개) → 같은 격자
  { const R2=rng(o.seed+100); const fake=spots.concat(Array.from({length:15},()=>({x:40+R2()*(o.w-80),y:40+R2()*(o.h-80),I:50+R2()*60,snr:6+2*R2()}))); const L2=DP.fitLattice(fake,cen,{size:o.w});
    ok(L2&&Math.abs(L2.len[0]-tl[0])<0.1&&Math.abs(L2.len[1]-tl[1])<0.1,`[${o.zone}] 거짓 스폿 15 개 섞어도 같은 격자 (${L2&&L2.len.map(x=>x.toFixed(2)).join(" / ")})`); }
  // 불확도가 정직한지: 실제 오차 ≤ 3σ
  ok(L&&Math.abs(L.len[0]-tl[0])<=3*L.sig[0]+0.01&&Math.abs(L.len[1]-tl[1])<=3*L.sig[1]+0.01,`[${o.zone}] σ(R1,R2) ${L&&L.sig.map(x=>x.toFixed(4)).join(", ")} px 가 실제 오차를 덮음`);
}

console.log("== 링: Au (타원 3 %, 빔스토퍼)");
{ const R=rng(SEED+9); const w=600, h=600, cx=311.4, cy=288.9, scale=0.02, ell=1.03, phi=0.6; const rings=DC.ringList(PH.Au,10); const img=new Float32Array(w*h);
  const cs=Math.cos(phi), sn=Math.sin(phi);
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const dx=x-cx, dy=y-cy; const u=dx*cs+dy*sn, v=-dx*sn+dy*cs; const r=Math.hypot(u/ell,v); // 장축(u) 방향으로 ell 배 늘어난 상
    let I=8+400*Math.exp(-r/60); for(const g of rings){ const rr=g.g/scale; if(rr>420) continue; I+=Math.min(3000,g.m*g.I/200)*Math.exp(-((r-rr)**2)/(2*1.8*1.8))/(rr/50); } img[y*w+x]=poisson(R,I); }
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const dx=x-cx, dy=y-cy; if(dx*dx+dy*dy<22*22||(dx<0&&Math.abs(dy)<8)) img[y*w+x]=poisson(R,0.3); }
  const bs=DP.findBeamStop(img,w,h); const r1=rings[0].g/scale; // 첫 링 위치를 대충 안다고 두고(사용자가 찍음) 광선 탐색
  const pts=DP.ringPoints(img,w,h,{x:cx+4,y:cy-3},r1*0.85,r1*1.15,{mask:bs.mask}); const el=DP.fitEllipse(pts);
  ok(el&&Math.hypot(el.cx-cx,el.cy-cy)<0.3,`타원 중심 오차 ${el&&Math.hypot(el.cx-cx,el.cy-cy).toFixed(3)} px (점 ${pts.length})`);
  ok(el&&Math.abs(el.a/el.b-ell)<0.003,`장단축비 ${el&&(el.a/el.b).toFixed(4)} (정답 ${ell})`);
  let dphi=el?Math.abs(el.phi-phi)%Math.PI:9; dphi=Math.min(dphi,Math.PI-dphi); ok(dphi<0.03,`장축 방향 오차 ${(dphi*180/Math.PI).toFixed(2)}°`);
  ok(el&&Math.abs(el.b-r1)<0.2,`단축 반지름 ${el&&el.b.toFixed(2)} (정답 ${r1.toFixed(2)})`);
  const az=DP.azimuthal(img,w,h,el,{mask:bs.mask,rmax:420}); const bg=DP.background1d(az.I,24); const pk=DP.peaks1d(az.r,az.I,bg,{k:6}).peaks;
  let bad=0, nm=0, worst=0; for(const g of rings){ const rr=g.g/scale; if(rr>400||rr<30) continue; nm++; let b=1e9; for(const p of pk) b=Math.min(b,Math.abs(p.r-rr)); if(b>0.4) bad++; worst=Math.max(worst,b); }
  ok(!bad,`방위 적분 봉우리 ${pk.length} · Au 링 ${nm} 개 모두 ±0.4 px (최악 ${worst.toFixed(3)} px)`);
}
console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
