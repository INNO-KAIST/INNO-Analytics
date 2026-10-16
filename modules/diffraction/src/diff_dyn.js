/* ============================================================
   INNO Analytics — Diffraction: 동역학 회절 스폿 세기 (Bloch 파, 정대축 0 차 Laue 영역). DOM 없음, Worker 가능.
   · 관례 (De Graef 2003 · Spence & Zuo 1992): 파동 exp(+2πi k·r), 2π 없는 역공간, 길이 nm.
     Schrödinger [∇² + 4π²(K0² + U(r))]ψ = 0,  U(r) = Σ U_g exp(2πi g·r),  U_g = (2 m e / h²) V_g.
     V_g = (h²/(2π m0 e)) conj(F_g)/Ω  (F_g = Σ f_e e^{+2πi g·x}, f_e Å — IT C 4.3.2.2, 비상대론 값)  ⇒  U_g = γ conj(F_g) / (π Ω)  [Å⁻²].
   · 구조 행렬 (고에너지 근사, 표면 ⟂ 빔):  H_gh = s_g δ_gh + U_{g−h}/(2K),  γ^j C^j = H C^j,  K² = K0² + U_0.
     ψ_g(t) = Σ_j C_0^{j*} C_g^j exp(2πi γ^j t),  I_g = |ψ_g|².  s_g = −(2 k_t·g + g²)/(2K) (정대축 k_t = 0 이면 −g²/2K).
   · 빔: 강한 빔 (|g| 껍질 단위, 기본 ≤ 200 개·≤ 35 1/nm) 은 정확히, 그 밖 ≤ 50 1/nm 는 Bethe 퍼텐셜 H_gh += Σ_w H_gw H_wh/(−s_w)
     (Zuo & Weickenmeier 1995). 출력 범위 안인데 강한 빔이 아닌 것 (큰 세포) 은 Bloch 파에 대한 1 차 섭동.
   · 흡수 (선택): 비례 모형 U'_g = κ U_g. 1 차 섭동은 비례 모형에서 크게 틀려 (Si [011] κ 0.1 에서 48 %) 쓰지 않고,
     비 Hermite H 를 exp(2πiδH) (축척·제곱 + Taylor) 걸음으로 정확히 푼다.
   · 고유값 풀이: 실대칭 Householder 3 대각화 + QL (JAMA tred2/tql2). 복소 Hermite 는 2N 실대칭 [[A,−B],[B,A]] 로.
     반전 중심이 있으면 원점을 옮겨 U 를 실수로 (세기 불변) → N 실대칭.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const DC=G.DIFF_CORE||(typeof require!=="undefined"?require("./diff_core.js"):null);
const DD={};
const M0C2=510998.95;           // eV (CODATA 2018)
const HC=12398.419843320026;    // eV·Å
const TWO_PI=2*Math.PI;
const now=()=>(typeof performance!=="undefined"&&performance.now)?performance.now():Date.now();
const {dot,cross,norm}=DC.vec;

/* 상대론 파장 λ (Å), γ = 1 + E/m0c² */
DD.wavelength=kV=>{ const E=kV*1e3; return HC/Math.sqrt(E*(2*M0C2+E)); };
DD.gammaRel=kV=>1+kV*1e3/M0C2;

/* ---------------- 실대칭 고유값 (JAMA tred2 + tql2) ----------------
   A: Float64Array n×n (행 우선, 대칭). 반환 d (오름차순), C: 열 배열 C[j] = 고유벡터 j */
function symEig(A,n){
  const C=new Array(n); for(let j=0;j<n;j++){ const c=new Float64Array(n); for(let i=0;i<n;i++) c[i]=A[i*n+j]; C[j]=c; }
  const d=new Float64Array(n), e=new Float64Array(n);
  // V(r,c) = C[c][r]
  for(let j=0;j<n;j++) d[j]=C[j][n-1];
  for(let i=n-1;i>0;i--){ let scale=0, h=0; for(let k=0;k<i;k++) scale+=Math.abs(d[k]);
    if(scale===0){ e[i]=d[i-1]; for(let j=0;j<i;j++){ d[j]=C[j][i-1]; C[j][i]=0; C[i][j]=0; } }
    else { for(let k=0;k<i;k++){ d[k]/=scale; h+=d[k]*d[k]; }
      let f=d[i-1], g=Math.sqrt(h); if(f>0) g=-g; e[i]=scale*g; h-=f*g; d[i-1]=f-g; for(let j=0;j<i;j++) e[j]=0;
      for(let j=0;j<i;j++){ f=d[j]; C[i][j]=f; const cj=C[j]; g=e[j]+cj[j]*f; for(let k=j+1;k<=i-1;k++){ g+=cj[k]*d[k]; e[k]+=cj[k]*f; } e[j]=g; }
      f=0; for(let j=0;j<i;j++){ e[j]/=h; f+=e[j]*d[j]; } const hh=f/(h+h); for(let j=0;j<i;j++) e[j]-=hh*d[j];
      for(let j=0;j<i;j++){ f=d[j]; g=e[j]; const cj=C[j]; for(let k=j;k<=i-1;k++) cj[k]-=(f*e[k]+g*d[k]); d[j]=cj[i-1]; cj[i]=0; } }
    d[i]=h; }
  for(let i=0;i<n-1;i++){ C[i][n-1]=C[i][i]; C[i][i]=1; const h=d[i+1]; const ci1=C[i+1];
    if(h!==0){ for(let k=0;k<=i;k++) d[k]=ci1[k]/h; for(let j=0;j<=i;j++){ const cj=C[j]; let g=0; for(let k=0;k<=i;k++) g+=ci1[k]*cj[k]; for(let k=0;k<=i;k++) cj[k]-=g*d[k]; } }
    for(let k=0;k<=i;k++) ci1[k]=0; }
  for(let j=0;j<n;j++){ d[j]=C[j][n-1]; C[j][n-1]=0; } C[n-1][n-1]=1; e[0]=0;
  // tql2
  for(let i=1;i<n;i++) e[i-1]=e[i]; e[n-1]=0; let f=0, tst1=0; const eps=Math.pow(2,-52);
  for(let l=0;l<n;l++){ tst1=Math.max(tst1,Math.abs(d[l])+Math.abs(e[l])); let m=l; while(m<n){ if(Math.abs(e[m])<=eps*tst1) break; m++; }
    if(m>l){ let iter=0; do{ if(++iter>200) throw new Error("고유값 풀이가 수렴하지 않습니다");
        let g=d[l]; let p=(d[l+1]-g)/(2*e[l]); let r=Math.hypot(p,1); if(p<0) r=-r; d[l]=e[l]/(p+r); d[l+1]=e[l]*(p+r); const dl1=d[l+1]; let h=g-d[l];
        for(let i=l+2;i<n;i++) d[i]-=h; f+=h; p=d[m]; let c=1, c2=c, c3=c; const el1=e[l+1]; let s=0, s2=0;
        for(let i=m-1;i>=l;i--){ c3=c2; c2=c; s2=s; g=c*e[i]; h=c*p; r=Math.hypot(p,e[i]); e[i+1]=s*r; s=e[i]/r; c=p/r; p=c*d[i]-s*g; d[i+1]=h+s*(c*g+s*d[i]);
          const a=C[i], b=C[i+1]; for(let k=0;k<n;k++){ const hb=b[k]; b[k]=s*a[k]+c*hb; a[k]=c*a[k]-s*hb; } }
        p=-s*s2*c3*el1*e[l]/dl1; e[l]=s*p; d[l]=c*p; } while(Math.abs(e[l])>eps*tst1); }
    d[l]+=f; e[l]=0; }
  // 정렬
  const idx=Array.from({length:n},(_,i)=>i).sort((a,b)=>d[a]-d[b]); const d2=new Float64Array(n), C2=new Array(n); for(let i=0;i<n;i++){ d2[i]=d[idx[i]]; C2[i]=C[idx[i]]; }
  return {d:d2,C:C2}; }
DD._symEig=symEig;

/* ---------------- 구조인자 (복소) ----------------
   원점 이동 r0 (분율), DW: opt.dw === false 이면 끔, opt.Bdefault: B = 0 인 원자에 넣을 값 (Å², 기본 0.5), opt.Bset: 모든 원자 B 고정 (검증용) */
function sfMaker(phase,opt){ const uc=DC.expandAtoms(phase); const cell=phase.cell; const r0=opt.r0||[0,0,0]; const dw=opt.dw!==false; const Bd=opt.Bdefault==null?0.5:+opt.Bdefault||0; const Bs=opt.Bset;
  const els=[...new Set(uc.map(a=>a.el))];
  const at=uc.map(a=>({ei:els.indexOf(a.el),x:a.x-r0[0],y:a.y-r0[1],z:a.z-r0[2],w:a.occ,B:!dw?0:Bs!=null?+Bs:(a.B||Bd)}));
  return function(h){ const g=DC.gLen(cell,h), s=g/20, s2=s*s; const f=els.map(el=>DC.fElectron(el,s)); let re=0, im=0;
    for(const a of at){ const w=a.w*f[a.ei]*(a.B?Math.exp(-a.B*s2):1); const ph=TWO_PI*(h[0]*a.x+h[1]*a.y+h[2]*a.z); re+=w*Math.cos(ph); im+=w*Math.sin(ph); }
    return [re,im]; }; }
/* 반전 중심 (R = −I 인 연산 x' = −x + t → 중심 t/2) */
function inversionCenter(ops){ for(const op of ops||[]){ const R=op.R; if(R[0]===-1&&R[4]===-1&&R[8]===-1&&!R[1]&&!R[2]&&!R[3]&&!R[5]&&!R[6]&&!R[7]) return op.t.map(x=>x/2); } return null; }
/* 격자 중심화 평행이동 (R = I, t ≠ 0) — 이 소멸 자리는 동역학적으로도 결코 들뜨지 않으므로 빔에서 뺀다 */
function centerings(ops){ const T=[]; for(const op of ops||[]){ const R=op.R; if(R[0]===1&&R[4]===1&&R[8]===1&&!R[1]&&!R[2]&&!R[3]&&!R[5]&&!R[6]&&!R[7]&&(op.t[0]||op.t[1]||op.t[2])) T.push(op.t); } return T; }
const nonInt=x=>Math.abs(x-Math.round(x))>1e-6;

/* 복소 보조 */
function cexpm1OverX(xr,xi){ // (e^x − 1)/x
  const m=Math.hypot(xr,xi); if(m<1e-4){ // 1 + x/2 + x²/6
    return [1+xr/2+(xr*xr-xi*xi)/6, xi/2+(2*xr*xi)/6]; }
  const ex=Math.exp(xr); const er=ex*Math.cos(xi)-1, ei=ex*Math.sin(xi); const d=xr*xr+xi*xi; return [(er*xr+ei*xi)/d,(ei*xr-er*xi)/d]; }

/* ---------------- 빔 고르기 ---------------- */
function zoneSetup(phase,uvw,opt){ const cell=phase.cell; const z=DC.reduceInt(uvw); const pb=DC.planeBasis(z); if(!pb) throw new Error("정대축 ["+uvw+"] 의 평면 기저를 못 찾았습니다");
  const [p,q]=pb; const Gp=DC.gvec(cell,p), Gq=DC.gvec(cell,q); const area=norm(cross(Gp,Gq)); const gmax=opt.gmax;
  const Mm=Math.ceil(gmax*norm(Gq)/area)+1, Nm=Math.ceil(gmax*norm(Gp)/area)+1; const T=centerings(phase.ops);
  const pts=[];
  for(let m=-Mm;m<=Mm;m++) for(let n=-Nm;n<=Nm;n++){ const h=[m*p[0]+n*q[0],m*p[1]+n*q[1],m*p[2]+n*q[2]]; const g=(m||n)?DC.gLen(cell,h):0; if(g>gmax*(1+1e-12)) continue;
    if(T.some(t=>nonInt(h[0]*t[0]+h[1]*t[1]+h[2]*t[2]))) continue; pts.push({h,m,n,g}); }
  // 2D 좌표: DC.zonePattern 과 같은 e1 (첫 허용 반사, 같은 정렬) 또는 opt.ref
  const nz=DC.rvec(cell,z); const zz=nz.map(x=>x/norm(nz)); let e1=null;
  if(opt.ref){ const g=DC.gvec(cell,opt.ref); e1=g.map(x=>x/norm(g)); }
  else { const al=pts.filter(r=>r.g>0&&!DC.absent(phase.ops,r.h)).sort((a,b)=>a.g-b.g||(b.h[0]-a.h[0])||(b.h[1]-a.h[1])||(b.h[2]-a.h[2])); const r=al[0]||pts.find(r=>r.g>0); const g=DC.gvec(cell,r.h); e1=g.map(x=>x/norm(g)); }
  const e2=cross(zz,e1);
  for(const r of pts){ const g=DC.gvec(cell,r.h); r.x=dot(g,e1); r.y=dot(g,e2); }
  pts.sort((a,b)=>a.g-b.g||(b.h[0]-a.h[0])||(b.h[1]-a.h[1])||(b.h[2]-a.h[2]));
  return {z,p,q,pts,e1,e2,axis:zz}; }

/* ---------------- 본 계산 ----------------
   zoneIntensities(phase, uvw, opt)
   opt: kV (200) · thickness (nm, 수 또는 배열; 기본 10–100 nm 10 간격) · gmax (출력 범위 1/nm, 12)
        maxBeams (강한 빔 상한 200, 껍질 단위) · gStrong (강한 빔 반지름 상한, 기본 max(gmax, 35)) · gBethe (Bethe 섭동 빔 반지름 50, 0 = 끔) · maxBethe (3000)
        tilt [tx,ty] mrad (입사 빔 접선 성분, 2D 좌표 e1·e2 방향) · absorption κ (0 = 없음; > 0 이면 행렬 지수로 정확히)
        weak (true: g_cut < g ≤ gmax 출력 빔은 1 차 섭동) · beams ([[h,k,l]…] 강한 빔 직접 지정 — 이때 weak·Bethe 는 기본 끔)
        dw (true) · Bdefault (B = 0 원자에 넣을 B, Å², 기본 0.5) · Bset (모든 원자 B, 검증용) · origin ('auto' | 'cif') · gm (true) · gmProbe (nm) · absStep (0.125 nm) · ref (e1 기준 hkl) */
DD.zoneIntensities=function(phase,uvw,opt){ const T0=now(); opt=Object.assign({},opt||{});
  if(!phase||!phase.atoms||!phase.atoms.length) throw new Error("원자 좌표가 없는 상은 동역학 세기를 계산할 수 없습니다");
  const kV=+opt.kV||200; const lamA=DD.wavelength(kV), gam=DD.gammaRel(kV); const K0=10/lamA; // nm⁻¹
  let thick=opt.thickness==null?[10,20,30,40,50,60,70,80,90,100]:(Array.isArray(opt.thickness)?opt.thickness.slice():[+opt.thickness]); thick=thick.map(Number);
  const explicit=!!opt.beams; const gmax=opt.gmax||12, maxBeams=opt.maxBeams||200; const kap=+opt.absorption||0;
  const gB=opt.gBethe!=null?+opt.gBethe:(explicit?0:50); const maxBethe=opt.maxBethe||3000;
  const warn=[]; if(opt.includeHOLZ) warn.push("HOLZ 는 아직 계산하지 않습니다 (0 차 Laue 영역만)");
  // 원점: 반전 중심이 있으면 거기로 (U 실수, 세기 불변)
  const inv=opt.origin==="cif"?null:inversionCenter(phase.ops); const r0=inv||[0,0,0];
  const SF=sfMaker(phase,{r0,dw:opt.dw,Bdefault:opt.Bdefault,Bset:opt.Bset});
  const Ufac=gam/(Math.PI*phase.cell.V*1000)*100; // F(Å) → U (nm⁻²),  U_g = γ conj(F)/(πΩ)
  const gS=opt.gStrong||Math.max(gmax,35); const zs=zoneSetup(phase,uvw,{gmax:Math.max(gmax,gB,gS),ref:opt.ref}); const pts=zs.pts; const p=zs.p, q=zs.q;
  // U 저장: 평면 정수 좌표 차 (dm,dn) 의 표
  let Mx=0, Nx=0; for(const r of pts){ Mx=Math.max(Mx,Math.abs(r.m)); Nx=Math.max(Nx,Math.abs(r.n)); }
  if(explicit) for(const h of opt.beams) for(const v of h) Mx=Math.max(Mx,4*Math.abs(v)), Nx=Math.max(Nx,4*Math.abs(v));
  const OM=2*Mx+1, ON=2*Nx+1, WN=2*ON+1; const UT=new Float64Array((2*OM+1)*WN*2).fill(NaN);
  const Uk=(dm,dn)=>{ const k=((dm+OM)*WN+(dn+ON))*2; if(UT[k]!==UT[k]){ const F=SF([dm*p[0]+dn*q[0],dm*p[1]+dn*q[1],dm*p[2]+dn*q[2]]); UT[k]=F[0]*Ufac; UT[k+1]=-F[1]*Ufac; } return k; }; // 표 위치 (할당 없음)
  const Uof=(dm,dn)=>{ const k=Uk(dm,dn); return [UT[k],UT[k+1]]; };
  const U0=Uof(0,0)[0]; const K=Math.sqrt(K0*K0+U0);
  // 강한 빔 · 출력 바깥 빔 · Bethe 빔
  let strong;
  if(explicit){ const key=h=>h.join(","); const have=new Set(); strong=[];
    for(const h of [[0,0,0]].concat(opt.beams)){ const k=key(h); if(have.has(k)) continue; have.add(k);
      const pp=dot(p,p), qq=dot(q,q), pq=dot(p,q), det=pp*qq-pq*pq, hp=dot(h,p), hq=dot(h,q); const m=Math.round((hp*qq-hq*pq)/det), n=Math.round((hq*pp-hp*pq)/det);
      if(m*p[0]+n*q[0]!==h[0]||m*p[1]+n*q[1]!==h[1]||m*p[2]+n*q[2]!==h[2]) throw new Error("빔 "+h+" 이 정대축 ["+zs.z+"] 의 0 차 Laue 영역에 없습니다");
      const g=DC.gvec(phase.cell,h); strong.push({h:h.slice(),m,n,g:norm(g),x:dot(g,zs.e1),y:dot(g,zs.e2)}); }
    strong.sort((a,b)=>a.g-b.g); }
  else { const cand=pts.filter(r=>r.g<=gS*(1+1e-12)); strong=cand.slice(0,shellCut(cand,maxBeams)); }
  const N=strong.length; const gCut=strong[N-1].g; const inS=new Set(strong.map(r=>r.m+","+r.n));
  const rest=pts.filter(r=>!inS.has(r.m+","+r.n));
  const weakB=(explicit?opt.weak===true:opt.weak!==false)?rest.filter(r=>r.g<=gmax*(1+1e-12)):[];
  let bethe=gB>0?rest.filter(r=>r.g<=gB*(1+1e-12)):[]; if(bethe.length>maxBethe) bethe=bethe.slice(0,shellCut(bethe,maxBethe));
  // 들뜸 오차
  const tl=opt.tilt||[0,0]; const kt=[K*(+tl[0]||0)/1000,K*(+tl[1]||0)/1000]; const tilted=!!(kt[0]||kt[1]);
  const sOfK=k=>r=>-(2*(k[0]*r.x+k[1]*r.y)+r.g*r.g)/(2*K);
  const ctx={strong,N,weakB,bethe,Uof,Uk,UT,K,kap,U0,absStep:opt.absStep||0.125};
  const probe=opt.gmProbe||[17,43,89];
  const sOf=sOfK(kt); let res, gmI=null, sol;
  if(!kap){ sol=solveEig(buildH(ctx,sOf),ctx); res=evalEig(ctx,sol,sOf,thick); if(opt.gm!==false){ if(tilted){ const s0=sOfK([0,0]); gmI=evalEig(ctx,solveEig(buildH(ctx,s0),ctx),s0,probe); } else gmI=evalEig(ctx,sol,sOf,probe); } }
  else { if(opt.gm!==false&&!tilted){ const r=stepExp(ctx,buildH(ctx,sOf),sOf,thick.concat(probe)); sol=r.sol; res={I:r.I.map(v=>v.slice(0,thick.length))}; gmI={I:r.I.map(v=>v.slice(thick.length))}; }
    else { const r=stepExp(ctx,buildH(ctx,sOf),sOf,thick); sol=r.sol; res={I:r.I}; if(opt.gm!==false){ const s0=sOfK([0,0]); gmI={I:stepExp(ctx,buildH(ctx,s0),s0,probe).I}; } } }
  // 운동학 |F|² (원점 무관)
  const all=strong.concat(weakB); const gOut=gmax*(1+1e-9); let Ikmax=0; for(const r of all){ if(!r.g) continue; const u=Uof(r.m,r.n); r.Ikin=(u[0]*u[0]+u[1]*u[1])/(Ufac*Ufac); if(r.g<=gOut&&r.Ikin>Ikmax) Ikmax=r.Ikin; }
  const byKey=new Map(all.map((r,i)=>[r.m+","+r.n,i]));
  const beams=[]; let direct=null;
  all.forEach((r,i)=>{ const I=res.I[i]; if(!r.g){ direct={h:[0,0,0],I}; return; } if(r.g>gOut) return; // 출력은 g ≤ gmax 만 (강한 빔은 더 멀리까지 계산)
    const absent=DC.absent(phase.ops,r.h); const kinRel=Ikmax>0?r.Ikin/Ikmax:0; const kinForbidden=absent||kinRel<1e-6; const u=Uof(r.m,r.n); const uA=Math.hypot(u[0],u[1]);
    const b={h:r.h.slice(),g:r.g,x:r.x,y:r.y,I,Ikin:r.Ikin,IkinRel:kinRel,absent,kinForbidden,pert:i>=N,s:sOf(r),U:[u[0],u[1]],xi:uA>0?K/uA:Infinity,ddRoute:[],nRoute:0,gmExtinct:false,dynZero:false};
    if(kinForbidden){ // 이중 회절 경로: 허용 g1 + 허용 g2 = g (출력 빔 집합 안)
      const R=[]; for(const a of all){ if(!a.g||a===r||a.g>gOut) continue; const j=byKey.get((r.m-a.m)+","+(r.n-a.n)); if(j==null) continue; const c=all[j]; if(!c.g||c.g>gOut) continue;
        const ka=a.Ikin/Ikmax, kc=c.Ikin/Ikmax; if(ka<1e-4||kc<1e-4) continue; if(DC.absent(phase.ops,a.h)||DC.absent(phase.ops,c.h)) continue; if(a.m*1e6+a.n>c.m*1e6+c.n) continue; R.push({p:[a.h.slice(),c.h.slice()],w:ka*kc}); }
      R.sort((x,y)=>y.w-x.w); b.ddRoute=R.slice(0,6).map(x=>x.p); b.nRoute=R.length;
      if(gmI){ let ok=true; for(let k=0;k<probe.length;k++){ let mx=0; for(let jj=0;jj<all.length;jj++){ if(!all[jj].g||all[jj].g>gOut) continue; if(gmI.I[jj][k]>mx) mx=gmI.I[jj][k]; } if(gmI.I[i][k]>1e-9*mx){ ok=false; break; } } b.dynZero=ok; b.gmExtinct=ok&&R.length>0; } }
    beams.push(b); });
  return {phase:phase.name||phase.formula||"",zone:zs.z,axis:zs.axis,e1:zs.e1,e2:zs.e2,thickness:thick,direct,beams,
    meta:{kV,lambda:lamA/10,gamma:gam,K,K0,U0,absorption:kap,nBeams:N,nWeak:weakB.length,nBethe:bethe.length,gCut,gBethe:bethe.length?bethe[bethe.length-1].g:0,gmax,tilt:tl.slice(),origin:r0,centro:!!inv,
      solver:sol.mode,dim:sol.dim,orthoErr:sol.orthoErr,ms:now()-T0,msSolve:sol.ms,warn}}; };

/* 껍질(같은 |g|) 을 자르지 않고 앞 n 개 이하 — 대칭 보존 (GM 소멸이 수치적으로도 정확히 0) */
function shellCut(L,n){ if(L.length<=n) return L.length; let k=n; while(k>1&&Math.abs(L[k-1].g-L[k].g)<=1e-9*L[k].g) k--; if(k<=1){ k=n; while(k<L.length&&Math.abs(L[k].g-L[k-1].g)<=1e-9*L[k].g) k++; } return k; }

/* 구조 행렬 H (복소 N×N, nm⁻¹): 대각 s_g + iκU0/2K, 비대각 (1+iκ)U_{g−h}/2K, + Bethe Σ_w H_gw H_wh / (0 − s̃_w) */
function buildH(ctx,sOf){ const {strong,N,bethe,Uof,Uk,UT,K,kap,U0}=ctx; const f=1/(2*K); const Hr=new Float64Array(N*N), Hi=new Float64Array(N*N); const sq=kap*U0*f;
  for(let i=0;i<N;i++) for(let j=0;j<N;j++){ if(i===j){ Hr[i*N+i]=sOf(strong[i]); Hi[i*N+i]=sq; continue; } const u=Uof(strong[i].m-strong[j].m,strong[i].n-strong[j].n); Hr[i*N+j]=(u[0]-kap*u[1])*f; Hi[i*N+j]=(u[1]+kap*u[0])*f; }
  if(bethe.length){ const cr=new Float64Array(N), ci=new Float64Array(N), rr=new Float64Array(N), ri=new Float64Array(N); let realU=!kap;
    if(realU) for(const w of bethe.slice(0,8)) for(let h=0;h<Math.min(N,30);h++){ if(Uof(strong[h].m-w.m,strong[h].n-w.n)[1]!==0){ realU=false; break; } }
    for(const w of bethe){ const sw=sOf(w); const dr=-sw, di=-sq; const dd=dr*dr+di*di; if(dd<1e-30) continue; let mx=0;
      for(let h=0;h<N;h++){ const ka=Uk(strong[h].m-w.m,strong[h].n-w.n), kb=Uk(w.m-strong[h].m,w.n-strong[h].n); const a0=UT[ka], a1=UT[ka+1], b0=UT[kb], b1=UT[kb+1]; cr[h]=(a0-kap*a1)*f; ci[h]=(a1+kap*a0)*f; rr[h]=(b0-kap*b1)*f; ri[h]=(b1+kap*b0)*f; const m2=cr[h]*cr[h]+ci[h]*ci[h]; if(m2>mx) mx=m2; }
      if(mx/Math.sqrt(dd)<1e-14) continue; if(realU) for(let h=0;h<N;h++) if(ci[h]||ri[h]){ realU=false; break; }
      if(!di&&realU){ const iv=1/dr; for(let g=0;g<N;g++){ const xr=cr[g]*iv; if(!xr) continue; const o=g*N; for(let h=0;h<N;h++) Hr[o+h]+=xr*rr[h]; } continue; } // 실수 빠른 길
      for(let g=0;g<N;g++){ const xr=(cr[g]*dr+ci[g]*di)/dd, xi=(ci[g]*dr-cr[g]*di)/dd; if(!xr&&!xi) continue; const o=g*N; for(let h=0;h<N;h++){ Hr[o+h]+=xr*rr[h]-xi*ri[h]; Hi[o+h]+=xr*ri[h]+xi*rr[h]; } } } }
  return {Hr,Hi}; }

/* Hermite H (κ = 0) → Bloch 파. 실수면 N 실대칭, 아니면 2N [[A,−B],[B,A]] */
function solveEig(H,ctx){ const T0=now(); const N=ctx.N; let mR=0, mI=0; for(let i=0;i<N*N;i++){ mR=Math.max(mR,Math.abs(H.Hr[i])); mI=Math.max(mI,Math.abs(H.Hi[i])); }
  const real=mI<=1e-13*mR; const n=real?N:2*N; const A=new Float64Array(n*n);
  for(let i=0;i<N;i++) for(let j=0;j<N;j++){ const re=(H.Hr[i*N+j]+H.Hr[j*N+i])/2, im=(H.Hi[i*N+j]-H.Hi[j*N+i])/2; // Hermite 부분 (반올림 대칭화)
    if(real) A[i*n+j]=re; else { A[i*n+j]=re; A[(i+N)*n+(j+N)]=re; A[i*n+(j+N)]=-im; A[(i+N)*n+j]=im; } }
  const E=symEig(A,n); const lam=E.d, C=E.C;
  let orthoErr=0; for(let a=0;a<Math.min(n,6);a++) for(let b=0;b<Math.min(n,6);b++){ let s=0; const ca=C[a], cb=C[b]; for(let k=0;k<n;k++) s+=ca[k]*cb[k]; orthoErr=Math.max(orthoErr,Math.abs(s-(a===b?1:0))); }
  // 복소 벡터 (행 = 빔): 2N 이면 v = x + i y, 가중 ½ (각 복소 방향이 두 번 나온다: v 와 iv)
  const M=n; const Vr=new Float64Array(N*M), Vi=new Float64Array(N*M);
  for(let k=0;k<M;k++){ const c=C[k]; for(let h=0;h<N;h++){ Vr[h*M+k]=c[h]; Vi[h*M+k]=real?0:c[h+N]; } }
  return {mode:real?"real":"complex",dim:n,M,w:real?1:0.5,Vr,Vi,lam,orthoErr,ms:now()-T0}; }

/* 두께마다 세기 (κ = 0): 강한 빔 = Bloch 합, 바깥 빔 = Bloch 파에 대한 1 차 섭동 */
function evalEig(ctx,S,sOf,thick){ const {strong,N,weakB,Uof,K}=ctx; const M=S.M, w=S.w; const nt=thick.length; const I=new Array(N+weakB.length);
  const i0=strong.findIndex(r=>!r.g); const ar=new Float64Array(M), ai=new Float64Array(M); for(let k=0;k<M;k++){ ar[k]=S.Vr[i0*M+k]; ai[k]=-S.Vi[i0*M+k]; } // a_k = conj(v_k[0])
  const Er=new Float64Array(nt*M), Ei=new Float64Array(nt*M); // w · a_k · exp(2πiγt)
  for(let it=0;it<nt;it++){ const t=thick[it]; for(let k=0;k<M;k++){ const ph=TWO_PI*S.lam[k]*t; const cr=Math.cos(ph)*w, ci=Math.sin(ph)*w; Er[it*M+k]=cr*ar[k]-ci*ai[k]; Ei[it*M+k]=cr*ai[k]+ci*ar[k]; } }
  for(let h=0;h<N;h++){ const out=new Array(nt); const o=h*M; for(let it=0;it<nt;it++){ let re=0, im=0; const e=it*M; for(let k=0;k<M;k++){ const vr=S.Vr[o+k], vi=S.Vi[o+k], cr=Er[e+k], ci=Ei[e+k]; re+=vr*cr-vi*ci; im+=vr*ci+vi*cr; } out[it]=re*re+im*im; } I[h]=out; }
  if(weakB.length){ const f=1/(2*K); const Wr=new Float64Array(M), Wi=new Float64Array(M);
    for(let b=0;b<weakB.length;b++){ const g=weakB[b]; Wr.fill(0); Wi.fill(0);
      for(let h=0;h<N;h++){ const u=Uof(g.m-strong[h].m,g.n-strong[h].n); const a=u[0]*f, c=u[1]*f; if(!a&&!c) continue; const o=h*M; for(let k=0;k<M;k++){ const vr=S.Vr[o+k], vi=S.Vi[o+k]; Wr[k]+=a*vr-c*vi; Wi[k]+=a*vi+c*vr; } }
      const sg=sOf(g); const out=new Array(nt);
      for(let it=0;it<nt;it++){ const t=thick[it]; let re=0, im=0; const esr=Math.cos(TWO_PI*sg*t), esi=Math.sin(TWO_PI*sg*t);
        // Φ_k = (e^{2πiγt} − e^{2πist})/(γ − s) = e^{2πist} · 2πit · (e^x − 1)/x,  x = 2πi(γ − s)t
        for(let k=0;k<M;k++){ const m=cexpm1OverX(0,TWO_PI*(S.lam[k]-sg)*t); let pr=-TWO_PI*t*m[1], pi=TWO_PI*t*m[0]; const qr=pr*esr-pi*esi, qi=pr*esi+pi*esr; pr=qr; pi=qi;
          const war=Wr[k]*ar[k]-Wi[k]*ai[k], wai=Wr[k]*ai[k]+Wi[k]*ar[k]; re+=w*(war*pr-wai*pi); im+=w*(war*pi+wai*pr); }
        out[it]=re*re+im*im; }
      I[N+b]=out; } }
  return {I}; }

/* 복소 행렬 곱 (N×N) */
function cmul(Ar,Ai,Br,Bi,N){ const Cr=new Float64Array(N*N), Ci=new Float64Array(N*N);
  for(let i=0;i<N;i++){ const o=i*N; for(let k=0;k<N;k++){ const ar=Ar[o+k], ai=Ai[o+k]; if(!ar&&!ai) continue; const ok=k*N; for(let j=0;j<N;j++){ const br=Br[ok+j], bi=Bi[ok+j]; Cr[o+j]+=ar*br-ai*bi; Ci[o+j]+=ar*bi+ai*br; } } }
  return [Cr,Ci]; }
/* exp(2πi δ H) — 축척·제곱 + Taylor 12 차 (‖A‖₁ ≤ 0.5 에서 오차 < 1e−14) */
function expmStep(H,N,delta){ const Ar=new Float64Array(N*N), Ai=new Float64Array(N*N); const c=TWO_PI*delta; for(let i=0;i<N*N;i++){ Ar[i]=-c*H.Hi[i]; Ai[i]=c*H.Hr[i]; }
  let nrm=0; for(let j=0;j<N;j++){ let s=0; for(let i=0;i<N;i++) s+=Math.hypot(Ar[i*N+j],Ai[i*N+j]); nrm=Math.max(nrm,s); }
  const sc=Math.max(0,Math.ceil(Math.log2(nrm/0.5))); const f=Math.pow(2,-sc); for(let i=0;i<N*N;i++){ Ar[i]*=f; Ai[i]*=f; }
  let Pr=new Float64Array(N*N), Pi=new Float64Array(N*N); for(let i=0;i<N;i++) Pr[i*N+i]=1;
  for(let k=12;k>=1;k--){ const [Qr,Qi]=cmul(Ar,Ai,Pr,Pi,N); for(let i=0;i<N*N;i++){ Qr[i]/=k; Qi[i]/=k; } for(let i=0;i<N;i++) Qr[i*N+i]+=1; Pr=Qr; Pi=Qi; }
  for(let s=0;s<sc;s++){ const [Qr,Qi]=cmul(Pr,Pi,Pr,Pi,N); Pr=Qr; Pi=Qi; }
  return [Pr,Pi]; }

/* 흡수 (κ > 0): 비 Hermite H 를 고유값 없이 — 두께 간격 δ 의 전파 행렬 exp(2πiδH) 로 한 걸음씩.
   바깥 빔: dψ_w/dz = 2πi(s̃_w ψ_w + Σ_h H_wh ψ_h), 원천을 걸음 사이 선형으로 보고 정확히 적분. 두께는 δ 의 배수로 맞춘다 */
function stepExp(ctx,H,sOf,thick){ const T0=now(); const {strong,N,weakB,Uof,K,kap,U0}=ctx; const d=ctx.absStep; const nt=thick.length;
  const idx=thick.map(t=>Math.max(0,Math.round(t/d))); const nMax=Math.max(0,...idx);
  const [Er,Ei]=expmStep(H,N,d); const msExp=now()-T0;
  const W=weakB.length; const f=1/(2*K); const sq=kap*U0*f;
  const Wr=new Float64Array(W*N), Wi=new Float64Array(W*N); const e1r=new Float64Array(W), e1i=new Float64Array(W), p1r=new Float64Array(W), p1i=new Float64Array(W), p2r=new Float64Array(W), p2i=new Float64Array(W);
  for(let b=0;b<W;b++){ const g=weakB[b]; for(let h=0;h<N;h++){ const u=Uof(g.m-strong[h].m,g.n-strong[h].n); // 2πi H_c[w,h]
      const hr=(u[0]-kap*u[1])*f, hi=(u[1]+kap*u[0])*f; Wr[b*N+h]=-TWO_PI*hi; Wi[b*N+h]=TWO_PI*hr; }
    // a = 2πi s̃ = 2πi(s + i sq) → aδ;  φ1 = (e^{aδ}−1)/a, φ2 = (e^{aδ} − 1 − aδ)/(a²δ)
    const ar=-TWO_PI*sq, ai=TWO_PI*sOf(g); const xr=ar*d, xi=ai*d; const ex=Math.exp(xr); e1r[b]=ex*Math.cos(xi); e1i[b]=ex*Math.sin(xi);
    const m1=cexpm1OverX(xr,xi); p1r[b]=m1[0]*d; p1i[b]=m1[1]*d; // φ1 = δ·(e^x−1)/x
    // φ2 = δ·(e^x − 1 − x)/x² = δ·(m1 − 1)/x
    const ur=m1[0]-1, ui=m1[1]; const mm=xr*xr+xi*xi; if(mm<1e-8){ p2r[b]=d*(0.5+xr/6); p2i[b]=d*(xi/6); } else { p2r[b]=d*(ur*xr+ui*xi)/mm; p2i[b]=d*(ui*xr-ur*xi)/mm; } }
  const I=new Array(N+W); for(let i=0;i<N+W;i++) I[i]=new Array(nt);
  let pr=new Float64Array(N), pi=new Float64Array(N); const i0=strong.findIndex(r=>!r.g); pr[i0]=1;
  const wr=new Float64Array(W), wi=new Float64Array(W); let Sr=new Float64Array(W), Si=new Float64Array(W);
  const src=(Sr,Si)=>{ for(let b=0;b<W;b++){ let re=0, im=0; const o=b*N; for(let h=0;h<N;h++){ const a=Wr[o+h], c=Wi[o+h]; re+=a*pr[h]-c*pi[h]; im+=a*pi[h]+c*pr[h]; } Sr[b]=re; Si[b]=im; } };
  if(W) src(Sr,Si);
  const rec=n=>{ for(let it=0;it<nt;it++) if(idx[it]===n){ for(let h=0;h<N;h++) I[h][it]=pr[h]*pr[h]+pi[h]*pi[h]; for(let b=0;b<W;b++) I[N+b][it]=wr[b]*wr[b]+wi[b]*wi[b]; } };
  rec(0);
  const nr=new Float64Array(N), ni=new Float64Array(N); const S2r=new Float64Array(W), S2i=new Float64Array(W);
  for(let n=1;n<=nMax;n++){ nr.fill(0); ni.fill(0); for(let i=0;i<N;i++){ let re=0, im=0; const o=i*N; for(let j=0;j<N;j++){ const a=Er[o+j], c=Ei[o+j]; re+=a*pr[j]-c*pi[j]; im+=a*pi[j]+c*pr[j]; } nr[i]=re; ni[i]=im; }
    pr.set(nr); pi.set(ni);
    if(W){ src(S2r,S2i); for(let b=0;b<W;b++){ const dSr=S2r[b]-Sr[b], dSi=S2i[b]-Si[b];
        const r1=e1r[b]*wr[b]-e1i[b]*wi[b]+(Sr[b]*p1r[b]-Si[b]*p1i[b])+(dSr*p2r[b]-dSi*p2i[b]);
        const i1=e1r[b]*wi[b]+e1i[b]*wr[b]+(Sr[b]*p1i[b]+Si[b]*p1r[b])+(dSr*p2i[b]+dSi*p2r[b]); wr[b]=r1; wi[b]=i1; }
      Sr.set(S2r); Si.set(S2i); }
    rec(n); }
  return {I,sol:{mode:"expm",dim:N,orthoErr:0,ms:now()-T0,msExp,step:d}}; }

/* 두께 평균: tRange [t1,t2] 안의 계산 두께 (사다리꼴 가중). 반환 {I:[빔마다], I0, n} */
DD.thicknessAverage=function(res,tRange){ const T=res.thickness; const lo=tRange?tRange[0]:-Infinity, hi=tRange?tRange[1]:Infinity; const idx=[]; T.forEach((t,i)=>{ if(t>=lo-1e-9&&t<=hi+1e-9) idx.push(i); });
  if(!idx.length) throw new Error("두께 범위 안 계산 점이 없습니다"); const w=new Array(idx.length).fill(1);
  if(idx.length>1){ for(let k=0;k<idx.length;k++){ const a=k>0?T[idx[k]]-T[idx[k-1]]:0, b=k<idx.length-1?T[idx[k+1]]-T[idx[k]]:0; w[k]=(a+b)/2; } }
  const W=w.reduce((s,x)=>s+x,0); const av=I=>idx.reduce((s,i,k)=>s+w[k]*I[i],0)/W;
  return {I:res.beams.map(b=>av(b.I)),I0:res.direct?av(res.direct.I):null,n:idx.length,tRange:[T[idx[0]],T[idx[idx.length-1]]]}; };

/* 두께 모를 때: 10–80 nm 균일 평균 (2 nm 간격) + 선택적 미세 기울임 평균 (tiltAvg mrad 반지름 원 위 nTilt 개 + 중심).
   결과 beams 에 Iavg (입사 1 대비), Irel (회절 빔 최대 = 1). 같은 상·정대축·조건은 캐시 */
const CACHE=new WeakMap();
DD.zoneAverage=function(phase,uvw,opt){ opt=Object.assign({},opt||{}); const tR=opt.tRange||[10,80], step=opt.tStep||2; const th=[]; for(let t=tR[0];t<=tR[1]+1e-9;t+=step) th.push(+t.toFixed(6));
  const key=JSON.stringify([DC.reduceInt(uvw),opt.kV||200,opt.gmax||12,opt.maxBeams||150,opt.gStrong||0,opt.gBethe,opt.maxBethe||0,opt.absorption||0,tR,step,opt.tiltAvg||0,opt.nTilt||6,opt.weak,opt.dw,opt.Bdefault==null?0.5:opt.Bdefault,opt.Bset,opt.ref||null]);
  let m=CACHE.get(phase); if(!m){ m=new Map(); CACHE.set(phase,m); } if(m.has(key)) return m.get(key);
  const T0=now(); const base=DD.zoneIntensities(phase,uvw,Object.assign({},opt,{thickness:th,tilt:[0,0]})); let av=DD.thicknessAverage(base).I; let nT=1;
  if(opt.tiltAvg>0){ const nt=opt.nTilt||6; const sum=av.slice();
    for(let k=0;k<nt;k++){ const a=TWO_PI*k/nt; const r=DD.zoneIntensities(phase,uvw,Object.assign({},opt,{thickness:th,tilt:[opt.tiltAvg*Math.cos(a),opt.tiltAvg*Math.sin(a)],gm:false}));
      const A=DD.thicknessAverage(r).I; for(let i=0;i<sum.length;i++) sum[i]+=A[i]; nT++; }
    av=sum.map(x=>x/nT); }
  let mx=0; for(const x of av) if(x>mx) mx=x; base.beams.forEach((b,i)=>{ b.Iavg=av[i]; b.Irel=mx>0?av[i]/mx:0; });
  base.meta.avg={tRange:tR,tStep:step,tiltAvg:opt.tiltAvg||0,nTilt:nT,ms:now()-T0}; m.set(key,base); return base; };

/* hkl → 빔 찾기 (Friedel 짝·같은 hkl) */
DD.beamOf=(res,h)=>res.beams.find(b=>b.h[0]===h[0]&&b.h[1]===h[1]&&b.h[2]===h[2])||null;

if(typeof module!=="undefined"&&module.exports) module.exports=DD;
G.DIFF_DYN=DD;
})();
