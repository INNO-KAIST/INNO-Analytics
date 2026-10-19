/* ============================================================
   INNO Analytics — Diffraction: 지수 붙이기·상 찾기·보정. DOM 없음.
   · matchSpots: 잰 격자(R1·R2·각) ↔ 상마다 정대축 표 → FOM 표. 보정 없이(모양: 길이비·각) 또는 보정 있게(길이도)
   · verifyZone: 고른 정대축으로 모든 스폿 예측 ↔ 관측 (빠진 강한 스폿, 운동학적 금지 자리 = 이중 회절 가능)
   · matchRings: 링 봉우리 ↔ 상마다 링 표, 보정 없이 κ 찾기
   · calibrateAffine / calibrateRings: 알려진 표준으로 화소 → 1/nm (평균 축척 + 타원 찌그러짐)
   길이: 관측 px, 보정 s (1/nm/px), 예측 g (1/nm). κ = g_예측 / (R·s) — 1 이면 보정이 맞음.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const DC=G.DIFF_CORE||(typeof require!=="undefined"?require("./diff_core.js"):null);
const DI={};
const D2R=Math.PI/180;

/* 잰 값 꾸러미: fit = DIFF_PATTERN.fitLattice 결과 (len·sig px, ang12, sigAng) 또는 손으로 {len:[R1,R2,R3], ang12, sig?, sigAng?} */
DI.measured=function(fit,opt){ opt=opt||{}; const L=fit.len, S=fit.sig||[0,0,0]; const sys=opt.distortion!=null?opt.distortion:0.01; // 렌즈 찌그러짐·보정 잔차 (상대)
  const r21=L[1]/L[0]; const sr=Math.hypot(r21*Math.hypot(S[0]/L[0],S[1]/L[1]),r21*sys*Math.SQRT2);
  const sa=Math.hypot(fit.sigAng||0,sys*180/Math.PI); // 찌그러짐 e 는 각을 최대 ~e rad 바꾼다
  return {R:L.slice(0,3),sigR:S.slice(0,3),r21,sr21:sr,ang:fit.ang12,sAng:sa,sys}; };

/* 정대축 표를 상마다 한 번만 */
const TABLE=new WeakMap();
DI.zoneTable=function(phase,opt){ opt=opt||{}; const k=(opt.maxIdx||5)+":"+(opt.gmax||12); let m=TABLE.get(phase); if(!m){ m={}; TABLE.set(phase,m); } if(!m[k]) m[k]=DC.zoneTable(phase,Object.assign({},opt,{refl:DI.refl3(phase,opt.gmax||12)})); return m[k]; };
/* 3D 반사(+|F|²) 를 상·gmax 마다 한 번만 — 정대축 표와 verifyZone 의 세기 기준(최댓값)이 같이 쓴다. gmax 는 GCAP 으로 막는다 */
const GCAP=25; // 1/nm (d ≥ 0.04 nm) — SAED 에서 쓸 만한 끝. 이보다 큰 반사 목록은 만들지 않는다 (가짜 격자 → 10⁹ 개 반사 → 탭 메모리 초과)
const R3=new WeakMap();
DI.refl3=function(phase,gmax){ const G=Math.min(GCAP,gmax||12); let m=R3.get(phase); if(!m){ m={}; R3.set(phase,m); } if(!m[G]){ const L=DC.reflections(phase,G,{F:!!(phase.atoms&&phase.atoms.length)}); const pm=new Float64Array(L.length); let mx=0; for(let i=0;i<L.length;i++){ if(L[i].I!=null&&L[i].I>mx) mx=L[i].I; pm[i]=mx; } L.prefMax=pm; m[G]=L; } return m[G]; };
/* g ≤ gmax 인 3D 반사 중 |F|² 최댓값 (refl3 의 누적 최댓값으로) */
function imaxUpTo(phase,g,G){ const L=DI.refl3(phase,G); let lo=0, hi=L.length; while(lo<hi){ const mid=(lo+hi)>>1; if(L[mid].g<=g) lo=mid+1; else hi=mid; } return lo?L.prefMax[lo-1]:0; }
/* 정대축 [uvw] 의 0 차 Laue 면 정수 격자점을 2D 로만 늘어놓는다 (|g| ≤ gmax, 원점 제외) — 3D 목록에서 거르는 것보다 수백 배 적다.
   순서는 DC.reflections 와 같다 (g 오름차순, 같은 g 는 h 내림차순). maxN 을 넘으면 null */
DI.zoneLattice=function(phase,uvw,gmax,maxN){ const cell=phase.cell; const z=DC.reduceInt(uvw); const pb=DC.planeBasis(z); if(!pb) return []; const [p,q]=pb;
  const gp=DC.gvec(cell,p), gq=DC.gvec(cell,q); const a=gp[0]*gp[0]+gp[1]*gp[1]+gp[2]*gp[2], b=gp[0]*gq[0]+gp[1]*gq[1]+gp[2]*gq[2], c=gq[0]*gq[0]+gq[1]*gq[1]+gq[2]*gq[2]; const det=a*c-b*b; if(!(det>0)) return [];
  const N1=Math.ceil(gmax*Math.sqrt(c/det))+1, N2=Math.ceil(gmax*Math.sqrt(a/det))+1; if((2*N1+1)*(2*N2+1)>(maxN||40000)*1.6+64) return null;
  const out=[]; for(let i=-N1;i<=N1;i++) for(let j=-N2;j<=N2;j++){ if(!i&&!j) continue; const h=[i*p[0]+j*q[0],i*p[1]+j*q[1],i*p[2]+j*q[2]]; const g=DC.gLen(cell,h); if(g>gmax) continue; out.push({h,g,d:1/g}); }
  if(maxN&&out.length>maxN) return null;
  out.sort((x,y)=>x.g-y.g||(y.h[0]-x.h[0])||(y.h[1]-x.h[1])||(y.h[2]-x.h[2])); return out; };

/* matchSpots(meas, phases, {scale, scaleErr, maxIdx, gmax, top})
   scale 이 있으면 R1·s 와 g1 도 비교 (scaleErr 상대 오차 = REDI "측정 오차 %"). 없으면 모양만.
   R1≈R2 일 때 순서가 바뀔 수 있어 두 순서 모두 본다. 결과: FOM 오름차순 (χ²/자유도) */
DI.matchSpots=function(meas,phases,opt){ opt=opt||{}; const out=[]; const s=opt.scale||null, se=opt.scaleErr!=null?opt.scaleErr:0.03;
  for(const ph of phases){ const T=DI.zoneTable(ph,opt);
    for(const z of T){ for(const swap of [false,true]){
        if(swap&&Math.abs(z.ratio21-1)>3*meas.sr21+0.05) continue; // 거의 같은 길이일 때만 뒤집어 본다
        const c1=(meas.r21-(swap?1/z.ratio21:z.ratio21))/meas.sr21; const c2=(meas.ang-z.ang12)/meas.sAng; let chi=c1*c1+c2*c2, dof=2;
        // R1 에 맞는 예측 g: 뒤집으면 R1 ↔ g2
        const g1=swap?z.g[1]:z.g[0], g2=swap?z.g[0]:z.g[1];
        const kappa=s?((g1/(meas.R[0]*s))+(g2/(meas.R[1]*s)))/2:null; const scaleImplied=(g1/meas.R[0]+g2/meas.R[1])/2; // 1/nm/px
        if(s){ const c3=(kappa-1)/Math.hypot(se,meas.sys,meas.sigR[0]/meas.R[0]); chi+=c3*c3; dof++; }
        const h1=swap?z.h2:z.h1, h2=swap?z.h1:z.h2; const zz=DC.zoneAxis(h1,h2); // 정대축 = h1 × h2 (R1·R2 순서와 같은 손잡이)
        out.push({phase:ph,phaseName:ph.name||ph.formula||"?",zone:zz,h1,h2,h3:z.h3,g:[g1,g2,z.g[2]],d:[1/g1,1/g2,1/z.g[2]],ang12:z.ang12,ratio21:g2/g1,chi2:chi,fom:chi/dof,kappa,scaleImplied,swap}); } } }
  out.sort((a,b)=>a.fom-b.fom||Math.max(...a.zone.map(Math.abs))-Math.max(...b.zone.map(Math.abs)));
  // 같은 상·같은 정대축(뒤집은 쌍) 은 좋은 것 하나만
  const seen=new Set(), res=[]; for(const r of out){ const zk=r.zone.join(",")>r.zone.map(x=>-x).join(",")?r.zone.join(","):r.zone.map(x=>-x).join(","); const k=(r.phase.id||r.phaseName)+"|"+zk; /* [uvw] 와 [-u-v-w] 는 같은 정대축 */ if(seen.has(k)) continue; seen.add(k); res.push(r); }
  return res.slice(0,opt.top||50); };

/* verifyZone(fit, cand, {gmax, scale?, rObsMax?, mask?, w?, h?})
   관측 기저 (R1,R2 px 벡터) 를 cand.h1·h2 에 맞추면 모든 관측 스폿 (n1,n2) → hkl = n1 h1 + n2 h2.
   예측: 정대축의 허용 반사(|F|²) 가 관측 범위 안에 있는데 관측이 없는 것 = 빠짐 (강도 순).  관측했는데 |F|²≈0 인 자리 = 이중 회절 가능 */
DI.verifyZone=function(fit,cand,opt){ opt=opt||{}; const ph=cand.phase; const cell=ph.cell; const hk=(n)=>[0,1,2].map(i=>n[0]*cand.h1[i]+n[1]*cand.h2[i]);
  const obs=fit.spots.filter(sp=>sp.indexed); const rmax=obs.length?Math.max(...obs.map(sp=>Math.hypot(sp.n[0]*fit.R1[0]+sp.n[1]*fit.R2[0],sp.n[0]*fit.R1[1]+sp.n[1]*fit.R2[1]))):0;
  const g1=DC.gLen(cell,cand.h1); const pxPerG=fit.len[0]/g1;
  // 예측 범위: 관측 반지름 → 1/nm. opt.gmax(정대축 표, 기본 15) 와 GCAP(25) 을 넘지 않게 — 가짜(얼룩) 격자에서 787 1/nm 같은 값이 나와 반사 수가 폭발했다
  const gCap=Math.min(GCAP,Math.max(opt.gmax||15,20)); const gNeed=(opt.rObsMax||rmax)/pxPerG*1.0001; const gmax=Math.min(gNeed,gCap);
  const Im=imaxUpTo(ph,gmax,Math.min(gCap,opt.gmax||15)); const hasF=!!(ph.atoms&&ph.atoms.length);
  const ZL=DI.zoneLattice(ph,cand.zone,gmax,opt.maxZonePts||20000); const inZone=[];
  if(ZL&&Im>0&&hasF) for(const r of ZL){ r.I=DC.structureFactor(ph,r.h).I; inZone.push(r); }
  const got=new Map(); for(const sp of obs){ const h=hk(sp.n); got.set(h.join(","),sp); }
  // 관측 자리의 예측 세기
  const rows=obs.map(sp=>{ const h=hk(sp.n); const F=Im>0?DC.structureFactor(ph,h).I:0; const absent=DC.absent(ph.ops,h); return {spot:sp,h,d:DC.dSpacing(cell,h),I:Im>0?F/Im:null,absent,weak:Im>0?F<1e-4*Im:absent}; });
  // 예측 강한 반사인데 관측 없음 — 관측 범위(반지름) 안, 가려진 곳 제외 (opt.posOf(h) → {x,y} 가 있으면 mask 확인)
  const missing=[]; if(Im>0) for(const r of inZone){ if(r.I<0.02*Im) continue; if(got.has(r.h.join(","))) continue; if(opt.isMasked&&opt.isMasked(r.h)) continue; missing.push({h:r.h,d:r.d,I:r.I/Im}); }
  missing.sort((a,b)=>b.I-a.I);
  const strongPred=Im>0?inZone.filter(r=>r.I>=0.02*Im).length:0; const ddN=rows.filter(r=>r.weak||r.absent).length;
  return {rows,missing,nObs:obs.length,nStrongPred:strongPred,completeness:strongPred?1-missing.length/strongPred:1,doubleDiff:ddN,pxPerG,gmax,clipped:gNeed>gCap,tooMany:!ZL}; };

/* 두 단계 순위: (1) 기하 FOM (matchSpots) 상위 nVerify 개에 (2) 세기 일치 — 관측 스폿 세기 ↔ 운동학 |F|² 의 Spearman ρ,
   운동학적으로 거의 0 인 자리에 놓인 관측 세기 비율 forb (이중 회절로 생길 수 있어 0 이 아닐 수 있다).
   score = FOM + 5·forb + (1 − ρ)  (작을수록 좋음). 운동학 세기는 영역축 동역학 회절에서 순서만 대략 맞으므로 ρ 는 순위만 쓴다 */
DI.spearman=function(a,b){ const n=a.length; if(n<3) return 0; const rk=v=>{ const o=v.map((x,i)=>[x,i]).sort((p,q)=>p[0]-q[0]); const r=new Array(n); let i=0; while(i<n){ let j=i; while(j+1<n&&o[j+1][0]===o[i][0]) j++; const m=(i+j)/2; for(let k=i;k<=j;k++) r[o[k][1]]=m; i=j+1; } return r; };
  const ra=rk(a), rb=rk(b); const ma=(n-1)/2; let sab=0,saa=0,sbb=0; for(let i=0;i<n;i++){ sab+=(ra[i]-ma)*(rb[i]-ma); saa+=(ra[i]-ma)**2; sbb+=(rb[i]-ma)**2; } return saa&&sbb?sab/Math.sqrt(saa*sbb):0; };
/* 세기 순위 상관: 포화 스폿은 세기 하한이라 뺀다 (4 개 미만이면 null → 점수에서 중립 0.5) */
DI.rhoUnsat=function(rows,Io,Ip){ const a=[], b=[]; rows.forEach((r,i)=>{ if(r.spot&&r.spot.saturated) return; a.push(Io[i]); b.push(Ip[i]); }); return a.length>=4?DI.spearman(a,b):null; };
const rhoTerm=r=>r==null?0.5:(1-r);
/* 얼룩(speckle) 격자인가: 기저가 이미지의 2.5 % 보다 짧거나, 지수 붙은 스폿이 기저의 25 배 넘게 퍼짐(링·결정립 얼룩을 가짜 격자로 맞춘 것) */
DI.isSpeckle=function(fit,size){ if(!fit||!fit.len) return false; const m=Math.min(fit.len[0],fit.len[1]); if(size&&m<0.025*size) return true;
  const sp=(fit.spots||[]).filter(s=>s.indexed&&s.n); if(sp.length){ let r=0; for(const s of sp) r=Math.max(r,Math.hypot(s.n[0]*fit.R1[0]+s.n[1]*fit.R2[0],s.n[0]*fit.R1[1]+s.n[1]*fit.R2[1])); if(r>25*m) return true; }
  return (fit.nIndexed||0)>300; };
/* rankSpots(fit, phases, opt)
   opt.scalePrior {scale, tol}: 파일 태그·스케일바처럼 대략 믿을 만한 보정 — κ 가 ±tol 밖인 후보를 뺀다 (모두 빠지면 빼지 않고 priorEmpty)
   opt.budgetMs: 시간 상한 (기본 25 s) — 넘으면 남은 후보는 검증하지 않고 버림 (truncated). opt.onProgress(0..1, 단계)
   얼룩 격자(DI.isSpeckle)는 순위를 매기지 않는다 → {refused:"speckle"} */
DI.rankSpots=function(fit,phases,opt){ opt=opt||{}; const meas=DI.measured(fit,opt); const t0=Date.now(); const budget=opt.budgetMs!=null?opt.budgetMs:25000; const prog=typeof opt.onProgress==="function"?opt.onProgress:null; let lastP=0;
  const tick=(fr,st)=>{ if(!prog) return; const t=Date.now(); if(t-lastP<120&&fr<1) return; lastP=t; try{ prog(fr,st); }catch(e){} };
  if(DI.isSpeckle(fit,Math.max(opt.w||0,opt.h||0))) return {meas,cands:[],refused:"speckle"};
  const nV=opt.nVerify||80; const SP=opt.scalePrior&&opt.scalePrior.scale>0?opt.scalePrior:null;
  let C=DI.matchSpots(meas,phases,Object.assign({},opt,{top:SP?100000:nV})); let priorEmpty=false;
  if(SP){ const tol=SP.tol||0.15; const Ck=C.filter(c=>Math.abs(c.scaleImplied/SP.scale-1)<=tol); if(Ck.length) C=Ck; else priorEmpty=true; C=C.slice(0,nV); }
  // 예측 위치가 이미지 밖·가린 곳이면 '빠짐' 으로 세지 않는다
  const W=opt.w, H=opt.h, M=opt.mask; const PF=opt.posFit||fit; const mk=c=>(W&&H)?(h=>{ const q=DI.predictPx(PF,c,h); const x=Math.round(q.x), y=Math.round(q.y); if(x<4||y<4||x>=W-4||y>=H-4) return true; return !!(M&&M[y*W+x]); }):null;
  let truncated=false; const done=[];
  for(let i=0;i<C.length;i++){ const c=C[i]; if(Date.now()-t0>budget){ truncated=true; break; } tick(0.7*i/C.length,"verify");
    const V=DI.verifyZone(fit,c,Object.assign({},opt,{isMasked:mk(c)})); const Io=V.rows.map(r=>Math.max(0,r.spot.I)), Ip=V.rows.map(r=>r.I); const tot=Io.reduce((a,b)=>a+b,0)||1; let forb=0; V.rows.forEach((r,i)=>{ if(r.I<0.01) forb+=Io[i]; });
    c.verify=V; c.hasF=!!(c.phase.atoms&&c.phase.atoms.length);
    // 빠진 강한 반사: 관측 범위 안 예측 세기 중 관측 안 된 비율 (세기 가중)
    const zI=V.missing.reduce((a,m)=>a+m.I,0), oI=V.rows.reduce((a,r)=>a+(r.I>=0.02?r.I:0),0); c.missFrac=zI+oI>0?zI/(zI+oI):0;
    if(c.hasF){ c.forb=forb/tot; c.rho=DI.rhoUnsat(V.rows,Io,Ip); c.score=c.fom+5*c.forb+rhoTerm(c.rho)+3*c.missFrac; } else { c.forb=null; c.rho=null; c.missFrac=null; c.score=c.fom+1; } done.push(c); }
  C=done; C.sort((a,b)=>a.score-b.score);
  // 동역학(Bloch-wave, ZOLZ, 두께 10–80 nm 평균) 세기로 상위 후보 다시 채점 — 이중 회절·Gjønnes–Moodie 소멸 반영 (docs-dyn.md §6)
  const DY=(typeof DIFF_DYN!=="undefined")?DIFF_DYN:(G.DIFF_DYN||null); const nDyn=opt.dynamic===false||!DY?0:(opt.nDyn||12);
  let nd=0; for(const c of C.slice(0,nDyn)){ if(Date.now()-t0>budget){ truncated=true; break; } tick(0.7+0.3*(nd++)/nDyn,"dynamic"); if(!c.hasF) continue; const V=c.verify; if(!V.rows.length) continue;
    let D=null; try{ const gObs=Math.min(60,Math.max(...V.rows.map(r=>1/r.d))*1.05); /* 동역학 빔 목록은 2D·200 빔 상한이라 60 1/nm 까지는 가볍다 (얼룩 격자는 위에서 거름) */ D=DY.zoneAverage(c.phase,c.zone,{kV:opt.kV||200,gmax:gObs,tRange:opt.tRange||[10,80]}); }catch(e){ c.dynError=String(e.message||e); continue; }
    const Io=V.rows.map(r=>Math.max(0,r.spot.I)), tot=Io.reduce((a,b)=>a+b,0)||1; const Id=[]; let forb=0;
    V.rows.forEach((r,i)=>{ const b=DY.beamOf(D,r.h); Id.push(b?b.Irel:0); const w=!b?1:(b.dynZero&&!b.gmExtinct)?1:b.gmExtinct?0.3:0; forb+=w*Io[i];
      r.Id=b?b.Irel:0; r.dd=!!(b&&b.kinForbidden&&b.nRoute>0&&!b.gmExtinct); r.gm=!!(b&&b.gmExtinct); r.dz=!b||(!!b.dynZero&&!b.gmExtinct); });
    forb/=tot; const rho=DI.rhoUnsat(V.rows,Io,Id); const got=new Set(V.rows.map(r=>r.h.join(","))); const isM=mk(c); let zI=0,oI=0; const dmiss=[];
    for(const b of D.beams){ if(b.Irel<0.02) continue; if(isM&&isM(b.h)) continue; if(got.has(b.h.join(","))) oI+=b.Irel; else { zI+=b.Irel; dmiss.push({h:b.h,d:1/b.g,I:b.Irel}); } }
    const miss=zI+oI>0?zI/(zI+oI):0; dmiss.sort((p,q)=>q.I-p.I);
    c.kin={rho:c.rho,forb:c.forb,miss:c.missFrac,score:c.score}; c.dyn={rho,forb,miss,missing:dmiss.slice(0,12),nBeams:D.meta.nBeams,kV:opt.kV||200,tRange:opt.tRange||[10,80]};
    c.rho=rho; c.forb=forb; c.missFrac=miss; c.score=c.fom+5*forb+rhoTerm(rho)+3*miss; }
  if(nDyn){ const head=C.slice(0,nDyn).sort((a,b)=>a.score-b.score); C.splice(0,nDyn,...head); }
  tick(1,"done");
  return {meas,cands:C,truncated,priorEmpty,ms:Date.now()-t0}; };

/* 화소 위치 예측: 관측 기저(px) 로 hkl → (x,y). h = n1 h1 + n2 h2 를 풀어 n 을 얻는다 */
DI.predictPx=function(fit,cand,h){ const a=cand.h1, b=cand.h2; // 최소제곱 (정수 해가 있으면 정확)
  const aa=a[0]*a[0]+a[1]*a[1]+a[2]*a[2], bb=b[0]*b[0]+b[1]*b[1]+b[2]*b[2], ab=a[0]*b[0]+a[1]*b[1]+a[2]*b[2], ha=h[0]*a[0]+h[1]*a[1]+h[2]*a[2], hb=h[0]*b[0]+h[1]*b[1]+h[2]*b[2]; const det=aa*bb-ab*ab;
  const n1=(ha*bb-hb*ab)/det, n2=(hb*aa-ha*ab)/det; return {x:fit.center.x+n1*fit.R1[0]+n2*fit.R2[0],y:fit.center.y+n1*fit.R1[1]+n2*fit.R2[1],n:[n1,n2]}; };

/* matchRings(peaks, phases, {scale, scaleErr, tol, gmax}) — peaks: [{r(px), h, fwhm}]. 반환: 상마다 최선의 κ·맞춤
   보정 없이: 가장 센 관측 봉우리를 예측 링 j (센 것 8 개) 에 맞춰 κ0 → 모든 봉우리 짝짓기 → 최소제곱 κ → 점수
   점수 = 설명한 관측 세기 비율 × (1 − 관측 범위 안 강한 예측 링 중 빠진 세기 비율)  · 평균 상대 어긋남이 작을수록 좋음 */
DI.matchRings=function(peaks,phases,opt){ opt=opt||{}; const tol=opt.tol||0.012; const P=peaks.filter(p=>p.r>0&&!p.broad&&!p.nearStop); if(!P.length) return [];
  // 높이로 무게 (겹친 봉우리의 반높이 폭은 불안정)
  const area=p=>Math.max(p.h,0); const Iobs=P.reduce((a,p)=>a+area(p),0); const strongest=P.slice().sort((a,b)=>area(b)-area(a)).slice(0,3);
  const rmax=Math.max(...P.map(p=>p.r)), rmin=Math.min(...P.map(p=>p.r)); const out=[];
  const GL=opt.gmaxLimit||20; // 1/nm — SAED 에서 쓸 만한 범위 (d ≥ 0.05 nm). 이보다 큰 축척을 요구하는 짝짓기는 버린다
  for(const ph of phases){ let best=null; const all=DI.ringCache(ph,GL);
    const SP=opt.scalePrior&&opt.scalePrior.scale>0&&!opt.noPrior?opt.scalePrior:null; // 대략 믿는 보정 (태그·스케일바): κ 가 ±tol 밖인 짝짓기는 버린다
    const tryK=(k0)=>{ if(!(k0>0)||rmax*k0>GL) return; if(SP&&Math.abs(k0/SP.scale-1)>(SP.tol||0.15)+0.03) return; const gmax=rmax*k0*1.03; const rings=all.filter(r=>r.g<=gmax); if(!rings.length) return; let Ip=0; for(const r of rings) if(r.I!=null) Ip=Math.max(Ip,r.I);
      let k=k0; let pairs=[];
      for(let it=0;it<4;it++){ pairs=[]; for(const p of P){ const g=p.r*k; let b=null; for(const r of rings){ const e=Math.abs(r.g-g)/r.g; if(e<=tol&&(!b||e<b.e)) b={r,e}; } if(b) pairs.push({p,ring:b.r,e:b.e}); }
        if(!pairs.length) return; let num=0,den=0; for(const q of pairs){ const w=area(q.p); num+=w*q.ring.g*q.p.r; den+=w*q.p.r*q.p.r; } k=num/den; }
      const expl=pairs.reduce((a,q)=>a+area(q.p),0)/Iobs;
      // 관측 범위 안 강한 예측 링이 빠졌는가 (예측 세기 m|F|² 의 5 % 이상)
      let miss=0, tot=0; const used=new Set(pairs.map(q=>q.ring)); for(const r of rings){ const rp=r.g/k; if(rp<rmin*0.98||rp>rmax*1.02) continue; const I=r.I!=null?r.I/Ip:1; if(I<0.05) continue; tot+=I; if(!used.has(r)) miss+=I; }
      const dev=Math.sqrt(pairs.reduce((a,q)=>a+((q.ring.g-q.p.r*k)/q.ring.g)**2,0)/pairs.length);
      // 세기 순서: 방위 평균 높이 ∝ m|F|²/g² (링 둘레 ∝ g, Lorentz ∝ d) — 짝지은 봉우리의 Spearman ρ
      const rho=pairs.length>=3&&pairs.every(q=>q.ring.I!=null)?DI.spearman(pairs.map(q=>area(q.p)),pairs.map(q=>q.ring.I/(q.ring.g*q.ring.g))):0;
      const score=expl*(tot?1-miss/tot:1)*Math.exp(-dev/tol)*(1+rho)/2;
      const kErr=dev/Math.sqrt(Math.max(1,pairs.length-1));
      const c={phase:ph,phaseName:ph.name||ph.formula||"?",rho,scale:k,scaleErr:k*kErr,kappa:opt.scale?k/opt.scale:null,explained:expl,missing:tot?miss/tot:0,dev,score,pairs:pairs.map(q=>({r:q.p.r,g:q.p.r*k,d:1/(q.p.r*k),hkl:q.ring.hkls[0],m:q.ring.m,gPred:q.ring.g,e:(q.p.r*k-q.ring.g)/q.ring.g})),rings};
      if(SP&&Math.abs(k/SP.scale-1)>(SP.tol||0.15)) return;
      if(opt.scale){ const se=opt.scaleErr!=null?opt.scaleErr:0.03; c.score*=Math.exp(-0.5*((k/opt.scale-1)/se)**2); }
      if(!best||c.score>best.score) best=c; };
    if(opt.scale&&opt.fixScale){ tryK(opt.scale); }
    else { const rings0=all; const strongR=rings0.filter(r=>r.I!=null).sort((a,b)=>b.I-a.I).slice(0,8).concat(rings0.slice(0,4));
      for(const sp of strongest) for(const r of strongR) tryK(r.g/sp.r);
      if(opt.scale) tryK(opt.scale); if(SP) tryK(SP.scale); }
    if(best) out.push(best); }
  out.sort((a,b)=>b.score-a.score);
  if(opt.scalePrior&&!opt.noPrior&&!out.length&&phases.length){ const R=DI.matchRings(peaks,phases,Object.assign({},opt,{noPrior:true})); R.priorEmpty=true; return R; } // 창 안에 맞는 상이 없으면 창 없이 (보정이 틀렸을 수 있음)
  return out; };

const RC=new WeakMap(); DI.ringCache=function(ph,gl){ let m=RC.get(ph); if(!m){ m={}; RC.set(ph,m); } if(!m[gl]) m[gl]=DC.ringList(ph,gl); return m[gl]; };

/* 아핀 보정: 관측 R1·R2 (px 벡터) ↔ 예측 g1·g2 (1/nm, 사이각) → T (px → 1/nm 2D). T = Q·S (극분해), S 의 고윳값 = 장·단 방향 축척.
   반환 {scale (기하 평균 1/nm/px), smax, smin, ellip = smax/smin, phi (smin 방향… 화면 각, rad), S} */
DI.calibrateAffine=function(R1,R2,g1,g2,ang12){ const det=(a,b)=>a[0]*b[1]-a[1]*b[0]; const sgn=Math.sign(det(R1,R2))||1;
  const G1=[g1,0], G2=[g2*Math.cos(ang12*D2R),sgn*g2*Math.sin(ang12*D2R)]; // 관측과 같은 손잡이
  const D=det(R1,R2); const Ri=[[R2[1]/D,-R2[0]/D],[-R1[1]/D,R1[0]/D]]; // [R1 R2]⁻¹
  const T=[[G1[0]*Ri[0][0]+G2[0]*Ri[1][0],G1[0]*Ri[0][1]+G2[0]*Ri[1][1]],[G1[1]*Ri[0][0]+G2[1]*Ri[1][0],G1[1]*Ri[0][1]+G2[1]*Ri[1][1]]];
  // S = sqrt(TᵀT)
  const a=T[0][0]*T[0][0]+T[1][0]*T[1][0], b=T[0][0]*T[0][1]+T[1][0]*T[1][1], c=T[0][1]*T[0][1]+T[1][1]*T[1][1];
  const tr=a+c, dd=Math.sqrt(Math.max(0,(a-c)*(a-c)/4+b*b)); const l1=tr/2+dd, l2=tr/2-dd; const smax=Math.sqrt(l1), smin=Math.sqrt(Math.max(l2,0));
  const phiMax=Math.atan2(l1-a,b||1e-30); // 큰 고윳값 방향 (화면 px 축 기준)
  return {T,scale:Math.sqrt(smax*smin),smax,smin,ellip:smax/smin,phi:((phiMax%Math.PI)+Math.PI)%Math.PI}; };
/* 아핀 보정 적용: 관측 px 벡터 → 1/nm */
DI.applyAffine=function(cal,v){ if(cal.T) return [cal.T[0][0]*v[0]+cal.T[0][1]*v[1],cal.T[1][0]*v[0]+cal.T[1][1]*v[1]]; return [v[0]*cal.scale,v[1]*cal.scale]; };
/* 보정 결과(축척 + 찌그러짐 φ, ellip) 를 대칭 행렬로 — 회전 없이 저장·재사용 (smax 방향 φ) */
DI.affineFromParams=function(scale,ellip,phi){ const smax=scale*Math.sqrt(ellip), smin=scale/Math.sqrt(ellip); const c=Math.cos(phi), s=Math.sin(phi);
  return {scale,ellip,phi,smax,smin,T:[[smax*c*c+smin*s*s,(smax-smin)*c*s],[(smax-smin)*c*s,smax*s*s+smin*c*c]]}; };
/* 보정한 R1·R2 길이·각 (1/nm, 도) */
DI.correctLattice=function(fit,cal){ const a=DI.applyAffine(cal,fit.R1), b=DI.applyAffine(cal,fit.R2); const c=[a[0]+b[0],a[1]+b[1]]; const L=v=>Math.hypot(v[0],v[1]);
  const ang=Math.acos(Math.max(-1,Math.min(1,(a[0]*b[0]+a[1]*b[1])/(L(a)*L(b)))))/D2R; return {g:[L(a),L(b),L(c)],ang12:ang}; };

if(typeof module!=="undefined"&&module.exports) module.exports=DI;
G.DIFF_INDEX=DI;
})();
