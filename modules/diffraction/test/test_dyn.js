// 동역학 회절 (diff_dyn.js) 검증 — 해석해(두 빔·얇은 결정), 문헌 소멸 거리, 독립 Bloch(gemmi+scipy expm), 다중 슬라이스(abTEM 퍼텐셜),
// Gjønnes–Moodie·이중 회절, 수렴(Bethe), 속도.  기준 자료: python3 -I ref/gen_dyn_ref.py → test/ref_dyn.json
const fs=require("fs"), path=require("path");
require("../src/diff_sg.js"); require("../src/diff_scat.js"); const DC=require("../src/diff_core.js"); const PACK=require("../src/diff_pack.js"); const DD=require("../src/diff_dyn.js");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } if(!c||process.env.V) console.log((c?"  ✓ ":"  ✗ ")+m+(d?" — "+d:"")); };
const pk=id=>DC.phaseFromRecord(PACK.phases.find(p=>p.id===id));
const REF=JSON.parse(fs.readFileSync(path.join(__dirname,"ref_dyn.json"),"utf8"));
const ID={Si:"pack:Si",Au:"pack:Au",aFe:"pack:aFe",ZnO:"pack:ZnO",rutile:"pack:TiO2-rutile",bNaYF4:"pack:bNaYF4",Y2O3:"pack:Y2O3",LiYF4:"pack:LiYF4"};
const key=h=>h.join(",");
const mapOf=r=>{ const M=new Map(r.beams.map(b=>[key(b.h),b.I])); M.set("0,0,0",r.direct.I); return M; };
const fx=(x,n)=>Number(x).toFixed(n==null?3:n), pc=x=>(x*100).toFixed(2)+" %";
const SUM={}; // 문서용 숫자

console.log("== 상수·고유값 풀이");
{ const l=DD.wavelength(200), g=DD.gammaRel(200); ok(Math.abs(l-0.0250793)<2e-7&&Math.abs(g-1.391391)<1e-6,`200 kV: λ = ${l.toFixed(7)} Å, γ = ${g.toFixed(6)} (0.0250793 Å, 1.391391)`);
  const n=80; const A=new Float64Array(n*n); let s=7; const rnd=()=>{ s=(s*16807)%2147483647; return s/2147483647-0.5; };
  for(let i=0;i<n;i++) for(let j=0;j<=i;j++){ const v=rnd()+(i===j?0.01*i:0); A[i*n+j]=v; A[j*n+i]=v; }
  // 축퇴: 같은 값 두 줄
  const E=DD._symEig(A,n); let res=0, orth=0; for(let k=0;k<n;k++){ const c=E.C[k]; for(let i=0;i<n;i++){ let t=0; for(let j=0;j<n;j++) t+=A[i*n+j]*c[j]; res=Math.max(res,Math.abs(t-E.d[k]*c[i])); } for(let m=0;m<n;m++){ let t=0; const d=E.C[m]; for(let i=0;i<n;i++) t+=c[i]*d[i]; orth=Math.max(orth,Math.abs(t-(k===m?1:0))); } }
  ok(res<1e-12&&orth<1e-12,`tred2+tql2 80×80: 잔차 ${res.toExponential(1)}, 직교 ${orth.toExponential(1)}`);
  const D=new Float64Array(9); D[0]=2; D[4]=2; D[8]=1; D[1]=D[3]=0; const E2=DD._symEig(D,3); ok(Math.abs(E2.d[0]-1)<1e-14&&Math.abs(E2.d[2]-2)<1e-14,"축퇴 고유값 (2,2,1)"); }

console.log("== 두 빔 Pendellösung (해석해)");
{ const ph=pk("pack:Si"); const h=[2,2,0]; const r0=DD.zoneIntensities(ph,[0,0,1],{thickness:[1],beams:[h],gm:false}); const b0=DD.beamOf(r0,h); const K=r0.meta.K;
  const T=[]; for(let t=0;t<=300;t+=7.5) T.push(t);
  // 정확한 Bragg: k_t = −g/2 → s = 0
  const tilt=[-1000*b0.x/2/K,-1000*b0.y/2/K]; const r=DD.zoneIntensities(ph,[0,0,1],{thickness:T,beams:[h],tilt,gm:false}); const b=DD.beamOf(r,h); const xi=b.xi;
  let e=0; T.forEach((t,i)=>{ e=Math.max(e,Math.abs(b.I[i]-Math.sin(Math.PI*t/xi)**2)); }); ok(Math.abs(b.s)<1e-12&&e<1e-12,`Si 220 s = 0: I_g = sin²(πt/ξ), ξ = ${fx(xi,2)} nm, 최대 차 ${e.toExponential(1)}`);
  // 정대축 (s = −g²/2K): I = sin²(πt√(s²+1/ξ²))/(1+(sξ)²)
  const r2=DD.zoneIntensities(ph,[0,0,1],{thickness:T,beams:[h],gm:false}); const c=DD.beamOf(r2,h); const s=c.s; const w=s*xi; let e2=0; T.forEach((t,i)=>{ e2=Math.max(e2,Math.abs(c.I[i]-Math.sin(Math.PI*t*Math.hypot(s,1/xi))**2/(1+w*w))); });
  ok(e2<1e-12,`Si 220 s = ${s.toExponential(3)} nm⁻¹ (w = ${fx(w,2)}): 들뜸 오차 Pendellösung 최대 차 ${e2.toExponential(1)}`);
  const sum=T.map((_,i)=>c.I[i]+r2.direct.I[i]); ok(sum.every(x=>Math.abs(x-1)<1e-12),"두 빔 Σ I = 1");
  // 흡수: 두 빔 + 비례 흡수 κ (행렬 지수 걸음) ↔ 해석해
  const kap=0.08; const r3=DD.zoneIntensities(ph,[0,0,1],{thickness:T,beams:[h],tilt,absorption:kap,gm:false}); const d=DD.beamOf(r3,h); const q0=kap*r3.meta.U0/(2*K);
  // 고유값 i q0 ± (1+iκ)/(2ξ) → I_g = ½ e^{−4π q0 t} [cosh(2πκt/ξ) − cos(2πt/ξ)]
  let e3=0; T.forEach((t,i)=>{ const an=0.5*Math.exp(-4*Math.PI*q0*t)*(Math.cosh(2*Math.PI*kap*t/xi)-Math.cos(2*Math.PI*t/xi)); e3=Math.max(e3,Math.abs(d.I[i]-an)); });
  ok(e3<1e-9,`두 빔 + 흡수 κ = ${kap} (행렬 지수 걸음): 해석해와 최대 차 ${e3.toExponential(1)}`);
  SUM.xiSi220_200=xi; }

console.log("== 소멸 거리 ξ = K/|U_g| ↔ 문헌 (Hirsch et al. 1965 표, 100 kV)");
{ const H=[["pack:Al",[1,1,1],[1,-1,0],56.3],["pack:Al",[2,0,0],[0,0,1],68.5],["pack:Si",[1,1,1],[1,-1,0],60.2],["pack:Si",[2,2,0],[0,0,1],75.7],["pack:Si",[4,0,0],[0,0,1],126.8],
    ["pack:Au",[1,1,1],[1,-1,0],18.3],["pack:Au",[2,0,0],[0,0,1],20.2],["pack:Au",[2,2,0],[0,0,1],27.8],["pack:aFe",[1,1,0],[0,0,1],28.6],["pack:aFe",[2,0,0],[0,0,1],41.2],
    ["pack:Cu",[1,1,1],[1,-1,0],28.6],["pack:Cu",[2,0,0],[0,0,1],32.6],["pack:MgO",[2,0,0],[0,0,1],46.1],["pack:MgO",[2,2,0],[0,0,1],66.2]];
  let worst=0; const rows=[];
  for(const [id,h,z,lit] of H){ const r=DD.zoneIntensities(pk(id),z,{kV:100,thickness:[1],beams:[h],Bset:0,gm:false}); const x=DD.beamOf(r,h).xi; const e=x/lit-1; worst=Math.max(worst,Math.abs(e)); rows.push(`${id.slice(5)} ${h.join("")} ${fx(x,1)}/${lit}`); }
  ok(worst<0.06,`14 반사 (B = 0) 최대 차 ${pc(worst)}: ${rows.join(", ")}`);
  const al=DD.beamOf(DD.zoneIntensities(pk("pack:Al"),[1,-1,0],{kV:100,thickness:[1],beams:[[1,1,1]],Bset:0,gm:false}),[1,1,1]).xi; ok(Math.abs(al-56)<1.5,`Al 111, 100 kV: ξ = ${fx(al,1)} nm (Hirsch ≈ 56 nm)`);
  const si4=DD.beamOf(DD.zoneIntensities(pk("pack:Si"),[0,0,1],{kV:200,thickness:[1],beams:[[4,0,0]],Bset:0.46,gm:false}),[4,0,0]).xi; ok(si4>164&&si4<224,`Si 400, 200 kV (B 0.46): ξ = ${fx(si4,1)} nm (측정 194 ± 30 nm, Acta Phys. Sin. 71 (2022), doi:10.7498/aps.71.20212271)`);
  const s220=DD.beamOf(DD.zoneIntensities(pk("pack:Si"),[0,0,1],{kV:200,thickness:[1],beams:[[2,2,0]],Bset:0.46,gm:false}),[2,2,0]).xi; SUM.xi220=s220; SUM.al=al; SUM.si4=si4; SUM.xiWorst=worst; }

// Bethe 는 t ≫ 1/|s_w| (≈ 1 nm) 에서만 맞는 단열 소거라 t → 0 극한을 바꾼다 → 이 시험은 Bethe 없이 (모든 빔 강하게)
console.log("== 얇은 결정 극한 (t = 0.02 nm, Bethe 없이): I_g → (πt/ξ_g)² sinc²(π s_g t) ∝ |F_g|²");
for(const [id,z] of [["pack:Si",[0,1,1]],["pack:Au",[0,0,1]],["pack:Y2O3",[0,0,1]],["pack:ZnO",[1,1,0]]]){ const t=0.02; const r=DD.zoneIntensities(pk(id),z,{thickness:[t],gm:false,gBethe:0,gStrong:12,maxBeams:2000}); let e=0, n=0, eF=0;
  const L=r.beams.filter(b=>b.IkinRel>1e-2); const ref=L.slice().sort((a,b)=>a.g-b.g)[0];
  for(const b of L){ const pred=(Math.PI*t/b.xi)**2*(b.s?Math.sin(Math.PI*b.s*t)**2/(Math.PI*b.s*t)**2:1); e=Math.max(e,Math.abs(b.I[0]/pred-1)); n++;
    if(b.g<6){ const ratio=(b.I[0]/DD.beamOf(r,ref.h).I[0])/(b.Ikin/ref.Ikin); eF=Math.max(eF,Math.abs(ratio-1)); } }
  ok(e<0.02&&eF<0.02,`${id.slice(5)} [${z}]: ${n} 빔 (|F|² > 1 % 최대) 운동학 1 차식과 최대 ${pc(e)}, g < 6 nm⁻¹ 상대 세기 ↔ |F|² 최대 ${pc(eF)}`); }

console.log("== 독립 Bloch (gemmi c4322 + scipy expm, 같은 빔 목록)");
{ let worst=0; for(const c of REF.bloch){ const ph=pk(ID[c.phase]); const o={kV:c.kV,thickness:c.thickness,absorption:c.kappa,beams:c.beams,gmax:c.gB+1,gm:false,dw:c.dw,Bdefault:0};
    if(c.kt.some(x=>x)){ const r0=DD.zoneIntensities(ph,c.zone,Object.assign({},o,{thickness:[1]})); o.tilt=[1000*DC.vec.dot(c.kt,r0.e1)/r0.meta.K,1000*DC.vec.dot(c.kt,r0.e2)/r0.meta.K]; }
    const r=DD.zoneIntensities(ph,c.zone,o); const M=mapOf(r); let e=0;
    for(let k=0;k<c.thickness.length;k++){ let mx=0; c.I.forEach(I=>{ mx=Math.max(mx,I[k]); }); c.beams.forEach((h,i)=>{ e=Math.max(e,Math.abs(M.get(key(h))[k]-c.I[i][k])/mx); }); }
    const auto=DD.zoneIntensities(ph,c.zone,{kV:c.kV,thickness:[1],gmax:c.gB,gStrong:c.gB,maxBeams:9999,gBethe:0,gm:false}).meta.nBeams;
    // ZnO·β-NaYF4: CIF 좌표 0.3333/0.6667 을 gemmi 는 특수 자리 1/3·2/3 로 다르게 겹친다 → 1e-3 정도 차
    const tol=/ZnO|bNaYF4/.test(c.phase)?5e-3:1e-5; worst=Math.max(worst,tol===1e-5?e:0);
    ok(e<tol&&auto===c.N,`${c.phase} [${c.zone}] ${c.kV} kV κ=${c.kappa}${c.kt.some(x=>x)?" 기울임":""}: ${c.N} 빔 (자동 목록 ${auto}), ${r.meta.solver}, 최대 차/최대 ${e.toExponential(1)}`); }
  SUM.blochWorst=worst; }

console.log("== 원점·경로 불변 (반전 중심 이동 실수 경로 ↔ CIF 원점 복소 2N 경로)");
{ const ph=pk("pack:Si"); const a=DD.zoneIntensities(ph,[0,1,1],{thickness:[20,50],gm:false}), b=DD.zoneIntensities(ph,[0,1,1],{thickness:[20,50],gm:false,origin:"cif"}); let e=0; a.beams.forEach((x,i)=>{ x.I.forEach((v,k)=>{ e=Math.max(e,Math.abs(v-b.beams[i].I[k])); }); });
  ok(a.meta.solver==="real"&&b.meta.solver==="complex"&&e<1e-9,`Si [011]: ${a.meta.solver} N=${a.meta.dim} ↔ ${b.meta.solver} 2N=${b.meta.dim}, 최대 차 ${e.toExponential(1)}`);
  const c=DD.zoneIntensities(ph,[0,1,1],{thickness:[10,40,90],weak:false,gm:false,gmax:36}); /* 출력에 강한 빔 전부 */ const s=c.thickness.map((_,k)=>c.direct.I[k]+c.beams.reduce((t,x)=>t+x.I[k],0)); ok(s.every(x=>Math.abs(x-1)<1e-10),`흡수 없음 Σ I = 1 (Bethe 포함 Hermite): ${s.map(x=>x.toFixed(12)).join(", ")}`);
  const d=DD.zoneIntensities(ph,[0,1,1],{thickness:[10,40,90],weak:false,gm:false,gmax:36,absorption:0.1}); const s2=d.thickness.map((_,k)=>d.direct.I[k]+d.beams.reduce((t,x)=>t+x.I[k],0)); ok(s2[0]<1&&s2[1]<s2[0]&&s2[2]<s2[1],`흡수 κ = 0.1: Σ I 감소 ${s2.map(x=>fx(x,4)).join(" → ")}`);
  const zp=DC.zonePattern(ph,[0,1,1],12); let ex=0, nn=0; for(const sp of zp.spots){ const bb=DD.beamOf(a,sp.h); if(!bb) continue; nn++; ex=Math.max(ex,Math.hypot(bb.x-sp.x,bb.y-sp.y)); } ok(nn===zp.spots.length&&ex<1e-9,`2D 좌표 = DC.zonePattern (${nn} 스폿, 최대 차 ${ex.toExponential(1)})`); }

console.log("== 다중 슬라이스 대조 (abTEM 퍼텐셜, B = 0.5 Å²)");
{ const rows=[]; SUM.ms=[];
  for(const c of REF.ms){ const ph=pk(ID[c.phase]); if(c.tilt.some(x=>x)) continue;
    const o={thickness:c.thickness,Bset:c.B,gm:false};
    const conv=DD.zoneIntensities(ph,c.zone,Object.assign({},o,{gStrong:45,maxBeams:900,gBethe:60,maxBethe:6000})); const def=DD.zoneIntensities(ph,c.zone,o);
    const st=(r)=>{ const M=mapOf(r); let eS=0, eM=0, e0=0; for(let k=0;k<c.thickness.length;k++){ let mx=0; c.beams.forEach((h,i)=>{ if(h.some(x=>x)) mx=Math.max(mx,c.I[i][k]); });
        c.beams.forEach((h,i)=>{ const v=M.get(key(h)); if(!v) return; const d=Math.abs(v[k]-c.I[i][k]); if(!h.some(x=>x)) e0=Math.max(e0,d); else { eM=Math.max(eM,d/mx); if(c.I[i][k]>0.2*mx) eS=Math.max(eS,d/c.I[i][k]); } }); } return {eS,eM,e0}; };
    const a=st(conv), b=st(def); SUM.ms.push({c,a,b,Nc:conv.meta.nBeams,Nd:def.meta.nBeams});
    // zolz: 같은 0 차 Laue 모형의 독립 풀이 → 슬라이스·표본화 오차만 (Au 는 Δx 0.04 Å·Δz c/8 에서 수 %, Δx 0.025 Å 이면 반으로)
    // abtem: 3D 슬라이스 (HOLZ 포함) — Si 만 판정, 나머지는 HOLZ 영향 기록
    const lim=c.mode==="zolz"?(c.phase==="Au"?{eM:c.sampling<0.04?0.03:c.nsub>8?0.05:0.07,e0:0.02}:{eM:0.03,e0:0.01}):(c.phase==="Si"?{eM:0.02,e0:0.01}:{eM:Infinity,e0:Infinity});
    ok(a.eM<lim.eM&&a.e0<lim.e0,`${c.phase} [${c.zone}] ${c.mode}${c.nsub!==8?" nsub "+c.nsub:""}${c.sampling!==0.04?" Δx "+c.sampling:""}: 수렴 Bloch(N ${conv.meta.nBeams}) 회절빔 최대 차/최대 ${pc(a.eM)} · 강한 빔 상대 ${pc(a.eS)} · I0 ${a.e0.toExponential(1)} | 기본값(N ${def.meta.nBeams}) ${pc(b.eM)} · ${pc(b.eS)}`); } }

console.log("== 순위 상관: 3D 다중 슬라이스(HOLZ 포함, 3 두께 평균) 세기를 '실측' 으로 보고 ρ(동역학 ZOLZ) ↔ ρ(운동학 |F|²)");
{ const DI=require("../src/diff_index.js"); SUM.rho=[];
  for(const c of REF.ms){ if(c.mode!=="abtem") continue; const ph=pk(ID[c.phase]); const r=DD.zoneIntensities(ph,c.zone,{thickness:c.thickness,Bset:c.B,gm:false}); const M=new Map(r.beams.map(b=>[key(b.h),b]));
    const idx=c.beams.map((h,i)=>i).filter(i=>c.beams[i].some(x=>x)&&M.get(key(c.beams[i]))); const nt=c.thickness.length;
    const obs=idx.map(i=>c.I[i].reduce((a,b)=>a+b,0)/nt), dyn=idx.map(i=>M.get(key(c.beams[i])).I.reduce((a,b)=>a+b,0)/nt), kin=idx.map(i=>M.get(key(c.beams[i])).Ikin);
    const top=idx.map((_,j)=>j).sort((a,b)=>obs[b]-obs[a]).slice(0,20); const rd=DI.spearman(obs,dyn), rk=DI.spearman(obs,kin), td=DI.spearman(top.map(j=>obs[j]),top.map(j=>dyn[j])), tk=DI.spearman(top.map(j=>obs[j]),top.map(j=>kin[j]));
    SUM.rho.push([c.phase,c.zone,idx.length,rd,rk,td,tk]);
    ok(rd>=0.85&&rd>=rk-0.01,`${c.phase} [${c.zone}] ${idx.length} 스폿: ρ 동역학 ${fx(rd)} · 운동학 ${fx(rk)} | 센 20 개 ${fx(td)} · ${fx(tk)}`); } }

console.log("== 이중 회절 · Gjønnes–Moodie");
{ const T=[10,20,30,60]; const zr=DC.parseCIF(fs.readFileSync(path.join(__dirname,"cif","ZrO2_P21c.cif"),"utf8"))[0];
  const R=(ph,z,o)=>DD.zoneIntensities(ph,z,Object.assign({thickness:T},o||{})); const rel=(r,b,k)=>b.I[k]/r.beams.reduce((m,x)=>Math.max(m,x.I[k]),0);
  { const r=R(pk("pack:Si"),[0,1,1]); const b=DD.beamOf(r,[2,0,0]); ok(b.absent&&b.kinForbidden&&!b.gmExtinct&&b.nRoute>0&&b.I.slice(1).every(x=>x>1e-3),`Si [011] (200): 공간군 소멸·운동학 0 인데 동역학 세기 ${b.I.map(x=>x.toExponential(1)).join("/")} (t ${T}), 경로 ${b.nRoute} (예: ${b.ddRoute[0].map(key).join(" + ")}), GM 아님`);
    const c=DD.beamOf(r,[2,2,-2]); ok(c.kinForbidden&&!c.absent&&c.I[1]>1e-5,`Si [011] (22-2): 구조상 0 (공간군 허용) → 이중 회절 ${c.I[1].toExponential(1)}`); SUM.si200=b.I; }
  { const r=R(pk("pack:Si"),[0,0,1]); const b=DD.beamOf(r,[2,0,0]); ok(b.dynZero&&b.nRoute===0&&!b.gmExtinct&&Math.max(...b.I)<1e-20,`Si [001] (200): 0 차 Laue 영역 안 경로 없음 (허용 반사 h+k ≡ 0 mod 4) → ${Math.max(...b.I).toExponential(1)} (실측에서 약하게 보이면 HOLZ 경유)`); }
  { const r=R(pk("pack:TiO2-rutile"),[0,0,1]); const L=[[1,0,0],[0,1,0],[3,0,0],[0,3,0]].map(h=>DD.beamOf(r,h)); ok(L.every(b=>b.gmExtinct&&b.nRoute>0&&Math.max(...b.I)<1e-20),`루틸 P4₂/mnm [001] (100)(010)(300)(030): 경로 ${L.map(b=>b.nRoute).join("/")} 있지만 GM 소멸 (최대 ${Math.max(...L.map(b=>Math.max(...b.I))).toExponential(1)})`);
    const rt=R(pk("pack:TiO2-rutile"),[0,0,1],{tilt:[5,1.5]}); const b=DD.beamOf(rt,[1,0,0]); ok(b.gmExtinct&&b.I.slice(1).every(x=>x>1e-4),`루틸 5.2 mrad 기울임: (100) = ${b.I.map(x=>x.toExponential(1)).join("/")} (정대축 GM 표시는 유지)`);
    const ms=REF.ms.find(c=>c.phase==="rutile"&&!c.tilt.some(x=>x)&&c.mode==="zolz"), mt=REF.ms.find(c=>c.phase==="rutile"&&c.tilt.some(x=>x));
    const i1=ms.beams.findIndex(h=>key(h)==="1,0,0"), i2=mt.beams.findIndex(h=>key(h)==="1,0,0"); const mx=Math.max(...mt.I.map((I,i)=>mt.beams[i].some(x=>x)?I[0]:0));
    ok(Math.max(...ms.I[i1])<1e-4*mt.I[i2][0]&&mt.I[i2][0]>1e-3*mx,`다중 슬라이스도 같음: 정대축 (100) ${Math.max(...ms.I[i1]).toExponential(1)} (격자 표본화로 대칭이 조금 깨진 수치 잡음), 기울임 ${mt.I[i2][0].toExponential(1)} (회절 최대 대비 ${fx(mt.I[i2][0]/mx,3)})`); }
  { const r=R(pk("pack:ZnO"),[1,1,0]); const L=[[0,0,1],[0,0,3]].map(h=>DD.beamOf(r,h)); ok(L.every(b=>b.gmExtinct),`ZnO P6₃mc [110] (001)(003): 6₃ 나사축·c 미끄럼 → GM 소멸`);
    const p=DD.beamOf(r,[0,0,2]), m=DD.beamOf(r,[0,0,-2]); ok(Math.abs(p.I[2]/m.I[2]-1)>0.05,`ZnO 극성: Friedel 깨짐 I(002)/I(00-2) = ${p.I.map((x,k)=>fx(x/m.I[k],2)).join("/")} (운동학은 1)`); }
  { const r=R(zr,[1,0,0]); const L=[[0,0,1],[0,1,0],[0,0,3],[0,3,0]].map(h=>DD.beamOf(r,h)); ok(L.every(b=>b&&b.gmExtinct),`ZrO2 P2₁/c [100] (001)(010)(003)(030): c 미끄럼·2₁ → GM 소멸`);
    const rt=R(zr,[1,0,0],{tilt:[3,2]}); const b=DD.beamOf(rt,[0,0,1]), c=DD.beamOf(rt,[0,1,0]); ok(b.I[2]>1e-5&&c.I[2]>1e-5,`ZrO2 3.6 mrad 기울임: (001) ${b.I[2].toExponential(1)}, (010) ${c.I[2].toExponential(1)} (t = 30 nm)`); }
  { const r=R(pk("pack:Fe3O4"),[1,1,0]); const b=DD.beamOf(r,[0,0,2]); ok(b.absent&&!b.gmExtinct&&b.I[1]>1e-3,`Fe3O4 [110] (002): 소멸이지만 이중 회절 ${b.I[1].toExponential(1)} (GM 아님)`); } }

console.log("== 수렴: 기본값 (강한 ≤ 200 빔 + Bethe ≤ 50 nm⁻¹) ↔ 큰 빔 집합, 두께 10–80 nm 평균");
{ SUM.conv=[]; const T=[]; for(let t=10;t<=80;t+=5) T.push(t);
  for(const [id,z,lim] of [["pack:Si",[0,1,1],0.01],["pack:Au",[0,0,1],0.02],["pack:Au",[0,1,1],0.03],["pack:aFe",[1,1,1],0.02],["pack:LiYF4",[1,0,0],0.02]]){
    const ph=pk(id); const ref=DD.zoneIntensities(ph,z,{thickness:T,gStrong:45,maxBeams:1000,gBethe:60,maxBethe:6000,gm:false}); const def=DD.zoneIntensities(ph,z,{thickness:T,gm:false});
    const A=DD.thicknessAverage(ref).I, B=DD.thicknessAverage(def).I; const M=new Map(ref.beams.map((b,i)=>[key(b.h),A[i]])); let mx=Math.max(...A), e=0; def.beams.forEach((b,i)=>{ e=Math.max(e,Math.abs(B[i]-M.get(key(b.h)))/mx); });
    SUM.conv.push([id,z,def.meta.nBeams,def.meta.nBethe,ref.meta.nBeams,e,def.meta.ms]);
    ok(e<lim,`${id.slice(5)} [${z}]: 기본 N ${def.meta.nBeams} + Bethe ${def.meta.nBethe} ↔ N ${ref.meta.nBeams}: 평균 세기 최대 차/최대 ${pc(e)} (${fx(def.meta.ms,0)} ms)`); }
  // 바깥 빔 1 차 섭동: Si [011] 강한 ≤ 12, 출력 ≤ 16
  const ph=pk("pack:Si"); const r=DD.zoneIntensities(ph,[0,1,1],{thickness:T,gmax:16,gStrong:12,gm:false}); const ref=DD.zoneIntensities(ph,[0,1,1],{thickness:T,gmax:16,gStrong:45,maxBeams:1000,gBethe:60,gm:false});
  const A=DD.thicknessAverage(ref).I, B=DD.thicknessAverage(r).I; const M=new Map(ref.beams.map((b,i)=>[key(b.h),A[i]])); const mx=Math.max(...A); let e=0, n=0; r.beams.forEach((b,i)=>{ if(!b.pert) return; n++; e=Math.max(e,Math.abs(B[i]-M.get(key(b.h)))/mx); });
  ok(n>0&&e<0.01,`바깥 빔 섭동 (Si [011] 12 < g ≤ 16 nm⁻¹, ${n} 빔): 평균 세기 최대 차/최대 ${pc(e)}`); }

console.log("== 두께 평균 · 미세 기울임 평균 · 캐시");
{ const ph=pk("pack:Au"); const a=DD.zoneAverage(ph,[0,0,1]); const mx=Math.max(...a.beams.map(b=>b.Irel)); ok(Math.abs(mx-1)<1e-12&&a.thickness.length===36&&a.meta.avg.tRange[1]===80,`zoneAverage: 10–80 nm 2 nm 간격 ${a.thickness.length} 점, Irel 최대 1`);
  ok(DD.zoneAverage(ph,[0,0,1])===a,"같은 조건 다시 부르면 캐시");
  const t=DD.zoneAverage(ph,[0,0,1],{tiltAvg:1,nTilt:6}); ok(t.meta.avg.nTilt===7&&t!==a,`기울임 평균 (1 mrad 원 6 + 중심): ${fx(t.meta.avg.ms,0)} ms`);
  const r=DD.zoneIntensities(ph,[0,0,1],{thickness:[10,20,30],gm:false}); const av=DD.thicknessAverage(r,[10,30]); const b=r.beams[0]; ok(Math.abs(av.I[0]-(b.I[0]*5+b.I[1]*10+b.I[2]*5)/20)<1e-15,"사다리꼴 두께 평균"); }

console.log("== 속도");
{ const ph=pk("pack:Si"); DD.zoneIntensities(ph,[0,1,1],{maxBeams:120,gBethe:0,gm:false}); const t0=Date.now(); const r=DD.zoneIntensities(ph,[0,1,1],{maxBeams:120,gBethe:0,gm:false}); const ms=Date.now()-t0;
  ok(ms<300,`Si [011] 120 빔 (실수, 10 두께): ${ms} ms (풀이 ${fx(r.meta.msSolve,0)} ms)`); SUM.t120=ms;
  const rc=DD.zoneIntensities(ph,[0,1,1],{maxBeams:120,gBethe:0,gm:false,origin:"cif"}); ok(rc.meta.ms<600,`같은 계산 복소 2N 경로: ${fx(rc.meta.ms,0)} ms`); SUM.t120c=rc.meta.ms;
  let worst=0, wn=""; const L=[]; for(const rec of PACK.phases){ const p=DC.phaseFromRecord(rec); const zt=DC.zoneTable(p,{maxIdx:2}).slice(0,2);
    for(const z of zt){ const t1=Date.now(); const q=DD.zoneIntensities(p,z.zone,{}); const dt=Date.now()-t1; L.push(dt); if(dt>worst){ worst=dt; wn=`${rec.id.slice(5)} [${z.zone}] N ${q.meta.nBeams} ${q.meta.solver} Bethe ${q.meta.nBethe}`; } } }
  L.sort((a,b)=>a-b); ok(worst<1500,`내장 21 상 × 낮은 정대축 2 개 (기본값, GM 판정 포함): 중앙값 ${L[L.length>>1]} ms, 최악 ${worst} ms (${wn})`); SUM.speed={med:L[L.length>>1],worst,wn};
  const ta=Date.now(); const k=DD.zoneIntensities(pk("pack:Au"),[0,0,1],{absorption:0.1,gm:false}); ok(Date.now()-ta<3000,`흡수 κ = 0.1 (행렬 지수, N ${k.meta.nBeams}, 0.125 nm 걸음 800): ${Date.now()-ta} ms`); }

console.log("== 입력 경계");
{ let e=null; try{ DD.zoneIntensities(DC.phaseFromParams({sgNumber:225,a:0.4}),[0,0,1]); }catch(x){ e=x; } ok(e&&/원자/.test(e.message),"원자 없는 상 → 오류");
  const r=DD.zoneIntensities(pk("pack:Au"),[0,0,1],{thickness:[10],includeHOLZ:true,gm:false}); ok(r.meta.warn.length===1,"includeHOLZ → 경고 (미구현)");
  let e2=null; try{ DD.zoneIntensities(pk("pack:Au"),[0,0,1],{beams:[[1,1,1]]}); }catch(x){ e2=x; } ok(e2&&/Laue/.test(e2.message),"정대축 밖 빔 지정 → 오류"); }

if(process.env.SUMMARY) fs.writeFileSync(process.env.SUMMARY,JSON.stringify(SUM,(k,v)=>k==="c"?{phase:v.phase,zone:v.zone,mode:v.mode,nsub:v.nsub,sampling:v.sampling,thickness:v.thickness}:v,1));
console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
