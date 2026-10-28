/* 실측 대조: 사용자 실측 자료의 Avantage·CasaXPS 출력 ↔ INNO (SPEC_real_validation.md)
   node test/test_real_kara.js   (modules/xps 에서, 약 10 초)
   자료: test/data/real/kara/ (비공개 — 공개 소스 백업에서 빠짐). 폴더가 없으면 모든 항목을 건너뛰고 0 으로 끝난다.
   A) avantage_231107/{1..4}.xls (+ .vms): Avantage Peak Table + 영역별 'Backgnd.'(Smart) 열 — Kratos 계 장비(Hybrid lens, Slot), TXFN = 1
   B) casa_240213/: CasaXPS 2.3.25 Integrated.vms (CASA region/comp 주석) + 'Fitting' 내보내기(성분 곡선·배경·포락선)            */
const path=require("path"), fs=require("fs");
const C=require("../src/xps_core.js"), IO=require("../src/xps_io.js");
const ROOT=path.join(__dirname,"data","real","kara");
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
function info(name,v){ console.log("  info "+name+(v?"  ("+v+")":"")); }
const med=a=>{ const s=a.slice().sort((x,y)=>x-y); return s.length?s[s.length>>1]:NaN; };
const mx=a=>a.reduce((m,v)=>Math.max(m,v),-Infinity);
const e1=v=>v.toExponential(1), p2=v=>(100*v).toFixed(2)+" %";
const rd=f=>new Uint8Array(fs.readFileSync(f));

if(!fs.existsSync(ROOT)){ console.log("test/data/real/kara 없음 — 실측 대조를 건너뜀 (skip)"); process.exit(0); }

(async()=>{
/* ================= A) Avantage ================= */
const AV=path.join(ROOT,"avantage_231107");
if(fs.existsSync(AV)){
  console.log("== A0 읽기: Kratos→Avantage .xls (위쪽 'Axis | Kinetic energy' 줄) · 시트 이름 → Peak Table 이름 · .vms 와 같은 세기");
  const regs=[];   // {f, name, sp, i0, i1, T, row}
  for(const f of ["1","2","3","4"]){
    let X; try{ X=await IO.readFile(f+".xls",rd(path.join(AV,f+".xls"))); }catch(e){ ok(false,`${f}.xls 읽기`,e.message); continue; }
    const pt=X.tables.find(t=>/peak\s*table/i.test(t.name)); const H=pt.rows.find(r=>r.some(v=>/^name$/i.test(String(v).trim()))).map(v=>String(v).trim()); const col=n=>H.indexOf(n);
    const rows=pt.rows.slice(pt.rows.findIndex(r=>r.some(v=>/^name$/i.test(String(v).trim())))+1).filter(r=>typeof r[col("Start BE")]==="number");
    let named=0; rows.forEach(r=>{ const sp=X.spectra.find(s=>s.name===String(r[col("Name")]).trim()); if(!sp) return; named++;
      const b=sp.extra.avantageCols[Object.keys(sp.extra.avantageCols).find(k=>/backgnd/i.test(k))]; const st=r[col("Start BE")], en=r[col("End BE")];
      let i0=0,i1=0; for(let i=0;i<sp.n;i++){ if(Math.abs(sp.be[i]-en)<Math.abs(sp.be[i0]-en)) i0=i; if(Math.abs(sp.be[i]-st)<Math.abs(sp.be[i1]-st)) i1=i; }
      let nz0=-1,nz1=-1; for(let i=0;i<sp.n;i++) if(b[i]){ if(nz0<0) nz0=i; nz1=i; }
      regs.push({f,name:sp.name,sp,i0,i1,T:b,row:Object.fromEntries(H.map((h,k)=>[h,r[k]])),bgRange:nz0===i0&&nz1===i1}); });
    ok(named===rows.length&&X.spectra.length===8,`${f}.xls: 스펙트럼 8 개, Peak Table ${rows.length} 줄 모두 이름 붙음`,`${named}/${rows.length}`);
    const V=await IO.readFile(f+".vms",rd(path.join(AV,f+".vms")));
    let worst=0, nm=0; X.spectra.forEach(x=>{ const v=V.spectra.find(v=>v.n===x.n&&Math.abs(v.be[0]-x.be[0])<0.02); if(!v) return; nm++; for(let i=0;i<x.n;i++) worst=Math.max(worst,Math.abs(x.y[i]-v.y[i])/Math.max(1,Math.abs(v.y[i]))); });
    ok(nm===8&&worst<1e-5,`${f}: .xls 세기 = .vms 세기/(dwell·scans) — Avantage 는 Kratos 투과 함수 열을 쓰지 않음`,`최대 상대차 ${e1(worst)}, .vms 투과 열 ${V.spectra.every(s=>s.tfArr)?"있음":"없음"}`);
  }
  ok(regs.length===26&&regs.every(r=>r.bgRange),"Avantage 배경 열의 0 아닌 구간 = Start/End BE 에 가장 가까운 채널 (26/26)",`${regs.filter(r=>r.bgRange).length}/${regs.length}`);
  const net=(r,B)=>{ const E=r.sp.be, y=r.sp.y; let A=0,H=0; for(let i=r.i0;i<=r.i1;i++){ if(i>r.i0) A+=0.5*((y[i]-B[i-r.i0])+(y[i-1]-B[i-1-r.i0]))*(E[i]-E[i-1]); H=Math.max(H,y[i]-B[i-r.i0]); } return {A,H}; };
  const Tof=r=>r.T.slice(r.i0,r.i1+1);

  console.log("== A1 끝점 = 9 점 중심 평균 (창 밖 4 채널 포함)");
  let ne=0; regs.forEach(r=>{ const y=r.sp.y; const av=c=>{ let s=0; for(let j=c-4;j<=c+4;j++) s+=y[j]; return s/9; };
    const okk=[[r.i0,r.T[r.i0]],[r.i1,r.T[r.i1]]].every(([c,t])=>Math.abs(av(c)-t)<=Math.max(0.06,6e-6*Math.abs(t))); if(okk) ne++; });
  ok(ne===26,"Avantage Smart 끝값 = 끝 채널 ±4 의 평균 (6 자리 인쇄 오차 안)",`${ne}/26`);

  console.log("== A2 Smart 배경: INNO (Avantage 규약) vs Avantage 'Backgnd.' — 최대|ΔB|/봉우리 높이, 넓이(P) 상대차");
  const rel={av:[],inno:[]}, ar={av:[],inno:[]}, arAbs={av:[],inno:[]};
  regs.forEach(r=>{ const E=r.sp.be.slice(r.i0,r.i1+1), y=r.sp.y.slice(r.i0,r.i1+1), T=Tof(r), ref=net(r,T);
    const Bav=C.bgStatic(E,y,{type:"smart",conv:"avantage",_pad:C.avantagePad(r.sp.y,r.i0,r.i1)}), Bin=C.bgStatic(E,y,{type:"smart",n:5});
    [["av",Bav],["inno",Bin]].forEach(([k,B])=>{ let m=0; for(let i=0;i<E.length;i++) m=Math.max(m,Math.abs(B[i]-T[i])); rel[k].push(m/ref.H); const a=net(r,B).A/ref.A-1; ar[k].push(a); arAbs[k].push(Math.abs(a)); }); });
  info("INNO 기본 Smart (n 5, 평활 데이터로 점별 자름)",`배경 중앙 ${e1(med(rel.inno))} · 최대 ${e1(mx(rel.inno))}; 넓이 |Δ| 중앙 ${p2(med(arAbs.inno))} · 최대 ${p2(mx(arAbs.inno))}`);
  ok(med(rel.av)<6e-3&&mx(rel.av)<0.06,"Avantage 규약 Smart 배경 ≈ Avantage",`중앙 ${e1(med(rel.av))} · 최대 ${e1(mx(rel.av))}`);
  ok(med(arAbs.av)<5e-3&&mx(arAbs.av)<0.05,"Avantage 규약 Smart 넓이(P) ≈ Avantage",`|Δ| 중앙 ${p2(med(arAbs.av))} · 최대 ${p2(mx(arAbs.av))}`);
  ok(med(arAbs.av)<med(arAbs.inno)/3&&mx(arAbs.av)<mx(arAbs.inno)/3,"Avantage 규약이 기본 Smart 보다 Avantage 넓이에 가깝다",`최대 ${p2(mx(arAbs.av))} vs ${p2(mx(arAbs.inno))}`);
  const fe=regs.map((r,i)=>[r,ar.av[i]]).filter(([r])=>/^Fe/.test(r.name)).map(([,a])=>a);
  info("남은 차이가 큰 영역: Fe 2p (넓은 창 23.8 eV, 큰 계단)",fe.map(p2).join(", ")+" — Avantage Shirley 고정점 식의 차이(미확인, SPEC_real_validation §A2)");
  const st=C.makeProblem({be:regs[0].sp.be,y:regs[0].sp.y,c:null},{lo:regs[0].sp.be[regs[0].i0],hi:regs[0].sp.be[regs[0].i1],bg:{type:"smart",conv:"avantage"},comps:[],cons:{}});
  ok(st.bg._pad&&st.bg._pad.lo.length===4&&st.bg._pad.hi.length===4,"makeProblem 이 Avantage Smart 에 창 밖 4 채널을 넘김");

  console.log("== A3 넓이(P) 적분 규약");
  const dA=regs.map(r=>Math.abs(net(r,Tof(r)).A/r.row["Area (P) CPS.eV"]-1));
  ok(mx(dA)<1e-3,"Area (P) = Avantage 배경 위 사다리꼴 적분 (끝 채널 = 가장 가까운 채널)",`중앙 ${e1(med(dA))} · 최대 ${e1(mx(dA))}`);

  console.log("== A4 ECF (Area(N) = Area(P)/(SF·TXFN·ECF)) · at%");
  const hv=1486.69;   // .vms 의 Al (mono) 1486.69
  const rows=regs.map(r=>({f:r.f,name:r.name,KE:hv-r.row["Peak BE"],ecf:r.row["Area (P) CPS.eV"]/(r.row["SF Al Scof"]*r.row["TXFN"]*r.row["Area (N) TPP-2M"]),row:r.row}));
  ok(rows.every(r=>r.row["TXFN"]===1),"모든 줄 TXFN = 1 — 이 자료로는 Avantage 투과 함수 형태를 정할 수 없음");
  const lsq=(xs,ys)=>{ const n=xs.length, mx_=xs.reduce((a,b)=>a+b)/n, my=ys.reduce((a,b)=>a+b)/n; let sxy=0,sxx=0; xs.forEach((x,i)=>{ sxy+=(x-mx_)*(ys[i]-my); sxx+=(x-mx_)*(x-mx_); }); const b=sxy/sxx; return {b,a:my-b*mx_,res:mx(xs.map((x,i)=>Math.abs(ys[i]-(my+b*(x-mx_)))))}; };
  const pw=lsq(rows.map(r=>Math.log(r.KE)),rows.map(r=>Math.log(r.ecf)));
  info("ECF ∝ KE^n 맞춤",`n = ${pw.b.toFixed(4)}, 최대 log 잔차 ${e1(pw.res)} (Area(N) 2 자리 반올림 한계 6.4e-4 보다 큼 → 거듭제곱 아님)`);
  const ratioDev=mat=>{ const l=rows.map(r=>Math.log(r.ecf/C.imfpValue({model:"tpp2m",mat},r.KE))); const m=l.reduce((a,b)=>a+b)/l.length; return mx(l.map(v=>Math.abs(v-m))); };
  const MAT_FIT={rho:6.22,Nv:10,M:183.3,Eg:0}, MAT_DEF={rho:5,Nv:24,M:150,Eg:3};
  ok(ratioDev(MAT_FIT)<7e-4,"ECF = 상수 × TPP-2M λ(KE = hν − Peak BE) — 물질 ρ 6.22, Nv/M 0.0546 (Eg 0) 로 반올림 한계 안",`최대 log 편차 ${e1(ratioDev(MAT_FIT))}`);
  info("INNO 기본 물질(ρ 5, Nv 24, M 150, Eg 3) 의 TPP-2M",`ECF 비 최대 log 편차 ${e1(ratioDev(MAT_DEF))}`);
  const atErr=imfp=>{ let m=0; ["1","2","3","4"].forEach(f=>{ const rr=rows.filter(r=>r.f===f); const q=C.quantify(rr.map(r=>({label:r.name,el:r.name,area:r.row["Area (P) CPS.eV"],S:r.row["SF Al Scof"],BE:r.row["Peak BE"],hv,PE:null,tf:{model:"none"}})),{imfp});
    q.forEach((x,i)=>{ m=Math.max(m,Math.abs(x.at-rr[i].row["Atomic %"])); }); }); return m; };
  ok(atErr({model:"tpp2m",mat:MAT_FIT})<0.01,"INNO quantify (Avantage Area(P)·SF, TPP-2M 맞춘 물질) at% = Avantage Atomic %",`최대 |Δ| ${atErr({model:"tpp2m",mat:MAT_FIT}).toFixed(4)} at%`);
  ok(atErr({model:"tpp2m",mat:MAT_DEF})<0.03,"INNO 기본 TPP-2M 물질로도 at% 차 < 0.03",`최대 |Δ| ${atErr({model:"tpp2m",mat:MAT_DEF}).toFixed(4)} at%`);
  info("KE^0.6 (INNO 기본 IMFP)",`최대 |Δ| ${atErr({model:"power",n:0.6}).toFixed(3)} at% — Avantage 'TPP-2M' 정량과 맞추려면 IMFP 를 TPP-2M 으로`);

  console.log("== A5 표 안의 일관성: Atomic % ↔ Area(N), PP At.% ↔ PP Hgt(N)");
  let dAt=0, dPP=0; ["1","2","3","4"].forEach(f=>{ const rr=regs.filter(r=>r.f===f).map(r=>r.row); const s=rr.reduce((a,r)=>a+r["Area (N) TPP-2M"],0), sp=rr.reduce((a,r)=>a+r["PP Hgt (N)"],0);
    rr.forEach(r=>{ dAt=Math.max(dAt,Math.abs(100*r["Area (N) TPP-2M"]/s-r["Atomic %"])); dPP=Math.max(dPP,Math.abs(100*r["PP Hgt (N)"]/sp-r["PP At. %"])); }); });
  ok(dAt<0.015&&dPP<0.015,"Atomic % = Area(N)/ΣArea(N), PP At.% = PP Hgt(N)/ΣPP Hgt(N)",`최대 ${dAt.toFixed(3)} · ${dPP.toFixed(3)} at%`);
  const k=rows.map(r=>r.row["PP Height CPS"]/(r.row["SF Al Scof"]*r.ecf*r.row["PP Hgt (N)"]));
  info("PP Hgt(N) = PP Height/(SF·ECF·k)",`k = 5.00 (1 줄 영역 ${k.filter((v,i)=>!/2p/.test(rows[i].name)).map(v=>v.toFixed(2)).filter((v,i,a)=>a.indexOf(v)===i).join("/")}), 2p 영역 ${k.filter((v,i)=>/2p/.test(rows[i].name)).map(v=>v.toFixed(2)).join(" ")} — 2p 의 k 는 미확인`);

  console.log("== A6 Height · Peak BE · FWHM (Avantage 'Standard' 봉우리) — 참고");
  const dh=[], db=[], df=[]; regs.forEach(r=>{ const E=r.sp.be.slice(r.i0,r.i1+1), T=Tof(r), nt=Float64Array.from(r.sp.y.slice(r.i0,r.i1+1),(v,i)=>v-T[i]); const pk=C.peakMaxPos(E,nt,E[0],E[E.length-1],0.15); if(!pk) return;
    dh.push(Math.abs(pk.y/r.row["Height CPS"]-1)); db.push(Math.abs(pk.be-r.row["Peak BE"])); });
  info("평활 최댓값(INNO peakMaxPos 0.15 eV) vs Avantage",`Height |Δ| 중앙 ${p2(med(dh))} 최대 ${p2(mx(dh))}; Peak BE |Δ| 중앙 ${med(db).toFixed(3)} eV 최대 ${mx(db).toFixed(3)} eV — Sr 3d 쌍에서 Avantage Peak BE 가 두 봉우리 사이(정의 미확인)`);
} else console.log("A) avantage_231107 없음 — 건너뜀");

/* ================= B) CasaXPS ================= */
const CS=path.join(ROOT,"casa_240213");
if(fs.existsSync(CS)){
  const parseExport=t=>{ const L=t.split(/\r?\n/); const H={}; let i=0; for(;i<L.length;i++){ const a=L[i].split("\t"); if(/^B\.E\./.test(a[0])) break; if(a[0]) H[a[0].trim()]=a.slice(2).filter(x=>x!==""); }
    const hdr=L[i].split("\t").filter(x=>x!==""); const rows=[]; for(i++;i<L.length;i++){ const a=L[i].split("\t").filter(x=>x!=="").map(Number); if(a.length>=hdr.length) rows.push(a); } return {H,rows}; };
  console.log("== B0 CasaXPS VAMAS 주석(CASA region / comp / Calib) 읽기");
  const V=await IO.readFile("Integrated.vms",rd(path.join(CS,"Integrated.vms")));
  const withM=V.spectra.filter(s=>s.meta.casaModel&&s.meta.casaModel.regions.length);
  ok(withM.length===6,"모형이 있는 블록 6 개 (Ba 4d · C 1s · O 1s × 2 시료)",withM.map(s=>s.meta.blockId).join(", "));
  const c0=withM[0].meta.casaModel;
  ok(c0.regions[0].bg==="Shirley"&&c0.regions[0].aw===2&&c0.comps.length===4&&c0.comps[1].area.link===0&&Math.abs(c0.comps[1].area.k-0.66666667)<1e-9&&c0.comps[1].fwhm.link===0&&c0.comps[1].pos.link===-1,
    "Ba 4d: Shirley, Av.Width 칸 2, 성분 4, 3/2 넓이 = 0.6667 × 5/2, 폭 연결, 위치 자유");
  const cases=[["Bare BCT (Ba4d) _ Fitting","[Bare] Ba4d.txt"],["Bare BCT (O1s) _ Fitting","[Bare] O1s.txt"],["D.I. on BCT (Ba4d) _ Fitting","[D.I.] Ba4d.txt"],["D.I. on BCT (O1s) _ Fitting","[D.I.] O1s.txt"]];
  const B1=[], G=[], LA=[], LA701=[], rms=[];
  for(const [fn,bid] of cases){
    const X=parseExport(fs.readFileSync(path.join(CS,fn),"latin1")); const sp=V.spectra.find(s=>s.meta.blockId===bid); const cm=sp.meta.casaModel; const nc=X.H.Name.length;
    const shiftX=X.rows[0][0]-sp.be[sp.n-1];
    ok(Math.abs(shiftX-cm.shift)<1e-6,`${fn}: 내보내기 BE 축 = 측정 BE + Calib 이동 (주석 Calib 줄의 합)`,`${shiftX.toFixed(4)} eV`);
    const ex=X.rows.map(r=>({be:r[0]-cm.shift,y:r[1],comp:r.slice(2,2+nc),bg:r[2+nc],env:r[3+nc]})).reverse();
    const M=C.casaModelToINNO(cm,sp.hv,{step:sp.step}); const P=C.makeProblem(sp,M); const E=P.E,y=P.y,n=E.length; const ctx={emin:E[0],emax:E[n-1]};
    const at=e=>ex.find(r=>Math.abs(r.be-e)<1e-3);
    const inr=ex.map(r=>Math.abs(r.bg-r.y)>1e-9); const a=inr.indexOf(true), b=inr.lastIndexOf(true);
    ok(Math.abs(E[1]-ex[a].be)<1e-6&&Math.abs(E[n-2]-ex[b].be)<1e-6,`${fn}: 변환된 창의 두 번째·끝에서 두 번째 채널 = CasaXPS 배경 구간 끝`,`${ex[a].be.toFixed(2)}–${ex[b].be.toFixed(2)} eV`);
    const B=C.bgStatic(E,y,P.bg); let mb=0,H=0; for(let i=1;i<n-1;i++){ mb=Math.max(mb,Math.abs(B[i]-at(E[i]).bg)); H=Math.max(H,y[i]-B[i]); } B1.push(mb/H);
    P.comps.forEach((c,j)=>{ const f=C.compEval(c,E,ctx); let md=0,m=0; for(let i=1;i<n-1;i++){ const r=at(E[i]); const v=r.comp[j]-r.bg; md=Math.max(md,Math.abs(f[i]-v)); m=Math.max(m,v); }
      (c.shape==="LA"?LA:G).push(md/m);
      if(c.shape==="LA"){ const c2=C.clone(c); c2.casaN=701; c2.F=c.F*C.casaWidthFactor("LA",{alpha:1,beta:1,casaN:c.casaN,wG:0})/C.casaWidthFactor("LA",{alpha:1,beta:1,casaN:701,wG:0}); const f2=C.compEval(c2,E,ctx); let d2=0; for(let i=1;i<n-1;i++){ const r=at(E[i]); d2=Math.max(d2,Math.abs(f2[i]-(r.comp[j]-r.bg))); } LA701.push(d2/m); } });
    const mu=new Float64Array(n); P.comps.forEach(c=>{ const f=C.compEval(c,E,ctx); for(let i=0;i<n;i++) mu[i]+=f[i]; });
    let sI=0,sX=0; for(let i=1;i<n-1;i++){ sI+=(y[i]-mu[i]-B[i])**2; sX+=(y[i]-at(E[i]).env)**2; }
    const P2=C.makeProblem(sp,M); const rr=C.fitOnce(P2,{initAreas:false}); const st=C.fitStats(P2,rr);
    rms.push({fn,casa:Math.sqrt(sX),inno0:Math.sqrt(sI),fit:st.casa.rms,areas:P2.comps.map((c,j)=>c.area/cm.comps[j].area.v-1)});
  }
  console.log("== B1 CasaXPS Shirley (주석 Av.Width 2 → ±1 채널, 영역 안 첫·끝 채널, 배타 누적합) vs 내보내기 Background");
  ok(mx(B1)<5e-4,"배경 최대|Δ|/봉우리 높이",B1.map(e1).join(" · "));
  console.log("== B2 Casa 매개변수에서 INNO 성분 곡선 vs 내보내기 성분 열");
  ok(mx(G)<5e-4,`GL(30) (${G.length} 성분): CasaXPS fwhm = 실제 FWHM 으로 F 환산`,`최대 ${e1(mx(G))} (fwhm 을 +1.7e-4 키우면 ≈ 2e-4 — SPEC_casa_v2 §2.2 와 같은 방향)`);
  ok(mx(LA)<2e-3&&mx(LA701)<2e-3,`LA(50) (${LA.length} 성분): LA(1,1,n), n = 700.5 (Fairley SI Eq. 6) 또는 701 (Cookbook p. 19) — 이 자료로는 구별 안 됨`,`700.5: ${e1(mx(LA))} · 701: ${e1(mx(LA701))}; 최소는 n ≈ 705–710 에서 1.1e-3 — SPEC_casa_v2 §4 의 Lorentzian 수치 차이`);
  console.log("== B3 포락선 RMS (√Σ(y − 포락선)², 채널 영역 안) · INNO 다시 맞춤");
  rms.forEach(r=>ok(Math.abs(r.inno0/r.casa-1)<0.01,`${r.fn}: Casa 매개변수에서 INNO RMS = Casa 내보내기 포락선 RMS`,`${r.inno0.toFixed(1)} / ${r.casa.toFixed(1)}`));
  rms.forEach(r=>ok(r.fit<=r.casa*1.0001,`${r.fn}: INNO(Casa 규약, 무가중) 다시 맞춤 RMS ≤ Casa 저장 모형 RMS`,`${r.fit.toFixed(1)} vs ${r.casa.toFixed(1)}; 넓이 변화 ${r.areas.map(v=>(100*v).toFixed(1)+"%").join(" ")}`));
  info("결론","저장된 CasaXPS 모형은 최소제곱 최소점이 아니다(INNO 가 같은 목적함수에서 더 낮은 RMS) → 넓이·폭 비교는 CasaXPS 에서 끝까지 맞춘 출력이 있어야 가능");
  const VA=path.join(CS,"various");
  if(fs.existsSync(VA)){
    console.log("== B4 D.I. Ba 4d 다른 맞춤 내보내기 (GL(30) #1 · LA(50) #1) — 배경");
    const sp=V.spectra.find(s=>s.meta.blockId==="[D.I.] Ba4d.txt");
    for(const nm of ["GL(30)","LA(50)"]){ const X=parseExport(fs.readFileSync(path.join(VA,"[D.I.] Ba4d _ "+nm+" #1"),"latin1")); const nc=X.H.Name.length; const sh=X.rows[0][0]-sp.be[sp.n-1];
      const ex=X.rows.map(r=>({be:r[0]-sh,y:r[1],bg:r[2+nc],env:r[3+nc]})).reverse(); const inr=ex.map(r=>Math.abs(r.bg-r.y)>1e-9); const a=inr.indexOf(true), b=inr.lastIndexOf(true);
      const st=C.makeProblem(sp,{lo:ex[a].be-sp.step-1e-6,hi:ex[b].be+sp.step+1e-6,bg:{type:"shirley",conv:"casa",aw:1},comps:[],cons:{}}); const B=C.bgStatic(st.E,st.y,st.bg);
      let m=0,H=0; for(let i=1;i<st.E.length-1;i++){ const r=ex.find(r=>Math.abs(r.be-st.E[i])<1e-3); m=Math.max(m,Math.abs(B[i]-r.bg)); H=Math.max(H,st.y[i]-B[i]); }
      ok(m/H<5e-4,`${nm} #1: 배경 (영역 ${ex[a].be.toFixed(2)}–${ex[b].be.toFixed(2)} eV, Calib ${sh.toFixed(4)})`,e1(m/H)); }
  }
} else console.log("B) casa_240213 없음 — 건너뜀");

/* ================= C) CasaXPS Cookbook 2026 정의 (p. 67) ================= */
console.log("== C GL(m)T(k) = (GL + (1 − GL)·E(1,k)), E(1,k) = exp(k·x) (x < 0 = 높은 BE, x 단위 F) — INNO AGT(tu 'fwhm', TM 0)");
{ const F=1.2, m=0.3, k=2; let d=0; for(let x=-3;x<=3;x+=0.01){ const t=x/F, gl=Math.exp(-4*Math.LN2*(1-m)*t*t)/(1+4*m*t*t), cb=gl+(x>0?(1-gl)*Math.exp(-k*t):0);
    d=Math.max(d,Math.abs(C.aglRaw(x,{F,m,ET:k,TM:0,CT:0,tu:"fwhm"})-cb)); }
  ok(d<1e-12,"INNO aglRaw = Cookbook 식",e1(d)); }

console.log(`\n${pass} 통과 · ${fail} 실패`); if(fail){ console.log(fails.map(s=>" - "+s).join("\n")); process.exit(1); }
})().catch(e=>{ console.error(e); process.exit(2); });
