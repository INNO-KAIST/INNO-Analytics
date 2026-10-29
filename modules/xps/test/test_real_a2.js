/* 실측 대조 a2: Thermo K-Alpha 실측 2 시료의 Avantage 출력 ↔ INNO Avantage 규약 (test/ref/avantage_a2_patch.md)
   node test/test_real_a2.js            (modules/xps 에서, 약 5 초)
   자료: test/data/real/a2/<시료 폴더 2 개>/ (비공개 — 공개 소스 백업에서 빠짐). 폴더가 없으면 모든 항목을 건너뛰고 0 으로 끝난다.
     Group/<Region> Scan.VGD · Survey.VGD · "[1] … _ Avantage.xlsx" (PeakTable_Smart / PeakTable_Shirley + <Region>_Smart / <Region>_Shirley 'Backgnd.' 곡선)
   Avantage 설정(사용자): Smart · Shirley 모두 Background Average at Start&End 0.50 eV, Shirley Iterations 50 · Convergence 0.01,
     K-Alpha Al Kα mono, PE 50, step 0.05 eV, TXF 는 측정 때 제거('TXF Shape Correction, retaining Intensity') → TXFN = 1.
   코어: v1.9.2 부터 패치가 src 에 합쳐졌으므로 src/xps_core.js 를 그대로 쓴다. 먼저 avantage_a2_apply.py --check 로
     네 단계가 모두 들어 있는지(아무것도 할 것이 없는지) 확인하고, 아니면 실패로 끝낸다. 다른 코어를 쓰려면 A2_CORE=/path/to/xps_core.js.   */
const path=require("path"), fs=require("fs"), cp=require("child_process");
const IO=require("../src/xps_io.js");
const ROOT=path.join(__dirname,"data","real","a2");
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
function info(name,v){ console.log("  info "+name+(v?"  ("+v+")":"")); }
const med=a=>{ const s=a.slice().sort((x,y)=>x-y); return s.length?s[s.length>>1]:NaN; };
const mx=a=>a.reduce((m,v)=>Math.max(m,v),-Infinity);
const e1=v=>v.toExponential(1), p3=v=>(100*v).toFixed(3)+" %";
const rd=f=>new Uint8Array(fs.readFileSync(f));

if(!fs.existsSync(ROOT)){ console.log("test/data/real/a2 없음 — 실측 대조를 건너뜀 (skip)"); process.exit(0); }

/* ---- 코어: src (패치가 합쳐져 있어야 함) ---- */
function liveCore(){
  const file=process.env.A2_CORE?path.resolve(process.env.A2_CORE):path.join(__dirname,"..","src","xps_core.js");
  const r=cp.spawnSync("python3",["-I",path.join(__dirname,"ref","avantage_a2_apply.py"),"--check",file],{encoding:"utf8"});
  const applied=r.status===0&&/nothing to do/.test(r.stdout||"");
  return {file,applied,note:(process.env.A2_CORE?"A2_CORE, ":"src, ")+(applied?"avantage_a2 패치 4 단계 모두 들어 있음":"패치 확인 실패: "+((r.stdout||"")+(r.stderr||"")).trim().split("\n").slice(-2).join(" / "))};
}
const CP=liveCore(); const C=require(CP.file);
console.log("코어: "+CP.file+"  ("+CP.note+")");
ok(CP.applied,"코어에 avantage_a2 패치가 한 번씩 모두 들어 있음 (avantage_a2_apply.py --check: nothing to do)",CP.note);

(async()=>{
const AVG_EV=0.5, HV=1486.68;
const MAT={rho:6.22,Nv:10,M:183.3,Eg:0};   // SPEC_real_validation §A3 (Kratos 자료에서 맞춘 TPP-2M 물질) — a2 에서도 같음
const sets=fs.readdirSync(ROOT).filter(d=>fs.statSync(path.join(ROOT,d)).isDirectory()).sort();
const regs=[];   // {set, bg, name, sp(VGD), i0, i1, T(full-length), row}
const tables={};
console.log("== A0 읽기: VGD 세기 = xlsx counts/s, 배경 곡선의 0 아닌 구간 = Start/End BE");
for(const set of sets){
  const dir=path.join(ROOT,set), xf=fs.readdirSync(dir).find(f=>/Avantage\.xlsx$/i.test(f)); if(!xf){ ok(false,set+": Avantage xlsx 없음"); continue; }
  const X=await IO.readFile(xf,rd(path.join(dir,xf)));
  for(const bg of ["Smart","Shirley"]){
    const pt=X.tables.find(t=>t.name==="PeakTable_"+bg); const H=pt.rows[0].map(v=>String(v).trim()); const col=n=>H.indexOf(n);
    const rows=pt.rows.slice(1).filter(r=>typeof r[col("Start BE")]==="number").map(r=>Object.fromEntries(H.map((h,k)=>[h,r[k]])));
    tables[set+"/"+bg]=rows;
    for(const row of rows){
      const nm=String(row.Name).trim(); const xs=X.spectra.find(s=>s.name===nm+"_"+bg); const vf=path.join(dir,"Group",nm+" Scan.VGD");
      if(!xs||!fs.existsSync(vf)){ ok(false,`${set} ${bg} ${nm}: 곡선 시트 또는 VGD 없음`); continue; }
      const sp=(await IO.readFile(nm+" Scan.VGD",rd(vf))).spectra[0];
      const bcol=xs.extra.avantageCols[Object.keys(xs.extra.avantageCols).find(k=>/backgnd/i.test(k))];
      const T=new Float64Array(sp.n); let dy=0;
      for(let k=0;k<xs.n;k++){ let b=0; for(let i=0;i<sp.n;i++) if(Math.abs(sp.be[i]-xs.be[k])<Math.abs(sp.be[b]-xs.be[k])) b=i; T[b]=bcol[k]; dy=Math.max(dy,Math.abs(sp.y[b]-xs.y[k])/Math.max(1,Math.abs(xs.y[k]))); }
      let i0=-1,i1=-1; for(let i=0;i<sp.n;i++) if(T[i]!==0){ if(i0<0) i0=i; i1=i; }
      const near=e=>{ let b=0; for(let i=0;i<sp.n;i++) if(Math.abs(sp.be[i]-e)<Math.abs(sp.be[b]-e)) b=i; return b; };
      regs.push({set,bg,name:nm,sp,i0,i1,T,row,dy,rangeOK:near(row["End BE"])===i0&&near(row["Start BE"])===i1});
    }
  }
}
const nR=regs.length;
ok(nR===24,"영역 24 개 (2 시료 × 6 영역 × Smart/Shirley)",""+nR);
ok(regs.every(r=>r.dy<1e-5),"VGD 세기 = xlsx counts/s (상대 1e-5 안)","최대 "+e1(mx(regs.map(r=>r.dy))));
ok(regs.every(r=>r.rangeOK),"Backgnd 0 아닌 구간 = Start/End BE 에 가장 가까운 채널",`${regs.filter(r=>r.rangeOK).length}/${nR}`);
{ const tf=regs[0].sp.meta.tf; info("VGD 투과 계수 (meta.tf, 이 자료에서는 Avantage 가 쓰지 않음)",tf?tf.coef.map(v=>v.toPrecision(4)).join(", "):"없음");
  // VGD prop 415/416 (Acq_TXFN_Correction · 모드) → xps_io 의 meta.tfRemoved · tfMode (v1.9.2)
  ok(regs.every(r=>r.sp.meta.tfRemoved===true&&r.sp.meta.tfMode===3),"VGD prop 415 (TXFN 보정 켜짐) = true, 416 (모드) = 3 → meta.tfRemoved (모든 영역) — 정량에서 TF 를 다시 곱하면 안 됨",`${regs.filter(r=>r.sp.meta.tfRemoved===true).length}/${regs.length}`);
  ok(regs.every(r=>r.row.TXFN===1),"Peak Table TXFN = 1 (모든 줄)"); }

const net=(r,B)=>{ const E=r.sp.be, y=r.sp.y; let A=0,Hh=0; for(let i=r.i0;i<=r.i1;i++){ if(i>r.i0) A+=0.5*((y[i]-B[i])+(y[i-1]-B[i-1]))*Math.abs(E[i]-E[i-1]); Hh=Math.max(Hh,y[i]-B[i]); } return {A,H:Hh}; };
const half=r=>C.avantageHalf({avgEV:AVG_EV},r.sp.be.slice(r.i0,r.i1+1));

console.log("== A1 끝값 = 끝 채널 ±h 평균 (h = round(0.50 eV / (2·0.05 eV)) = 5, 창 밖 채널 포함, 스펙트럼 끝에서 자름)");
{ let n=0, w=0; regs.forEach(r=>{ const h=half(r), y=r.sp.y, N=r.sp.n; const av=c=>{ let s=0,k=0; for(let j=Math.max(0,c-h);j<=Math.min(N-1,c+h);j++){ s+=y[j]; k++; } return s/k; };
    [r.i0,r.i1].forEach(c=>{ const d=Math.abs(av(c)-r.T[c]); w=Math.max(w,d/Math.abs(r.T[c])); if(d<=Math.max(0.06,6e-6*Math.abs(r.T[c]))) n++; }); });
  ok(n===2*nR,"끝값 48/48 이 ±5 채널 평균 (인쇄 자리수 안)",`${n}/${2*nR}, 최대 상대차 ${e1(w)}`);
  ok(regs.every(r=>half(r)===5),"avantageHalf({avgEV:0.5}) = 5 (step 0.05 eV)"); }

console.log("== A2 Area(P) = 내보낸 배경 위 사다리꼴 적분 (창 채널)");
{ const d=regs.map(r=>Math.abs(net(r,r.T).A/r.row["Area (P) CPS.eV"]-1)), w=regs[d.indexOf(mx(d))];
  const g=regs.map(r=>{ let G=0; const E=r.sp.be,y=r.sp.y; for(let i=r.i0+1;i<=r.i1;i++) G+=0.5*Math.abs((y[i]-r.T[i])+(y[i-1]-r.T[i-1]))*Math.abs(E[i]-E[i-1]); return Math.abs(net(r,r.T).A-r.row["Area (P) CPS.eV"])/G; });
  ok(mx(g)<5e-4,"Area(P) = ∫(y − Backgnd) 사다리꼴 (24/24; 차/∫|y − B| < 5e-4)",`중앙 ${e1(med(g))} · 최대 ${e1(mx(g))}`);
  info("순 넓이 상대차",`중앙 ${e1(med(d))} · 최대 ${e1(mx(d))} (${w.set} ${w.bg} ${w.name} — 양·음 상쇄로 순 넓이가 작음)`); }
console.log("== A3 Smart 접점 = 평활 S(±h, 창 안, 낮은 BE 끝 채널 제외)");
{ let nT=0, nAbove=0, nTouch=0; regs.filter(r=>r.bg==="Smart").forEach(r=>{ const y=Array.from(r.sp.y.slice(r.i0,r.i1+1)), S=C.avantageSmooth(y,half(r));
    for(let k=1;k<y.length-1;k++){ const d=r.T[r.i0+k]-S[k]; nT++; if(Math.abs(d)<=0.06+6e-6*S[k]) nTouch++; else if(d>0) nAbove++; } });
  info("Smart 내부 채널: 접점(배경 = S) · S 위로 지나는 채널",`${nTouch} · ${nAbove} / ${nT} (위로 지나는 채널은 모두 창 끝 근처 — 패치 md §2.5)`);
  ok(nAbove<=20,"Avantage Smart 가 S 위로 지나는 채널 ≤ 20 (12 영역 합)",""+nAbove); }

const runBg=(r,type,extra)=>{ const E=r.sp.be.slice(r.i0,r.i1+1), y=r.sp.y.slice(r.i0,r.i1+1);
  const bg=Object.assign({type,conv:"avantage",avgEV:AVG_EV},extra||{}); bg._pad=C.avantagePad(r.sp.y,r.i0,r.i1,C.avantageHalf(bg,E));
  const b=C.bgStatic(E,y,bg); const B=new Float64Array(r.sp.n); for(let k=0;k<b.length;k++) B[r.i0+k]=b[k]; return {B,b}; };
const cmp=(r,B)=>{ const ref=net(r,r.T); let m=0; for(let i=r.i0;i<=r.i1;i++) m=Math.max(m,Math.abs(B[i]-r.T[i])); return {dB:m/ref.H,dA:net(r,B).A/r.row["Area (P) CPS.eV"]-1}; };
const modelArea={};
for(const bg of ["Smart","Shirley"]){
  console.log(`== A4 ${bg}: INNO conv:"avantage" (avgEV 0.5) vs Avantage — 최대|ΔB|/봉우리 높이 · Area(P) 상대차`);
  const rs=regs.filter(r=>r.bg===bg), res=[];
  rs.forEach(r=>{ const {B,b}=runBg(r,bg==="Smart"?"smart":"shirley"); const c=cmp(r,B); res.push(Object.assign({r,b},c)); modelArea[r.set+"/"+bg+"/"+r.name]=net(r,B).A;
    console.log(`        ${r.set.padEnd(8)} ${r.name.padEnd(5)} ΔB/H ${e1(c.dB).padStart(8)}  ΔA ${p3(c.dA).padStart(9)}`+(b.splits?`  접점 ${b.splits.length}`:"")+(b.aboveFrac!=null?`  S(±2h) 위 ${(100*b.aboveFrac).toFixed(0)} %`:"")+(b.fallback?"  직선":"")); });
  const good=res.filter(x=>x.dB<=1.5e-3&&Math.abs(x.dA)<=1e-3);
  info(`${bg} |ΔA| 중앙 · 최대`,`${p3(med(res.map(x=>Math.abs(x.dA))))} · ${p3(mx(res.map(x=>Math.abs(x.dA))))}; ΔB/H 중앙 ${e1(med(res.map(x=>x.dB)))}`);
  if(bg==="Smart"){
    ok(good.length>=10,"Smart: ΔB/H ≤ 1.5e-3 그리고 |ΔA| ≤ 0.1 % 인 영역 ≥ 10/12",`${good.length}/12 (남는 것: 두 시료의 O1s ΔB — md §2.5)`);
    ok(mx(res.map(x=>Math.abs(x.dA)))<=3e-3,"Smart: 모든 영역 |ΔA| ≤ 0.3 %",p3(mx(res.map(x=>Math.abs(x.dA)))));
  } else {
    const co=res.reduce((m,x)=>Math.abs(x.dA)>Math.abs(m.dA)?x:m), rest=res.filter(x=>x!==co);   // 병적 영역 = |ΔA| 최대 (a2: S2 Co3p, md §4)
    ok(rest.every(x=>x.dB<=1e-3&&Math.abs(x.dA)<=1e-3),"Shirley: |ΔA| 최대 영역(병적, md §4)을 뺀 11 영역 ΔB/H ≤ 1e-3 · |ΔA| ≤ 0.1 %",`${rest.filter(x=>x.dB<=1e-3&&Math.abs(x.dA)<=1e-3).length}/11`);
    if(co){ ok(Math.abs(co.dA)<=0.035&&co.b.aboveFrac>0.5,`Shirley ${co.r.name} (병적: 배경이 평활 데이터 위 > 50 %) — |ΔA| ≤ 3.5 %, 경고 플래그 aboveFrac > 0.5`,`ΔA ${p3(co.dA)}, aboveFrac ${co.b.aboveFrac.toFixed(2)}`); }
  }
}

console.log("== A5 방향 · 옛 규약 · makeProblem");
{ let w=0; regs.forEach(r=>{ const E=r.sp.be.slice(r.i0,r.i1+1), y=r.sp.y.slice(r.i0,r.i1+1); const bg={type:r.bg==="Smart"?"smart":"shirley",conv:"avantage",avgEV:AVG_EV}; bg._pad=C.avantagePad(r.sp.y,r.i0,r.i1,5);
    const a=C.bgStatic(E,y,bg), b=C.bgStatic(E.slice().reverse(),y.slice().reverse(),Object.assign({},bg,{_pad:{lo:bg._pad.hi.slice().reverse(),hi:bg._pad.lo.slice().reverse()}}));
    for(let k=0;k<a.length;k++) w=Math.max(w,Math.abs(a[k]-b[a.length-1-k])); });
  ok(w<1e-6,"BE 내림차순 배열도 같은 배경 (뒤집기 불변)",e1(w)); }
{ const r=regs.find(q=>q.bg==="Smart"); const st=C.makeProblem({be:r.sp.be,y:r.sp.y,c:null},{lo:r.sp.be[r.i0],hi:r.sp.be[r.i1],bg:{type:"smart",conv:"avantage",avgEV:AVG_EV},comps:[],cons:{}});
  const st2=C.makeProblem({be:r.sp.be,y:r.sp.y,c:null},{lo:r.sp.be[r.i0],hi:r.sp.be[r.i1],bg:{type:"shirley",conv:"avantage",avgEV:AVG_EV},comps:[],cons:{}});
  ok(st.bg._pad&&st2.bg._pad&&st.bg._pad.lo.length===Math.min(5,r.i0)&&st.bg._pad.hi.length===Math.min(5,r.sp.n-1-r.i1),"makeProblem 이 Smart · Shirley 모두에 창 밖 min(h, 남은 채널) 을 넘김",`lo ${st.bg._pad.lo.length}, hi ${st.bg._pad.hi.length}`);
  const E=r.sp.be.slice(r.i0,r.i1+1), y=r.sp.y.slice(r.i0,r.i1+1), pad=C.avantagePad(r.sp.y,r.i0,r.i1);
  const v18=C.bgStatic(E,y,{type:"smart",conv:"avantage",rule:"v18",_pad:pad}), dflt=C.bgStatic(E,y,{type:"smart",conv:"avantage",_pad:pad}), dsh=C.bgStatic(E,y,{type:"shirley",conv:"avantage",n:5,it:10,cl:0});
  ok(v18.rule==="v18"&&dflt.rule==="v18"&&dsh.rule!=="v19","옛 규약(v1.8)은 기본(avgEV · ah 없음)과 bg.rule \"v18\" 로 그대로 — 기존 프로젝트·SPEC §2c 시험 불변"); }

console.log("== A6 Area(N) · Atomic % (Area(N) = Area(P)/(SF·TXFN·ECF), ECF = TPP-2M λ(KE = hν − Peak BE))");
{ let dev=0, self=0; for(const k of Object.keys(tables)){ const rows=tables[k]; const rat=rows.map(r=>r["Area (P) CPS.eV"]/(r.SF*r.TXFN*r["Area (N)"])/C.tpp2m(HV-r["Peak BE"],MAT)); const m=rat.reduce((a,b)=>a+b)/rat.length; dev=Math.max(dev,mx(rat.map(v=>Math.abs(Math.log(v/m)))));
    const s=rows.reduce((a,r)=>a+r["Area (N)"],0); self=Math.max(self,mx(rows.map(r=>Math.abs(100*r["Area (N)"]/s-r["Atomic %"])))); }
  ok(dev<5e-4,"ECF = 상수 × TPP-2M λ (ρ 6.22, Nv 10, M 183.3, Eg 0 — Kratos 자료와 같은 물질)",`최대 log 편차 ${e1(dev)}`);
  ok(self<0.01,"Atomic % = Area(N)/ΣArea(N)",`최대 ${self.toFixed(4)} at%`); }
const atErr=(k,tf,areaOf)=>{ const rows=tables[k]; const q=C.quantify(rows.map(r=>({label:r.Name,el:r.Name,area:areaOf(r),S:r.SF,BE:r["Peak BE"],hv:HV,PE:50,tf})),{imfp:{model:"tpp2m",mat:MAT}}); return mx(q.map((x,i)=>Math.abs(x.at-rows[i]["Atomic %"]))); };
{ const tfV=regs[0].sp.meta.tf; let e0=0, eM=0, eT=0, eTa=0;
  for(const k of Object.keys(tables)){ const [set,bg]=k.split("/");
    e0=Math.max(e0,atErr(k,{model:"none"},r=>r["Area (P) CPS.eV"]));
    eM=Math.max(eM,atErr(k,{model:"none"},r=>modelArea[set+"/"+bg+"/"+String(r.Name).trim()]));
    eT=Math.max(eT,atErr(k,tfV,r=>r["Area (P) CPS.eV"])); eTa=Math.max(eTa,atErr(k,Object.assign({},tfV,{model:"avantage"}),r=>r["Area (P) CPS.eV"])); }
  ok(e0<0.01,"INNO quantify (Avantage Area(P), SF, TF 없음, TPP-2M) at% = Avantage Atomic %",`최대 |Δ| ${e0.toFixed(4)} at%`);
  ok(eM<0.1,"INNO quantify (INNO Avantage 규약 Area(P)) at% ≈ Avantage",`최대 |Δ| ${eM.toFixed(3)} at% (병적 Shirley 영역 넓이 차 포함)`);
  ok(eT>1,"VGD 투과 계수를 다시 곱하면(UI tfFor 'auto' 경로) at% 가 크게 틀어짐 — 이중 보정 확인",`최대 |Δ| ${eT.toFixed(2)} at% (thermo 형태), ${eTa.toFixed(2)} at% (Avantage log10 형태)`); }

console.log("== A7 Kratos 자료(test/data/real/kara, 있을 때): v1.8 규약 vs v1.9.2 규약(h = 4, bg.ah 4) — Smart Area(P)");
{ const AV=path.join(__dirname,"data","real","kara","avantage_231107");
  if(!fs.existsSync(AV)) info("kara 자료 없음 — 건너뜀");
  else { const res={v18:[],v19:[]};
    for(const f of ["1","2","3","4"]){ const X=await IO.readFile(f+".xls",rd(path.join(AV,f+".xls")));
      const pt=X.tables.find(t=>/peak\s*table/i.test(t.name)); const hi=pt.rows.findIndex(r=>r.some(v=>/^name$/i.test(String(v).trim()))); const H=pt.rows[hi].map(v=>String(v).trim()); const col=n=>H.indexOf(n);
      pt.rows.slice(hi+1).filter(r=>typeof r[col("Start BE")]==="number").forEach(r=>{ const sp=X.spectra.find(s=>s.name===String(r[col("Name")]).trim()); if(!sp) return;
        const T=sp.extra.avantageCols[Object.keys(sp.extra.avantageCols).find(k=>/backgnd/i.test(k))], st=r[col("Start BE")], en=r[col("End BE")];
        let i0=0,i1=0; for(let i=0;i<sp.n;i++){ if(Math.abs(sp.be[i]-en)<Math.abs(sp.be[i0]-en)) i0=i; if(Math.abs(sp.be[i]-st)<Math.abs(sp.be[i1]-st)) i1=i; }
        const E=sp.be.slice(i0,i1+1), y=sp.y.slice(i0,i1+1); const A=B=>{ let a=0; for(let k=1;k<E.length;k++) a+=0.5*((y[k]-B[k])+(y[k-1]-B[k-1]))*(E[k]-E[k-1]); return a; }; const Aref=A(T.slice(i0,i1+1));
        res.v18.push(Math.abs(A(C.bgStatic(E,y,{type:"smart",conv:"avantage",_pad:C.avantagePad(sp.y,i0,i1,4)}))/Aref-1));
        res.v19.push(Math.abs(A(C.bgStatic(E,y,{type:"smart",conv:"avantage",ah:4,_pad:C.avantagePad(sp.y,i0,i1,4)}))/Aref-1)); }); }
    info("v1.8 |ΔA| 중앙 · 최대",`${p3(med(res.v18))} · ${p3(mx(res.v18))}`); info("v1.9.2 (ah 4) |ΔA| 중앙 · 최대",`${p3(med(res.v19))} · ${p3(mx(res.v19))}`);
    ok(res.v19.length===26&&med(res.v19)<med(res.v18)/3&&mx(res.v19)<0.05,"Kratos 26 영역: v1.9.2 규약(h = 4)이 v1.8 보다 중앙값 3 배 이상 작고 최대 < 5 %",`${p3(med(res.v19))} vs ${p3(med(res.v18))}`); } }

console.log(`\n${pass} 통과 · ${fail} 실패`);
if(fail){ console.log("실패:\n  "+fails.join("\n  ")); process.exit(1); }
})().catch(e=>{ console.error(e); process.exit(1); });
