/* CasaXPS 2.3.27 Round-3 실측 대조 (SPEC_casa_v3.md · test/ref/casa_readings_v3.md · test/ref/ref_casa_v3.py → ref_casa_v3.json)
   node test/test_casa_v3.js   (modules/xps 에서, 약 1 분)
   1) INNO CasaXPS 규약 피팅(weight "casa", Linear conv casa Av.W 1, 영역 280.03/297.13) 으로 K_*.vms 및 V2 행을 맞춤 →
      (a) 독립 scipy 기준(ref_casa_v3.json) 과 넓이·위치·폭·RMS, (b) CasaXPS 판독과 넓이·fwhm·위치·RMS·DoF·asym
   2) 표기 왕복: compToCasa(맞춘 성분) 표기 = CasaXPS 문자열; casaToComp(문자열, CasaXPS pos·fwhm) → 맞춘 INNO 성분 (폭·위치)
   3) casaModelToINNO: DS fwhm 연결 비 = wL, 위치 한계는 CasaXPS 위치 − posShift */
const path=require("path"), fs=require("fs");
const D=p=>path.join(__dirname,p);
const C=require("../src/xps_core.js"), IO=require("../src/xps_io.js");
const REF=JSON.parse(fs.readFileSync(D("ref/ref_casa_v3.json"),"utf8"));
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
const rel=(a,b)=>Math.abs(a-b)/Math.abs(b), e1=x=>x.toExponential(1);
const I0=19, I1=361;
// 판정 기준 (SPEC_casa_v3 §1·§3, casa_v3_patch §6): CasaXPS fwhm 허용 — 꼭대기 매끈(정함) / 경험 δ / ST(미정)
const KIND=s=>/^ST/.test(s)?"ST":/T\(/.test(s)?"tail":/^LA\(([\d.]+),([\d.]+),0\)/.test(s)&&(()=>{ const m=s.match(/^LA\(([\d.]+),([\d.]+)/); return m[1]!==m[2]; })()?"laKink":"smooth";
const FW_TOL={smooth:5e-4,laKink:1e-3,tail:1e-3,ST:5e-3}, ASYM_TOL={smooth:1.2e-3,laKink:2e-3,tail:1.2e-3,ST:3e-4};

(async()=>{
  const rd=async f=>(await IO.readFile(f+".vms",new Uint8Array(fs.readFileSync(D("data/casaval/"+f+".vms"))))).spectra[0];
  console.log("== 1 CasaXPS 규약 피팅 — 기준 구현(scipy) 과 CasaXPS 2.3.27 판독");
  const fitted={};
  for(const t of REF.fits){ const k=t.read, s=await rd(t.file);
    // 시작값: 판독 fwhm 을 조금 흔든 값에서 (casaToComp 로 F 환산) — 맞춤이 판독값을 그대로 받지 않도록
    const r0=C.casaToComp(t.casa,{be:285.05,area:12000,fwhm:k.fwhm*1.06,strict:true}); const comp=r0.comp; comp.id="a";
    const cons={a:{}}; C.paramNames(comp).forEach(q=>{ if(["be","area","F"].includes(q)||(comp.shape==="DS"&&q==="wL")||C.derivedParam(comp,q)) return; cons.a[q]={mode:"fixed"}; });
    const P=C.makeProblem({be:s.be,y:s.y,c:1},{lo:s.be[I0],hi:s.be[I1],bg:{type:"linear",conv:"casa",aw:1},comps:[comp],cons,weight:"casa"});
    const r=C.fitOnce(P,{}); const st=C.fitStats(P,r); const c=P.comps[0], cz=C.compToCasa(c), W=c.shape==="DS"?c.wL:c.F; fitted[t.file]={c,cz,st};
    // (a) 기준 구현: 같은 정의·같은 CasaXPS 피팅 규약 → 숫자 일치 (DS 넓이 = INNO 창 넓이 ≠ 기준의 진폭 매개변수 → 넓이 제외)
    const aOK=t.area==null||rel(c.area,t.area)<1e-4, rOK=t.rms<1?st.casa.rms<0.05:rel(st.casa.rms,t.rms)<1e-4;
    ok(aOK&&Math.abs(c.be-t.pos)<1e-4&&rel(W,t.F)<1e-4&&rOK&&st.casa.dof===t.dof_rule,`${t.file} ${t.casa}: INNO CasaXPS-mode fit = scipy reference (area, pos, width, RMS ≤ 1e-4; DoF rule)`,
      `${t.area==null?"area — (DS)":"area "+e1(rel(c.area,t.area))} · pos ${e1(Math.abs(c.be-t.pos))} eV · ${c.shape==="DS"?"wL":"F"} ${e1(rel(W,t.F))} · RMS ${st.casa.rms.toFixed(3)} / ${t.rms.toFixed(3)} · DoF ${st.casa.dof}/${t.dof_rule}`);
    // (b) CasaXPS 판독
    const kind=KIND(t.casa), sh=cz.posShift, cpos=c.be+(sh||0);
    const vs=`area ${c.area.toFixed(1)} vs ${k.area}, fwhm ${cz.fwhm!=null?cz.fwhm.toFixed(6):"—"} vs ${k.fwhm} (${cz.fwhm!=null?e1(rel(cz.fwhm,k.fwhm)):"—"}), pos ${cpos.toFixed(4)}${sh?` (INNO ${c.be.toFixed(4)} + ${sh.toFixed(4)})`:""} vs ${k.pos}, RMS ${st.casa.rms.toFixed(3)} vs ${k.rms}, DoF ${st.casa.dof} vs ${k.dof}, asym ${cz.asym!=null?cz.asym.toFixed(4):"—"} vs ${k.asym}`;
    ok(st.casa.dof===k.dof,`${t.file}: CasaXPS DoF = #{channels with y ≠ background} − 3`,`${st.casa.dof} vs ${k.dof}`);
    const areaOK=t.file.startsWith("K_DS")?true:rel(c.area,k.area)<1e-3;
    const posOK=kind==="ST"?true:Math.abs(cpos-k.pos)<2e-4, rmsOK=rel(st.casa.rms,k.rms)<0.035||Math.abs(st.casa.rms-k.rms)<20;   // CasaXPS RMS floor on exact data: up to 22.9 (LA(2,3,0)), SPEC_casa_v3 §2
    const fwOK=cz.fwhm!=null&&rel(cz.fwhm,k.fwhm)<FW_TOL[kind], asOK=cz.asym!=null&&Math.abs(cz.asym-k.asym)<ASYM_TOL[kind];
    const exOK=kind==="smooth"?cz.fwhmExact===true:cz.fwhmExact===false;
    ok(areaOK&&posOK&&rmsOK&&fwOK&&asOK&&exOK&&cz.notation===t.casa,`${t.file} vs CasaXPS ${t.casa}: ${t.file.startsWith("K_DS")?"":"area ≤ 1e-3, "}${kind==="ST"?"":"pos ≤ 2e-4 eV, "}RMS ≤ 3.5 % or 20 CPS, fwhm ≤ ${FW_TOL[kind]}, asym ≤ ${ASYM_TOL[kind]} (${kind==="smooth"?"exact":"empirical — marked approximate"})`,vs);
    if(t.file.startsWith("K_DS")) console.log(`     not pinned: DS area — CasaXPS ${k.area} / INNO window area ${c.area.toFixed(1)} = ${(k.area/c.area).toFixed(4)} (SPEC_casa_v3 §4)`);
    if(kind==="ST") console.log(`     not pinned: ST position — CasaXPS ${k.pos} − INNO origin ${c.be.toFixed(4)} = ${(k.pos-c.be).toFixed(4)} eV; fwhm ${cz.fwhm.toFixed(6)} (true FWHM of the ST shape) vs ${k.fwhm} (SPEC_casa_v3 §5)`);
  }
  console.log("== 1b asym index only (V2_GL30T2 data fitted as GL(30)T(1), GL(0)T(2) in CasaXPS): empirical δ = 0.011 F");
  for(const t of REF.shape_only){ const c=C.casaToComp(t.shape,{fwhm:1,strict:true}).comp; c.F=1; const cz=C.compToCasa(c);
    ok(Math.abs(cz.asym-t.asym_read)<2e-3&&cz.fwhmExact===false,`${t.shape}: asym ${cz.asym.toFixed(4)} vs CasaXPS ${t.asym_read} (≤ 2e-3, empirical)`,`fwhm/F ${cz.fwhm.toFixed(5)} (ref δ-implied ${t.w_with_delta.toFixed(5)})`); }

  console.log("== 2 표기 왕복: INNO 성분 → CasaXPS (표기·fwhm·pos) → INNO 성분");
  for(const t of REF.fits){ const {c,cz}=fitted[t.file], k=t.read, kind=KIND(t.casa);
    const back=C.casaToComp(cz.notation,{be:cz.posShift!=null?c.be+cz.posShift:c.be,area:c.area,fwhm:cz.fwhm,strict:true}).comp, W=c.shape==="DS"?"wL":"F";
    const same=["shape","mod","tail","casaN","stMu","stN","alpha","beta","mu","m","tk","tj"].every(q=>(c[q]==null&&back[q]==null)||c[q]===back[q]||Math.abs((+c[q]||0)-(+back[q]||0))<1e-9);
    ok(same&&rel(back[W],c[W])<1e-5&&Math.abs(back.be-c.be)<1e-6,`${t.casa}: INNO → CasaXPS → INNO is the identity (width, origin)`,`${W} ${e1(rel(back[W],c[W]))} · be ${e1(Math.abs(back.be-c.be))}`);
    // CasaXPS 판독(pos, fwhm) → INNO: 맞춘 INNO 성분의 폭·원점을 되찾는가
    const fr=C.casaToComp(t.casa,{be:k.pos,area:k.area,fwhm:k.fwhm,strict:true}), fc=fr.comp;
    const tolW=FW_TOL[kind], posTol=kind==="ST"?null:2e-4;
    ok(rel(fc[W],c[W])<tolW&&(posTol==null||Math.abs(fc.be-c.be)<posTol)&&fr.exact===(kind==="smooth"&&!t.file.startsWith("K_DS")),`CasaXPS ${t.casa} (pos ${k.pos}, fwhm ${k.fwhm}) → INNO ${W} ≤ ${tolW}${posTol?", origin ≤ 2e-4 eV":""}; exact = ${kind==="smooth"&&!t.file.startsWith("K_DS")}`,
      `${W} ${fc[W].toFixed(6)} vs fitted ${c[W].toFixed(6)} (${e1(rel(fc[W],c[W]))}), be ${fc.be.toFixed(4)} vs ${c.be.toFixed(4)}`); }
  { const r=C.casaToComp("DS(0.1,100)",{be:285.0523,fwhm:0.907537}); ok(Math.abs(r.comp.be-284.9996)<2e-4&&r.posShift>0.05&&/not pinned/.test(r.note),"casaToComp('DS(0.1,100)', pos 285.0523, fwhm 0.907537) → INNO origin 284.9996 ± 2e-4; area normalisation flagged 'not pinned'",`be ${r.comp.be.toFixed(5)}, shift ${r.posShift.toFixed(5)}`); }
  { const a=C.casaToComp("ST(1.7)LA(1,1,200)",{be:285,fwhm:1}), b=C.casaToComp("ST(1.7)LA(1,1,200)",{be:285,fwhm:1,innoUnits:true}), tf=C.trueFWHM("LA",C.shapeParams(b.comp,1)).fwhm;
    ok(Math.abs(tf-1)<1e-6&&a.comp.be===285&&a.posShift===null&&!a.exact&&/0\.0984/.test(a.note),"ST: width from the ST-modified shape (true FWHM = given width), position offset reported 'not pinned'",`true FWHM ${tf.toFixed(6)}`); }
  { const g=C.casaToComp("GL(30)T(2)",{be:285,fwhm:1.18051}).comp; ok(rel(g.F,1.2)<1e-3,"GL(30)T(2), CasaXPS fwhm 1.18051 → F 1.2 (≤ 1e-3; old rule F = fwhm was −1.6 %)",`F ${g.F.toFixed(5)}`); }

  console.log("== 3 casaModelToINNO: DS fwhm ratio uses wL; position bounds/offsets on CasaXPS positions");
  { const hv=1486.6, ke=b=>hv-b;
    const cm={shift:0,regions:[{name:"C 1s",bg:"Linear",keStart:ke(297.13),keEnd:ke(280.03),aw:2,cross:[]}],comps:[
      {name:"A",shape:"DS(0.1,100)",area:{v:20000,lo:0,hi:1e9,link:-1,k:1},fwhm:{v:0.9,lo:0.5,hi:1.5,link:-1,k:1},pos:{v:ke(285.05),lo:ke(285.5),hi:ke(284.5),link:-1,k:0}},
      {name:"B",shape:"DS(0.1,100)",area:{v:5000,lo:0,hi:1e9,link:-1,k:1},fwhm:{v:0.9,lo:0.5,hi:1.5,link:0,k:1},pos:{v:ke(286.05),lo:ke(286.5),hi:ke(285.5),link:-1,k:0}},
      {name:"G",shape:"GL(30)",area:{v:5000,lo:0,hi:1e9,link:-1,k:1},fwhm:{v:1.2,lo:0.5,hi:2,link:-1,k:1},pos:{v:ke(287),lo:ke(287.5),hi:ke(286.5),link:0,k:1.95}}]};
    const M=C.casaModelToINNO(cm,hv,{step:0.05}); const A=M.comps[0], q=M.cons[A.id], qB=M.cons[M.comps[1].id], sh=A.wL*Math.abs(1/Math.tan(Math.PI/1.9));
    const kA=A.wL/0.9;
    ok(q.wL&&Math.abs(q.wL.lo-0.5*kA)<1e-9&&Math.abs(q.wL.hi-1.5*kA)<1e-9&&qB.wL&&qB.wL.mode==="link"&&!q.F,"DS: MFWHM bounds converted with wL/fwhm (not 1), FWHM link → wL link",`wL ${A.wL.toFixed(5)}, bounds ${q.wL.lo.toFixed(4)}–${q.wL.hi.toFixed(4)}`);
    ok(Math.abs(A.be-(285.05-sh))<1e-9&&Math.abs(q.be.lo-(284.5-sh))<1e-9&&Math.abs(q.be.hi-(285.5-sh))<1e-9,"DS: INNO origin and position bounds = CasaXPS values − wL·|cot(π/(2−α))|",`shift ${sh.toFixed(5)} eV`);
    ok(M.notes.some(n=>/G: position link not reproduced.*differs/.test(n)),"position link to a component with a different CasaXPS–INNO shift: noted",M.notes.find(n=>/^G:/.test(n))); }

  console.log(`\n${pass} 통과 · ${fail} 실패`); if(fail){ console.log(fails.map(s=>" - "+s).join("\n")); process.exit(1); }
})().catch(e=>{ console.error(e); process.exit(2); });
