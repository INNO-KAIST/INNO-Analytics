/* CasaXPS 2.3.27 Round-3–6 실측 대조 (SPEC_casa_v3.md · SPEC_casa_v4.md · test/ref/casa_readings_v3.md · ref_casa_v3.py → ref_casa_v3.json · ref_casa_v4.py → ref_casa_v4.json)
   node test/test_casa_v3.js   (modules/xps 에서, 약 1 분)
   1) INNO CasaXPS 규약 피팅(weight "casa", Linear conv casa Av.W 1, 영역 280.03/297.13) 으로 K_*.vms 및 V2 행을 맞춤 →
      (a) 독립 scipy 기준(ref_casa_v3.json) 과 넓이·위치·폭·RMS, (b) CasaXPS 판독과 넓이·fwhm·위치·RMS·DoF·asym
   2) 표기 왕복: compToCasa(맞춘 성분) 표기 = CasaXPS 문자열; casaToComp(문자열, CasaXPS pos·fwhm) → 맞춘 INNO 성분 (폭·위치)
   3) casaModelToINNO: DS fwhm 연결 비 = wL, 위치 한계는 CasaXPS 위치 − posShift, DS 넓이·한계·연결 = CasaXPS 넓이 / dsCasaAreaFactor
   4) Round 4–6 (SPEC_casa_v4): DS 넓이(±600·wL)·위치, 폭 숫자 규칙(격자 F/40 포물선, 2·min(l,h)), ST 모양(계단 D/2 늦춤)·위치(가장 큰 마디 k·D) */
const path=require("path"), fs=require("fs");
const D=p=>path.join(__dirname,p);
const C=require("../src/xps_core.js"), IO=require("../src/xps_io.js");
const REF=JSON.parse(fs.readFileSync(D("ref/ref_casa_v3.json"),"utf8")), REF4=JSON.parse(fs.readFileSync(D("ref/ref_casa_v4.json"),"utf8"));
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
const rel=(a,b)=>Math.abs(a-b)/Math.abs(b), e1=x=>x.toExponential(1);
const I0=19, I1=361;
// 판정 기준 (SPEC_casa_v4 §3–§4, casa_v4_patch §7): CasaXPS fwhm·asym 은 하나의 규칙으로 정함 — ST 의 바탕에 ⊗G(n > 0)가 있으면(STG) 근사
const stBaseN=s=>{ const m=s.match(/^ST\([^)]*\)(TLA|LA|SGL|GL)\(([^)]*)\)/); if(!m) return 0; const a=m[2].split(",").map(Number); return +(m[1]==="LA"||m[1]==="TLA"?a[2]:a[1])||0; };
const KIND=s=>/^ST/.test(s)?(stBaseN(s)>0?"STG":"ST"):/T\(/.test(s)?"tail":/^LA\(([\d.]+),([\d.]+),0\)/.test(s)&&(()=>{ const m=s.match(/^LA\(([\d.]+),([\d.]+)/); return m[1]!==m[2]; })()?"laKink":"smooth";
const FW_TOL={smooth:6e-4,laKink:6e-4,tail:6e-4,ST:6e-4,STG:6e-3}, ASYM_TOL={smooth:1.2e-3,laKink:1e-3,tail:1e-3,ST:1.2e-3,STG:1.2e-2};
const isST=k=>k==="ST"||k==="STG";
const ST6=Object.fromEntries(REF4.st6.map(t=>[t.file,t]));   // ST 기준(ref_casa_v4 st6): CasaXPS 모양 = 계단 D/2 늦춤, D = K·F/40

(async()=>{
  const rd=async f=>(await IO.readFile(f+".vms",new Uint8Array(fs.readFileSync(D("data/casaval/"+f+".vms"))))).spectra[0];
  console.log("== 1 CasaXPS 규약 피팅 — 기준 구현(scipy) 과 CasaXPS 2.3.27 판독");
  const fitted={};
  for(const t of REF.fits){ const k=t.read, s=await rd(t.file);
    // 시작값: 판독 fwhm 을 조금 흔든 값에서 (casaToComp 로 F 환산) — 맞춤이 판독값을 그대로 받지 않도록
    const r0=C.casaToComp(t.casa,{be:285.05,area:12000,fwhm:k.fwhm*1.06,strict:true}); const comp=r0.comp; comp.id="a";
    const cons={a:{}}; C.paramNames(comp).forEach(q=>{ if(["be","area","F"].includes(q)||(comp.shape==="DS"&&q==="wL")||C.derivedParam(comp,q)) return; cons.a[q]={mode:"fixed"}; });
    const P=C.makeProblem({be:s.be,y:s.y,c:1},{lo:s.be[I0],hi:s.be[I1],bg:{type:"linear",conv:"casa",aw:1},comps:[comp],cons,weight:"casa"});
    const r=C.fitOnce(P,{}); const st=C.fitStats(P,r); const c=P.comps[0], win={emin:P.E[0],emax:P.E[P.E.length-1]}, cz=C.compToCasa(c,win), W=c.shape==="DS"?c.wL:c.F; fitted[t.file]={c,cz,st};
    // (a) 기준 구현: 같은 정의·같은 CasaXPS 피팅 규약 → 숫자 일치 (DS 넓이 = INNO 창 넓이 ≠ 기준의 진폭 매개변수 → 넓이 제외)
    //     ST: casaToComp 가 stCasaK 를 붙여 CasaXPS 모양(계단 D/2 늦춤)으로 맞춤 → ref_casa_v4 st6 와 대조 (LA 바탕은 창 밖 꼬리 때문에 넓이 2e-3)
    const t6=isST(KIND(t.casa))?ST6[t.file]:null, tr=t6||t, aTol=t6&&/LA\(/.test(t.casa)?2e-3:1e-4;
    const aOK=tr.area==null||rel(c.area,tr.area)<aTol, rOK=tr.rms<1?st.casa.rms<0.05:rel(st.casa.rms,tr.rms)<1e-4;
    ok(aOK&&Math.abs(c.be-tr.pos)<1e-4&&rel(W,tr.F)<1e-4&&rOK&&st.casa.dof===t.dof_rule,`${t.file} ${t.casa}: INNO CasaXPS-mode fit = scipy reference${t6?" (ref_casa_v4 st6, lagged step)":""} (area${aTol>1e-4?" ≤ "+aTol:""}, pos, width, RMS ≤ 1e-4; DoF rule)`,
      `${tr.area==null?"area — (DS)":"area "+e1(rel(c.area,tr.area))} · pos ${e1(Math.abs(c.be-tr.pos))} eV · ${c.shape==="DS"?"wL":"F"} ${e1(rel(W,tr.F))} · RMS ${st.casa.rms.toFixed(3)} / ${tr.rms.toFixed(3)} · DoF ${st.casa.dof}/${t.dof_rule}`);
    // (b) CasaXPS 판독
    const kind=KIND(t.casa), sh=cz.posShift, cpos=c.be+(sh||0);
    const vs=`area ${c.area.toFixed(1)} vs ${k.area}, fwhm ${cz.fwhm!=null?cz.fwhm.toFixed(6):"—"} vs ${k.fwhm} (${cz.fwhm!=null?e1(rel(cz.fwhm,k.fwhm)):"—"}), pos ${cpos.toFixed(4)}${sh?` (INNO ${c.be.toFixed(4)} + ${sh.toFixed(4)})`:""} vs ${k.pos}, RMS ${st.casa.rms.toFixed(3)} vs ${k.rms}, DoF ${st.casa.dof} vs ${k.dof}, asym ${cz.asym!=null?cz.asym.toFixed(4):"—"} vs ${k.asym}`;
    ok(st.casa.dof===k.dof,`${t.file}: CasaXPS DoF = #{channels with y ≠ background} − 3`,`${st.casa.dof} vs ${k.dof}`);
    const isDS=t.file.startsWith("K_DS"), areaOK=isDS?cz.area!=null&&rel(cz.area,k.area)<5e-4:rel(c.area,k.area)<(kind==="STG"?2e-3:1e-3);
    const posOK=Math.abs(cpos-k.pos)<2e-4, rmsOK=rel(st.casa.rms,k.rms)<0.035||Math.abs(st.casa.rms-k.rms)<20;   // CasaXPS RMS floor on exact data: up to 22.9 (LA(2,3,0)), SPEC_casa_v3 §2
    const fwOK=cz.fwhm!=null&&rel(cz.fwhm,k.fwhm)<FW_TOL[kind], asOK=cz.asym!=null&&Math.abs(cz.asym-k.asym)<ASYM_TOL[kind];
    const exOK=kind==="STG"?cz.fwhmExact===false:cz.fwhmExact===true;
    ok(areaOK&&posOK&&rmsOK&&fwOK&&asOK&&exOK&&cz.notation===t.casa,`${t.file} vs CasaXPS ${t.casa}: ${isDS?"area (× dsCasaAreaFactor) ≤ 5e-4":"area ≤ "+(kind==="STG"?"2e-3":"1e-3")}, pos ≤ 2e-4 eV, RMS ≤ 3.5 % or 20 CPS, fwhm ≤ ${FW_TOL[kind]}, asym ≤ ${ASYM_TOL[kind]} (${kind==="STG"?"⊗G base — marked approximate":"pinned"})`,
      vs+(isDS&&cz.area!=null?`, CasaXPS area ${cz.area.toFixed(1)} (×${cz.areaFactor.toFixed(5)})`:""));
  }
  console.log("== 1b asym index only (V2_GL30T2 data fitted as GL(30)T(1), GL(0)T(2) in CasaXPS): grid-parabola rule (SPEC_casa_v4 §3)");
  for(const t of REF.shape_only){ const c=C.casaToComp(t.shape,{fwhm:1,strict:true}).comp; c.F=1; const cz=C.compToCasa(c);
    ok(Math.abs(cz.asym-t.asym_read)<6e-4&&cz.fwhmExact===true,`${t.shape}: asym ${cz.asym.toFixed(4)} vs CasaXPS ${t.asym_read} (≤ 6e-4, pinned)`,`fwhm/F ${cz.fwhm.toFixed(5)} (v3 δ-implied ${t.w_with_delta.toFixed(5)})`); }

  console.log("== 2 표기 왕복: INNO 성분 → CasaXPS (표기·fwhm·pos) → INNO 성분");
  for(const t of REF.fits){ const {c,cz}=fitted[t.file], k=t.read, kind=KIND(t.casa);
    const back=C.casaToComp(cz.notation,{be:cz.posShift!=null?c.be+cz.posShift:c.be,area:c.area,fwhm:cz.fwhm,strict:true}).comp, W=c.shape==="DS"?"wL":"F";
    const same=["shape","mod","tail","casaN","stMu","stN","alpha","beta","mu","m","tk","tj"].every(q=>(c[q]==null&&back[q]==null)||c[q]===back[q]||Math.abs((+c[q]||0)-(+back[q]||0))<1e-9);
    ok(same&&rel(back[W],c[W])<1e-5&&Math.abs(back.be-c.be)<1e-6,`${t.casa}: INNO → CasaXPS → INNO is the identity (width, origin)`,`${W} ${e1(rel(back[W],c[W]))} · be ${e1(Math.abs(back.be-c.be))}`);
    // CasaXPS 판독(pos, fwhm) → INNO: 맞춘 INNO 성분의 폭·원점을 되찾는가
    const fr=C.casaToComp(t.casa,{be:k.pos,area:k.area,fwhm:k.fwhm,strict:true}), fc=fr.comp;
    // casaToComp 의 exact: 매끈·꺾인 LA·T(k) 모두 정함(SPEC_casa_v4 §3; compToCasa 와 같게) → 판독 그대로 맞춤; DS 는 창(넓이 환산)이 없으면 false; ST ⊗G 바탕은 근사
    const tolW=kind==="STG"?6e-3:FW_TOL[kind], posTol=kind==="STG"?1e-3:2e-4, exp=kind==="STG"||t.file.startsWith("K_DS")?false:true;
    ok(rel(fc[W],c[W])<tolW&&Math.abs(fc.be-c.be)<posTol&&fr.exact===exp,`CasaXPS ${t.casa} (pos ${k.pos}, fwhm ${k.fwhm}) → INNO ${W} ≤ ${tolW}, origin ≤ ${posTol} eV; exact = ${exp}`,
      `${W} ${fc[W].toFixed(6)} vs fitted ${c[W].toFixed(6)} (${e1(rel(fc[W],c[W]))}), be ${fc.be.toFixed(4)} vs ${c.be.toFixed(4)} (${e1(fc.be-c.be)})`); }
  { const r=C.casaToComp("DS(0.1,100)",{be:285.0523,fwhm:0.907537}), {c}=fitted.K_DS01_100_F1, r2=C.casaToComp("DS(0.1,100)",{be:285.0523,fwhm:0.907537,area:20811.2,emin:280.03,emax:297.13});
    ok(Math.abs(r.comp.be-284.9996)<2e-4&&r.posShift>0.05&&Math.abs(r.posShift-r.comp.wL*Math.abs(1/Math.tan(Math.PI/1.9)))<1e-12&&/pass the fit window/.test(r.note)&&!r.areaFactor,"casaToComp('DS(0.1,100)', pos 285.0523, fwhm 0.907537) → INNO origin 284.9996 ± 2e-4, posShift = wL·|cot(π/1.9)|; without a window the area is kept and the note asks for it",`be ${r.comp.be.toFixed(5)}, shift ${r.posShift.toFixed(5)}`);
    ok(r2.areaFactor>1&&rel(r2.comp.area,c.area)<5e-4&&r2.exact===true,"casaToComp('DS(0.1,100)', area 20811.2, window 280.03–297.13) → INNO window area = CasaXPS area / dsCasaAreaFactor ≈ INNO fitted area (≤ 5e-4)",`${r2.comp.area.toFixed(1)} vs ${c.area.toFixed(1)} (${e1(rel(r2.comp.area,c.area))}), factor ${r2.areaFactor.toFixed(5)}`);
    const q=C.compToCasa(r2.comp,{emin:280.03,emax:297.13}); ok(rel(q.area,20811.2)<1e-9,"DS area round trip: compToCasa(casaToComp(…), window).area = CasaXPS area",q.area.toFixed(3)); }
  { const a=C.casaToComp("ST(1.7)LA(1,1,200)",{be:285,fwhm:1}), b=C.casaToComp("ST(1.7)LA(1,1,200)",{be:285,fwhm:1,innoUnits:true}), tf=C.trueFWHM("LA",C.shapeParams(b.comp,1)).fwhm, cw=C.casaWidthAsym(a.comp);
    ok(Math.abs(tf-1)<1e-6&&b.comp.stCasaK==null&&cw.node!=null&&a.posShift===cw.node&&Math.abs(a.comp.be-(285-cw.node))<1e-12&&Math.abs(a.comp.stCasaK-1.99129)<1e-4&&Math.abs(cw.node/(a.comp.stCasaK*a.comp.F/40)-Math.round(cw.node/(a.comp.stCasaK*a.comp.F/40)))<1e-9,
      "ST: casaToComp sets stCasaK = base FWHM/F (1.99129 for LA(1,1,200)) and posShift = largest node k·D of CasaXPS's ST table; innoUnits: true FWHM = given width, no stCasaK",`K ${a.comp.stCasaK.toFixed(5)}, node ${cw.node.toFixed(5)} eV = ${(cw.node/(a.comp.stCasaK*a.comp.F/40)).toFixed(3)}·D`);
    ok(C.casaToComp("ST(1.7)GL(30)",{be:285,fwhm:1.1}).exact===true&&a.exact===false,"ST: casaToComp exact for a base without ⊗G (ST(1.7)GL(30)), approximate for LA(1,1,200)"); }
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
    ok(M.notes.some(n=>/G: position link not reproduced.*differs/.test(n)),"position link to a component with a different CasaXPS–INNO shift: noted",M.notes.find(n=>/^G:/.test(n)));
    // DS 넓이 (SPEC_casa_v4 §2): CasaXPS 넓이·한계 ÷ dsCasaAreaFactor(INNO 창), DS ↔ 다른 모양 넓이 연결에 인자가 실림
    const cm2=JSON.parse(JSON.stringify(cm)); cm2.comps[2].area={v:6000,lo:0,hi:1e9,link:0,k:0.5}; cm2.comps[1].area={v:5000,lo:100,hi:9e4,link:-1,k:1};
    const M2=C.casaModelToINNO(cm2,hv,{step:0.05}), A2=M2.comps[0], B2=M2.comps[1], fA=C.dsCasaAreaFactor(A2,M2.lo,M2.hi), fB=C.dsCasaAreaFactor(B2,M2.lo,M2.hi), qG=M2.cons[M2.comps[2].id], qB2=M2.cons[B2.id];
    ok(fA>1&&rel(A2.area,20000/fA)<1e-12&&rel(B2.area,5000/fB)<1e-12&&rel(qB2.area.lo,100/fB)<1e-12&&rel(qB2.area.hi,9e4/fB)<1e-12&&qG.area.mode==="link"&&rel(qG.area.k,0.5*fA)<1e-12&&M2.comps[2].area===6000,
      "DS: CasaXPS areas and area bounds → INNO window areas (÷ dsCasaAreaFactor over the INNO window); a GL area linked to a DS carries the factor (k × factor)",`factor A ${fA.toFixed(5)}, B ${fB.toFixed(5)}, link k ${qG.area.k.toFixed(5)}`); }

  console.log("== 4 Round 4–6 (SPEC_casa_v4, ref_casa_v4.json): new files — scipy reference and CasaXPS readings");
  const fitK=async(file,casa,fw,i1)=>{ const s=await rd(file), comp=C.casaToComp(casa,{be:285.05,area:12000,fwhm:fw*1.06,strict:true}).comp; comp.id="a";
    const cons={a:{}}; C.paramNames(comp).forEach(q=>{ if(["be","area","F"].includes(q)||(comp.shape==="DS"&&q==="wL")||C.derivedParam(comp,q)) return; cons.a[q]={mode:"fixed"}; });
    const P=C.makeProblem({be:s.be,y:s.y,c:1},{lo:s.be[I0],hi:s.be[i1||I1],bg:{type:"linear",conv:"casa",aw:1},comps:[comp],cons,weight:"casa"}); const r=C.fitOnce(P,{}), c=P.comps[0];
    return {c,st:C.fitStats(P,r),cz:C.compToCasa(c,{emin:P.E[0],emax:P.E[P.E.length-1]}),W:c.shape==="DS"?c.wL:c.F}; };
  const done=new Set(REF.fits.map(t=>t.file));
  for(const t of REF4.fits.filter(t=>t.kind!=="ST"&&(!done.has(t.file)||t.i1===220))){ const k=t.i1===220?{area:19946.5,fwhm:0.876718,pos:285.0489,rms:2202.23,dof:197,asym:0.1178}:t.read; if(t.i1===219) continue;
    const casa=t.casa.replace(/ region.*/,""), {c,st,cz,W}=await fitK(t.file,casa,k.fwhm,t.i1), nm=t.file+(t.i1?"@"+t.i1:"");
    ok(Math.abs(c.be-t.pos)<1e-5&&rel(W,t.F)<1e-5&&(rel(st.casa.rms,t.rms)<1e-4||Math.abs(st.casa.rms-t.rms)<0.01)&&st.casa.dof===t.dof_rule&&(t.area_pred==null||rel(cz.area,t.area_pred)<2e-4),`${nm} ${casa}: INNO fit = scipy reference (pos, width, RMS, DoF${t.area_pred!=null?"; CasaXPS area via factor":""})`,
      `pos ${e1(c.be-t.pos)} eV · W ${e1(rel(W,t.F))} · RMS ${e1(rel(st.casa.rms,t.rms))} · DoF ${st.casa.dof}/${t.dof_rule}${t.area_pred!=null?" · area "+e1(rel(cz.area,t.area_pred)):""}`);
    const cpos=c.be+(cz.posShift||0), isDS=c.shape==="DS", aOK=isDS?rel(cz.area,k.area)<5e-4:rel(c.area,k.area)<1e-3;
    ok(aOK&&Math.abs(cpos-k.pos)<2e-4&&rel(cz.fwhm,k.fwhm)<6e-4&&Math.abs(cz.asym-k.asym)<1e-3&&st.casa.dof===k.dof&&cz.fwhmExact===true,`${nm} vs CasaXPS ${casa}: area${isDS?" (× dsCasaAreaFactor) ≤ 5e-4":" ≤ 1e-3"}, pos ≤ 2e-4 eV, fwhm ≤ 6e-4, asym ≤ 1e-3, DoF`,
      `area ${(isDS?cz.area:c.area).toFixed(1)} vs ${k.area}, pos ${cpos.toFixed(4)} vs ${k.pos}, fwhm ${cz.fwhm.toFixed(6)} vs ${k.fwhm} (${e1(rel(cz.fwhm,k.fwhm))}), asym ${cz.asym.toFixed(4)} vs ${k.asym}, DoF ${st.casa.dof}/${k.dof}`);
    if(t.file==="K_LA243_153_0_F1"){ const cw=C.casaWidthAsym(c); ok(cw.side==="high"&&cz.asym<0&&k.asym<0,"LA(2.43,1.53,0) (tail on the low-BE side): CasaXPS fwhm = 2 × the high-BE half width (narrower side), asym = 1 − l/h < 0",`side ${cw.side}, asym ${cz.asym.toFixed(4)} vs ${k.asym}`); } }
  { const t=REF4.ds5, {c,st,cz}=await fitK(t.file,"DS(0.3,0)",t.read.fwhm), k=t.read, cpos=c.be+cz.posShift;
    ok(Math.abs(c.be-t.pos)<1e-5&&rel(c.wL,t.wL)<1e-5&&rel(st.casa.rms,t.rms)<1e-4&&rel(cz.area,k.area)<5e-4&&Math.abs(cpos-k.pos)<2e-4&&rel(cz.fwhm,k.fwhm)<6e-4&&Math.abs(cz.asym-k.asym)<1e-3,"K_DS03_0_F1 DS(0.3,0): = scipy reference; vs CasaXPS area (× factor) ≤ 5e-4, pos ≤ 2e-4 eV, fwhm ≤ 6e-4, asym ≤ 1e-3",
      `area ${cz.area.toFixed(1)} vs ${k.area} (×${cz.areaFactor.toFixed(5)}), pos ${cpos.toFixed(4)} vs ${k.pos}, fwhm ${cz.fwhm.toFixed(6)} vs ${k.fwhm}, asym ${cz.asym.toFixed(4)} vs ${k.asym}`); }
  for(const t of REF4.st6){ const kind=KIND(t.casa), {c,st,cz}=await fitK(t.file,t.casa,t.read.fwhm), k=t.read, cpos=c.be+cz.posShift, lab=/LA\(/.test(t.casa);
    ok(Math.abs(c.be-t.pos)<1e-5&&rel(c.F,t.F)<1e-5&&rel(c.area,t.area)<(lab?2e-3:1e-4)&&rel(st.casa.rms,t.rms)<1e-4&&Math.abs(cz.posShift/(c.stCasaK*c.F/40)-t.k)<1e-6,`${t.file} ${t.casa}: INNO fit (stCasaK, lagged step) = scipy reference st6; CasaXPS pos = origin + ${t.k}·D`,
      `pos ${e1(c.be-t.pos)} eV · F ${e1(rel(c.F,t.F))} · area ${e1(rel(c.area,t.area))} · RMS ${e1(rel(st.casa.rms,t.rms))} · K ${c.stCasaK.toFixed(5)}`);
    const aT=lab?2e-3:2e-4, fT=kind==="STG"?6e-3:5e-4, sT=kind==="STG"?1.2e-2:1.2e-3;
    ok(rel(c.area,k.area)<aT&&rel(cz.fwhm,k.fwhm)<fT&&Math.abs(cz.asym-k.asym)<sT&&Math.abs(cpos-k.pos)<2e-4&&st.casa.dof===k.dof&&cz.exact===(kind!=="STG"),`${t.file} vs CasaXPS ${t.casa}: area ≤ ${aT}, fwhm ≤ ${fT}, asym ≤ ${sT}, pos ≤ 2e-4 eV, DoF${kind==="STG"?" (⊗G base: approximate)":""}`,
      `area ${c.area.toFixed(1)} vs ${k.area} (${e1(rel(c.area,k.area))}), fwhm ${cz.fwhm.toFixed(6)} vs ${k.fwhm} (${e1(rel(cz.fwhm,k.fwhm))}), asym ${cz.asym.toFixed(4)} vs ${k.asym}, pos ${cpos.toFixed(4)} vs ${k.pos} (${e1(cpos-k.pos)}), DoF ${st.casa.dof}/${k.dof}`); }

  console.log(`\n${pass} 통과 · ${fail} 실패`); if(fail){ console.log(fails.map(s=>" - "+s).join("\n")); process.exit(1); }
})().catch(e=>{ console.error(e); process.exit(2); });
