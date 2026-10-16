// 내장 결정상 팩(src/diff_pack.js) 검사: 기록 검증·왕복, CIF 원문 대조, gemmi·pymatgen 기준값(test/ref_pack.json, ref/gen_pack.py), 화학 조성·밀도
const fs=require("fs"), path=require("path"), vm=require("vm");
require("../src/diff_sg.js"); require("../src/diff_scat.js"); const DC=require("../src/diff_core.js");
const PACK=require("../src/diff_pack.js");
const SCAT=globalThis.DIFF_SCAT;
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } if(!c||process.env.V) console.log((c?"  ✓ ":"  ✗ ")+m+(c||!d?"":" — "+d)); };
const REF=JSON.parse(fs.readFileSync(path.join(__dirname,"ref_pack.json"),"utf8"));
const PACKDIR=path.join(__dirname,"..","ref","pack");
const NA=6.02214076e23;

console.log("== 팩 전체");
{ const code=fs.readFileSync(path.join(__dirname,"..","src","diff_pack.js"),"utf8");
  const c1={self:{}}; vm.runInNewContext(code,c1); const c2={}; vm.runInNewContext(code,c2);
  ok(c1.self.DIFF_PACK&&c1.self.DIFF_PACK.phases.length===PACK.phases.length,"브라우저(self) 전역 DIFF_PACK");
  ok(c2.DIFF_PACK&&c2.DIFF_PACK.phases.length===PACK.phases.length,"self 없는 환경(globalThis) 전역 DIFF_PACK");
  ok(PACK.version==="2026-10-17"&&PACK.phases.length===Object.keys(REF.phases).length,`version ${PACK.version}, ${PACK.phases.length} 상 = 기준값 ${Object.keys(REF.phases).length}`);
  const ids=PACK.phases.map(r=>r.id); ok(new Set(ids).size===ids.length&&ids.every(i=>/^pack:[A-Za-z0-9-]+$/.test(i)),"id 유일·pack: 접두"); }

const sameNum=(a,b,t)=>Math.abs(a-b)<=(t||1e-12)*Math.max(1,Math.abs(a),Math.abs(b));
let nRefl=0, nF=0;
for(const rec of PACK.phases){
  const short=rec.id.slice(5), R=REF.phases[short]; console.log(`== ${rec.id}  ${rec.name}  (COD ${rec.source.codId}, ${rec.sg.symbol})`);
  if(!R){ ok(false,`${short}: 기준값 없음`); continue; }
  // 1. 기록 검증 · 왕복
  const v=DC.validateRecord(rec); ok(v.ok&&!v.errors.length,"validateRecord 오류 없음",v.errors.join("; "));
  if(v.warnings.length) console.log("    경고: "+v.warnings.join(" / "));
  ok(rec.verified===true&&rec.source.kind==="COD"&&/^\d{7}$/.test(rec.source.codId)&&rec.source.citation,"출처(COD id·인용)·verified");
  const ph=DC.phaseFromRecord(rec); const r2=DC.recordFromPhase(ph,{id:rec.id}); const ph2=DC.phaseFromRecord(r2);
  ok(["a","b","c","al","be","ga"].every(k=>r2.cell[k]===rec.cell[k]&&ph2.cell[k]===ph.cell[k])&&JSON.stringify(r2.symops)===JSON.stringify(rec.symops)&&ph2.ops.length===ph.ops.length
     &&r2.atoms.length===rec.atoms.length&&r2.atoms.every((a,i)=>["label","el","x","y","z","occ","B"].every(k=>a[k]===rec.atoms[i][k])),"phaseFromRecord ↔ recordFromPhase 왕복");
  ok(rec.symops.every(s=>DC.formatSymop(DC.parseSymop(s))===s),`대칭 연산 ${rec.symops.length} 개 문자열 = formatSymop 꼴`);
  // 2. 단위 (Å → nm) 와 CIF 원문
  const A=R.cellA; ok(sameNum(rec.cell.a*10,A[0],1e-12)&&sameNum(rec.cell.b*10,A[1],1e-12)&&sameNum(rec.cell.c*10,A[2],1e-12)&&rec.cell.a>0.15&&rec.cell.c<2.5,`격자 nm = CIF Å ÷ 10 (a ${rec.cell.a} nm = ${A[0]} Å)`);
  const cp=DC.parseCIF(fs.readFileSync(path.join(PACKDIR,R.cif),"utf8"))[0];
  const so=new Set(cp.ops.map(DC.formatSymop));
  ok(["a","b","c","al","be","ga"].every(k=>sameNum(cp.cell[k],ph.cell[k],1e-12))&&cp.sg.number===rec.sg.number&&so.size===rec.symops.length&&rec.symops.every(s=>so.has(s)),"DC.parseCIF(원문) 격자·공간군·연산 = 기록");
  const bSkip=/not used \(B = 0\)/.test(rec.note||"");
  ok(cp.atoms.length===rec.atoms.length&&cp.atoms.every((a,i)=>{ const b=rec.atoms[i]; return a.el===b.el&&a.label===b.label&&sameNum(a.x,b.x)&&sameNum(a.y,b.y)&&sameNum(a.z,b.z)&&sameNum(a.occ,b.occ)&&(bSkip?b.B===0:Math.abs(a.B-b.B)<1e-5); }),
     "DC.parseCIF(원문) 원자 자리 = 기록"+(bSkip?" (B 는 note 대로 0)":""));
  // 3. 기준값: 반사 목록·d·소멸칙·|F|²
  const list=DC.reflections(ph,REF.gmax,{allowed:false}); const key=h=>h.join(","); const mine=new Map(list.map(r=>[key(r.h),r]));
  const rows=R.rows; let Fmax=0; for(const r of rows) Fmax=Math.max(Fmax,r[5]);
  let miss=0, bd=0, ba=0, bf=0, bz=0, bzMax=0, worst=0;
  for(const [h,k,l,d,abs,F2] of rows){ const hv=[h,k,l]; const g=1/d; const m=mine.get(key(hv));
    if(!m){ if(Math.abs(g-REF.gmax)>1e-9) miss++; continue; }
    if(Math.abs(m.d-d)>1e-9*d) bd++;
    const a=DC.absent(ph.ops,hv); if(a!==abs) ba++;
    const I=DC.structureFactor(ph,hv).I;
    if(F2>=1e-6*Fmax){ const rel=Math.abs(I-F2)/F2; if(rel>2e-4){ bf++; worst=Math.max(worst,rel); } nF++; } else if(Math.abs(I-F2)>1e-6*Fmax) bf++;
    if(abs){ bzMax=Math.max(bzMax,I/Fmax); if(I>1e-6*Fmax) bz++; } }
  nRefl+=rows.length;
  ok(!miss&&list.length===rows.length,`g ≤ ${REF.gmax} 1/nm 반사 목록 = gemmi (${rows.length} hkl)`,`빠짐 ${miss}, DC ${list.length} vs ${rows.length}`);
  ok(!bd&&!ba,`d·소멸칙 = gemmi (허용 ${rows.filter(r=>!r[4]).length})`,`d ${bd} 소멸 ${ba}`);
  ok(!bf,`전자 |F|² = pymatgen 자리 × gemmi c4322 (상대 2e-4, 최대의 1e-6 이상)`,`다른 ${bf} (최악 ${worst.toExponential(2)})`);
  // 4 자리로 반올림된 특수 좌표(0.3333 등)는 소멸 반사에 아주 작은 잔여 세기를 남긴다 → 최대의 1e-6 아래면 0 으로 본다
  ok(!bz,`소멸 반사의 |F|² < 최대의 1e-6 (구조와 연산이 맞음, 최대 잔여 ${bzMax.toExponential(1)})`,`${bz} 개`);
  // 4. 화학 조성 · 밀도
  const uc=DC.expandAtoms(ph); const comp={}; for(const a of uc) comp[a.el]=(comp[a.el]||0)+a.occ;
  const fz=DC.parseFormula(R.formula_sum), fr=DC.parseFormula(rec.formula); const els=new Set([...Object.keys(comp),...Object.keys(fz)]);
  const cOK=[...els].every(e=>Math.abs((comp[e]||0)/R.Z-(fz[e]||0))<=1e-3*Math.max(1,fz[e]||0)&&Math.abs((fr[e]||0)-(fz[e]||0))<=1e-3*Math.max(1,fz[e]||0));
  const pOK=[...els].every(e=>Math.abs((comp[e]||0)-(R.comp_cell[e]||0))<=1e-6*Math.max(1,comp[e]||0));
  ok(cOK&&pOK,`단위세포 조성 ÷ Z ${R.Z} = CIF 화학식 '${R.formula_sum}' = 기록 '${rec.formula}' (pymatgen 자리 조성과 같음)`,JSON.stringify(comp));
  let M=0; for(const a of uc) M+=a.occ*SCAT[a.el][1]; const rho=M/(ph.cell.V*1e-21*NA);
  if(R.density_cif!=null){ const dv=(rho-R.density_cif)/R.density_cif; ok(Math.abs(dv)<0.03,`밀도 ${rho.toFixed(3)} g/cm³ vs CIF ${R.density_cif} (${(dv*100).toFixed(2)} %)`); }
  else { const dv=(rho-R.density_pmg)/R.density_pmg; ok(Math.abs(dv)<0.01,`밀도 ${rho.toFixed(3)} g/cm³ (CIF 값 없음 — pymatgen ${R.density_pmg.toFixed(3)} 과 ${(dv*100).toFixed(2)} %)`); }
}
console.log(`\n반사 ${nRefl} 개 (|F|² 상대 비교 ${nF} 개)`);
console.log(`${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
