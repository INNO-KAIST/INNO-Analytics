// Diffraction core 대조: gemmi(소멸칙·전자 구조인자, IT C 4.3.2.2)·pymatgen(정대축 반사 d) — 기준 자료 ref/gen_core_ref.py
const fs=require("fs"), path=require("path");
require("../src/diff_sg.js"); require("../src/diff_scat.js"); const DC=require("../src/diff_core.js");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } if(!c||process.env.V) console.log((c?"  ✓ ":"  ✗ ")+m+(c||!d?"":" — "+d)); };
const REF=JSON.parse(fs.readFileSync(path.join(__dirname,"ref_core.json"),"utf8"));
console.log("== CIF → 상, gemmi·pymatgen 대조");
for(const [name,r] of Object.entries(REF)){
  const ph=DC.parseCIF(fs.readFileSync(path.join(__dirname,"cif",r.cif),"utf8"))[0];
  let bd=0, ba=0, bf=0, worstF=0, nI=0; let Fmax=0; for(const row of r.rows) if(!row[4]) Fmax=Math.max(Fmax,row[5]);
  for(const [h,k,l,d,abs,F2] of r.rows){ const hv=[h,k,l]; const dd=DC.dSpacing(ph.cell,hv); if(Math.abs(dd-d)>1e-9*d) bd++;
    const a=DC.absent(ph.ops,hv); if(a!==abs) ba++;
    const F=DC.structureFactor(ph,hv).I; const rel=Math.abs(F-F2)/Math.max(F2,1e-6*Fmax); if(rel>2e-4&&Math.abs(F-F2)>1e-6*Fmax){ bf++; worstF=Math.max(worstF,rel); } nI++; }
  ok(!bd&&!ba&&!bf,`${name} (${r.hm}, ${r.rows.length} hkl): d·소멸칙·|F|² = gemmi`,`d ${bd} 소멸 ${ba} F ${bf} (최악 ${worstF.toExponential(2)})`);
  // 표준 설정 표(번호) 의 소멸칙 = CIF 연산의 소멸칙 (기본 설정과 같은 축일 때만 — 단사정 b 축, 능면체 육방 축)
  const e=DC.sgByNumber(r.sg); const ops=DC.sgOps(e); let bs=0; for(const row of r.rows){ if(DC.absent(ops,row.slice(0,3))!==row[4]) bs++; }
  ok(!bs,`${name}: 공간군 표 #${r.sg} (${e.short}) 소멸칙 = gemmi`,`다른 hkl ${bs}`);
  // pymatgen 정대축 반사 d
  for(const [zb,list] of Object.entries(r.zones)){ let b=0; for(const [h,d] of list){ if(Math.abs(DC.dSpacing(ph.cell,h)-d)>1e-6*d) b++; if(h[0]*+zb.split(",")[0]+h[1]*+zb.split(",")[1]+h[2]*+zb.split(",")[2]!==0) b++; } ok(!b,`${name} [${zb}]: 0 차 Laue 반사 ${list.length} 개 d = pymatgen`,`다른 ${b}`); }
}
console.log("== 공간군 표");
{ const L=DC.sgList(); ok(L.length===530&&new Set(L.map(e=>e.number)).size===230,"530 설정 · 230 공간군");
  const fcc=DC.sgOps(DC.sgByNumber(225)); ok(DC.absent(fcc,[1,0,0])&&!DC.absent(fcc,[1,1,1])&&!DC.absent(fcc,[2,0,0])&&DC.absent(fcc,[2,1,0]),"Fm-3m: 100·210 소멸, 111·200 허용");
  const dia=DC.sgOps(DC.sgByNumber(227)); ok(DC.absent(dia,[2,0,0])&&!DC.absent(dia,[2,2,2])&&!DC.absent(dia,[4,0,0]),"Fd-3m: 200 소멸, 222 공간군상 허용(구조상 0), 400 허용");
  ok(DC.sgBySymbol("P 1 21/c 1").number===14&&DC.sgBySymbol("P21/c").number===14&&DC.sgBySymbol("Fd-3m:2").choice==="2","기호로 찾기 (짧은·완전·설정)"); }
console.log("== 정대축·각도");
{ const fe=DC.phaseFromParams({sgNumber:227,a:0.8396}); const h1=[-2,2,0], h2=[-1,-3,1], h3=[-3,-1,1];
  const a13=DC.angleG(fe.cell,h1,h3), a23=DC.angleG(fe.cell,h2,h3), z=DC.zoneAxis(h1,h2);
  ok(Math.abs(a13-64.8)<0.05&&Math.abs(a23-50.5)<0.05,`REDI 설명서 Fe3O4 예: ∠R1·R3 ${a13.toFixed(2)}° ∠R2·R3 ${a23.toFixed(2)}° (64.8 · 50.5)`);
  ok(JSON.stringify(DC.reduceInt(z.map(x=>Math.abs(x))).sort())===JSON.stringify([1,1,4].sort()),`정대축 ${z} ~ [1 1 4]`); }
console.log("== CIF 읽기 경계");
{ const t=`data_x\n_cell_length_a 4.0782(3)\n_cell_length_b 4.0782\n_cell_length_c 4.0782\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n_symmetry_Int_Tables_number 225\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_U_iso_or_equiv\nAu1 0 0 0 0.0081\n`;
  const p=DC.parseCIF(t)[0]; ok(p.ops.length===192&&p.atoms[0].el==="Au"&&Math.abs(p.atoms[0].B-8*Math.PI**2*0.0081)<1e-9&&DC.expandAtoms(p).length===4,"연산 없는 CIF → 공간군 번호, 표준 불확도 괄호, Uiso→B, 단위세포 4 원자");
  const s=DC.parseSymop("-x+1/2, y+0.5 ,-z"); ok(s.R.join()==="-1,0,0,0,1,0,0,0,-1"&&s.t.join()==="0.5,0.5,0","대칭 연산 문자열"); }
console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
