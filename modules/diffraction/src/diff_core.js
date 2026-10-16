/* ============================================================
   INNO Analytics — Diffraction core (결정학 계산). 브라우저·Worker·Node 공용, DOM 없음.
   · 격자: 직접·역격자(2π 없는 결정학 관례 a*·a = 1), 면간거리, 각도, 정대축(Weiss 띠 법칙).
   · 공간군: DIFF_SG (spglib 데이터베이스 530 설정). 소멸칙 = 연산 (R,t) 에 대해 hR = h 이고 h·t 가 정수가 아니면 F(h) = 0.
   · 구조인자: 운동학, 전자 산란 인자 DIFF_SCAT (IT Vol. C 4.3.2.2, Peng 1996), Debye–Waller exp(−B s²), 점유율.
   · CIF: COD 등 — 격자상수·대칭 연산(xyz 문자열)·원자 자리(Uiso/Biso/occupancy)·공간군 번호·이름.
   단위: 길이 nm, 역공간 1/nm (화면·DB 단위). 산란 인자 표는 Å 기준이라 s 는 Å⁻¹ 로 바꿔 넣는다.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const SG=G.DIFF_SG||(typeof require!=="undefined"?require("./diff_sg.js"):null);
const SCAT=G.DIFF_SCAT||(typeof require!=="undefined"?require("./diff_scat.js"):null);
const DC={};
const D2R=Math.PI/180, R2D=180/Math.PI;
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const norm=a=>Math.sqrt(dot(a,a));
const gcd=(a,b)=>{ a=Math.abs(a); b=Math.abs(b); while(b){ const t=a%b; a=b; b=t; } return a; };
DC.vec={dot,cross,norm};

/* ---------------- 격자 ---------------- */
/* a,b,c (nm), α,β,γ (deg) → 직교 좌표의 직접 기저 A(행 = a,b,c), 역기저 B(행 = a*,b*,c*), 부피, 계량 */
DC.cell=function(a,b,c,al,be,ga){
  if(!(a>0&&b>0&&c>0)) throw new Error("격자상수는 양수여야 합니다");
  const ca=Math.cos(al*D2R), cb=Math.cos(be*D2R), cg=Math.cos(ga*D2R), sg=Math.sin(ga*D2R);
  const v2=1-ca*ca-cb*cb-cg*cg+2*ca*cb*cg; if(!(v2>1e-12)) throw new Error("각 α·β·γ 로 격자를 만들 수 없습니다");
  const V=a*b*c*Math.sqrt(v2);
  const A=[[a,0,0],[b*cg,b*sg,0],[c*cb,c*(ca-cb*cg)/sg,V/(a*b*sg)]];
  const B=[cross(A[1],A[2]).map(x=>x/V),cross(A[2],A[0]).map(x=>x/V),cross(A[0],A[1]).map(x=>x/V)];
  return {a,b,c,al,be,ga,V,A,B};
};
DC.gvec=(cell,h)=>{ const B=cell.B; return [h[0]*B[0][0]+h[1]*B[1][0]+h[2]*B[2][0],h[0]*B[0][1]+h[1]*B[1][1]+h[2]*B[2][1],h[0]*B[0][2]+h[1]*B[1][2]+h[2]*B[2][2]]; };
DC.rvec=(cell,u)=>{ const A=cell.A; return [u[0]*A[0][0]+u[1]*A[1][0]+u[2]*A[2][0],u[0]*A[0][1]+u[1]*A[1][1]+u[2]*A[2][1],u[0]*A[0][2]+u[1]*A[1][2]+u[2]*A[2][2]]; };
DC.gLen=(cell,h)=>norm(DC.gvec(cell,h));                // 1/nm
DC.dSpacing=(cell,h)=>1/DC.gLen(cell,h);                 // nm
DC.angleG=(cell,h1,h2)=>{ const g1=DC.gvec(cell,h1), g2=DC.gvec(cell,h2); const c=dot(g1,g2)/(norm(g1)*norm(g2)); return Math.acos(Math.max(-1,Math.min(1,c)))*R2D; };
/* 정대축: [uvw] = h1 × h2 (Weiss), 가장 작은 정수로 */
DC.zoneAxis=(h1,h2)=>{ let u=cross(h1,h2); const g=gcd(gcd(u[0],u[1]),u[2])||1; u=u.map(x=>x/g); return u; };
DC.reduceInt=v=>{ const g=gcd(gcd(v[0],v[1]),v[2])||1; return v.map(x=>x/g); };

/* ---------------- 공간군 ---------------- */
let SG_CACHE=null;
function sgList(){ if(SG_CACHE) return SG_CACHE; if(!SG) throw new Error("공간군 표(DIFF_SG)가 없습니다");
  SG_CACHE=SG.S.map(e=>({hall:e[0],number:e[1],full:e[2],short:e[3],choice:e[4],hallSymbol:e[5],opsRaw:e[6]})); return SG_CACHE; }
DC.sgList=sgList;
/* 연산: {R:[9], t:[3]} — x' = R x + t (분율 좌표) */
DC.sgOps=function(e){ if(e.ops) return e.ops; e.ops=e.opsRaw.map(o=>{ const r=Math.floor(o/13824), rem=o%13824; return {R:SG.R[r],t:[Math.floor(rem/576)/24,Math.floor((rem%576)/24)/24,(rem%24)/24]}; }); return e.ops; };
/* 기본 설정: 번호마다 spglib 의 첫 Hall 설정 (단사정 b 축 cell choice 1, 원점 choice 1, 능면체 = 육방 축) */
DC.sgByNumber=function(n){ const L=sgList(); const e=L.find(x=>x.number===n); if(!e) throw new Error("공간군 번호 "+n+" 이 없습니다 (1–230)"); return e; };
const normSym=s=>String(s||"").replace(/[\s_]/g,"").replace(/^([A-Z])/,"$1").toLowerCase();
/* 기호로 찾기: 완전·짧은 기호, choice(예 ':2', ':H') 를 붙여도 된다 */
DC.sgBySymbol=function(sym){ const L=sgList(); let s=String(sym||"").trim(); let ch=""; const m=s.match(/^(.*?)\s*:\s*([0-9A-Za-z-]+)$/); if(m){ s=m[1]; ch=m[2]; }
  const k=normSym(s); const hit=L.filter(e=>(normSym(e.full)===k||normSym(e.short)===k||normSym(e.hallSymbol)===k)&&(!ch||e.choice.toLowerCase()===ch.toLowerCase()));
  return hit[0]||null; };
/* 결정계: 번호로 */
DC.crystalSystem=n=>n<=2?"Triclinic":n<=15?"Monoclinic":n<=74?"Orthorhombic":n<=142?"Tetragonal":n<=167?"Trigonal":n<=194?"Hexagonal":"Cubic";
/* 소멸칙: 연산 목록으로 (중심화 평행이동 포함). h 는 정수 3 개 */
DC.absent=function(ops,h){
  for(const op of ops){ const R=op.R;
    // hR (행벡터 × R): (hR)_j = Σ_i h_i R_ij
    const a=h[0]*R[0]+h[1]*R[3]+h[2]*R[6], b=h[0]*R[1]+h[1]*R[4]+h[2]*R[7], c=h[0]*R[2]+h[1]*R[5]+h[2]*R[8];
    if(a===h[0]&&b===h[1]&&c===h[2]){ const ph=h[0]*op.t[0]+h[1]*op.t[1]+h[2]*op.t[2]; if(Math.abs(ph-Math.round(ph))>1e-6) return true; } }
  return false; };
/* Laue 군 궤도 (점군 회전 + 반전): hkl → 같은 족 */
DC.lauePoints=function(ops){ const seen=new Set(), L=[]; for(const op of ops){ const k=op.R.join(","); if(!seen.has(k)){ seen.add(k); L.push(op.R); } const m=op.R.map(x=>-x), k2=m.join(","); if(!seen.has(k2)){ seen.add(k2); L.push(m); } } return L; };
DC.family=function(rots,h){ const out=new Map(); for(const R of rots){ const v=[h[0]*R[0]+h[1]*R[3]+h[2]*R[6],h[0]*R[1]+h[1]*R[4]+h[2]*R[7],h[0]*R[2]+h[1]*R[5]+h[2]*R[8]]; out.set(v.join(","),v); } return [...out.values()]; };
/* 대칭 연산 문자열 'x,y,z' / '-x+1/2, y, -z' → {R,t} */
DC.parseSymop=function(s){ const parts=String(s).replace(/['"]/g,"").toLowerCase().split(","); if(parts.length!==3) throw new Error("대칭 연산 해석 실패: "+s);
  const R=[0,0,0,0,0,0,0,0,0], t=[0,0,0];
  parts.forEach((p,i)=>{ p=p.replace(/\s+/g,""); const re=/([+-]?)(\d*\.?\d*(?:\/\d+)?)?\*?([xyz])?/g; let m; let pos=0;
    const toks=p.match(/[+-]?[^+-]+/g)||[];
    for(let tk of toks){ let sign=1; if(tk[0]==="+") tk=tk.slice(1); else if(tk[0]==="-"){ sign=-1; tk=tk.slice(1); }
      const vm=tk.match(/^([xyz])$/)||tk.match(/^(\d*\.?\d*(?:\/\d+)?)\*?([xyz])$/);
      if(/^[xyz]$/.test(tk)){ R[i*3+"xyz".indexOf(tk)]+=sign; continue; }
      const mm=tk.match(/^(\d+(?:\.\d+)?)(?:\/(\d+(?:\.\d+)?))?\*?([xyz])?$/); if(!mm) throw new Error("대칭 연산 해석 실패: "+s);
      const val=+mm[1]/(mm[2]?+mm[2]:1); if(mm[3]) R[i*3+"xyz".indexOf(mm[3])]+=sign*val; else t[i]+=sign*val; } void re; void m; void pos; });
  return {R,t:t.map(v=>v-Math.floor(v+1e-9))}; };

/* ---------------- 산란 인자·구조인자 ---------------- */
/* f_e(s) [Å], s = sinθ/λ [Å⁻¹] */
DC.fElectron=function(el,s){ const p=SCAT&&SCAT[el]; if(!p) throw new Error("원소 "+el+" 의 전자 산란 인자가 없습니다"); const s2=s*s; let f=0; for(let i=0;i<5;i++) f+=p[2+i]*Math.exp(-p[7+i]*s2); return f; };
DC.elementOf=function(label){ const m=String(label||"").match(/^([A-Z][a-z]?)/); if(!m) return null; if(SCAT[m[1]]) return m[1]; if(SCAT[m[1][0]]) return m[1][0]; return null; };
/* 상(phase): {name, cell, ops, atoms:[{el,x,y,z,occ,B}]} → 단위세포 안 모든 원자 (대칭으로 펼치고 같은 자리는 하나로) */
DC.expandAtoms=function(phase){ if(phase._uc) return phase._uc; const out=[];
  for(const at of phase.atoms||[]){ const seen=[];
    for(const op of phase.ops){ const R=op.R; let x=R[0]*at.x+R[1]*at.y+R[2]*at.z+op.t[0], y=R[3]*at.x+R[4]*at.y+R[5]*at.z+op.t[1], z=R[6]*at.x+R[7]*at.y+R[8]*at.z+op.t[2];
      x-=Math.floor(x); y-=Math.floor(y); z-=Math.floor(z);
      if(seen.some(p=>{ const dx=Math.abs(p[0]-x), dy=Math.abs(p[1]-y), dz=Math.abs(p[2]-z); return Math.min(dx,1-dx)<2e-3&&Math.min(dy,1-dy)<2e-3&&Math.min(dz,1-dz)<2e-3; })) continue;
      seen.push([x,y,z]); out.push({el:at.el,x,y,z,occ:at.occ==null?1:at.occ,B:at.B||0}); } }
  phase._uc=out; return out; };
/* F(hkl) (Å 단위 산란 길이의 합) — 실수부·허수부·|F|² */
DC.structureFactor=function(phase,h){ const uc=DC.expandAtoms(phase); const g=DC.gLen(phase.cell,h); const s=g/2/10; // 1/nm → sinθ/λ in Å⁻¹ (s = g/2, g 는 1/nm → /10 해서 1/Å)
  let re=0, im=0; const cache={};
  for(const a of uc){ const f=cache[a.el]!=null?cache[a.el]:(cache[a.el]=DC.fElectron(a.el,s)); const dw=a.B?Math.exp(-a.B*s*s):1;
    const ph=2*Math.PI*(h[0]*a.x+h[1]*a.y+h[2]*a.z); const w=a.occ*f*dw; re+=w*Math.cos(ph); im+=w*Math.sin(ph); }
  return {re,im,I:re*re+im*im}; };

/* ---------------- 반사 목록 ---------------- */
/* g ≤ gmax (1/nm) 인 모든 hkl (0 제외). opts.allowed: 공간군 소멸칙으로 거른다. opts.F: 원자가 있으면 |F|² 도 */
DC.reflections=function(phase,gmax,opts){ opts=opts||{}; const cell=phase.cell; const out=[];
  const hm=Math.ceil(gmax/ norm(cell.B[0])*1.0001+1e-9)+0, km=Math.ceil(gmax/norm(cell.B[1])), lm=Math.ceil(gmax/norm(cell.B[2]));
  // 경사 격자에서 안전한 범위: |h| ≤ gmax·a (a = |a| 직접), 같은 방식으로 k,l
  const H=Math.ceil(gmax*cell.a), K=Math.ceil(gmax*cell.b), L=Math.ceil(gmax*cell.c); void hm; void km; void lm;
  for(let h=-H;h<=H;h++) for(let k=-K;k<=K;k++) for(let l=-L;l<=L;l++){ if(!h&&!k&&!l) continue; const v=[h,k,l]; const g=DC.gLen(cell,v); if(g>gmax) continue;
    if(opts.allowed!==false&&phase.ops&&DC.absent(phase.ops,v)) continue;
    const r={h:v,g,d:1/g}; if(opts.F&&phase.atoms&&phase.atoms.length){ const F=DC.structureFactor(phase,v); r.I=F.I; } out.push(r); }
  out.sort((a,b)=>a.g-b.g||(b.h[0]-a.h[0])||(b.h[1]-a.h[1])||(b.h[2]-a.h[2])); return out; };

/* 정대축 [uvw] 의 0 차 Laue 영역 회절도: 2D 좌표(e1 = 첫 기준 방향) */
DC.zonePattern=function(phase,uvw,gmax,opts){ opts=opts||{}; const cell=phase.cell; const n=DC.rvec(cell,uvw); const nn=norm(n); const z=n.map(x=>x/nn);
  const refl=DC.reflections(phase,gmax,{allowed:opts.allowed,F:opts.F}).filter(r=>r.h[0]*uvw[0]+r.h[1]*uvw[1]+r.h[2]*uvw[2]===0);
  let e1=null; if(opts.ref){ const g=DC.gvec(cell,opts.ref); e1=g.map(x=>x/norm(g)); } else if(refl.length){ const g=DC.gvec(cell,refl[0].h); e1=g.map(x=>x/norm(g)); }
  if(!e1) return {zone:uvw,spots:[]}; const e2=cross(z,e1);
  const spots=refl.map(r=>{ const g=DC.gvec(cell,r.h); return Object.assign({},r,{x:dot(g,e1),y:dot(g,e2)}); });
  return {zone:uvw,axis:z,e1,e2,spots}; };

/* 다결정 링: d 가 같은 반사를 한 링으로 (다중도 m, Σ|F|²) */
DC.ringList=function(phase,gmax,opts){ opts=opts||{}; const refl=DC.reflections(phase,gmax,{allowed:true,F:!!(phase.atoms&&phase.atoms.length)}); const rings=[]; const tol=opts.tol||1e-6;
  for(const r of refl){ const last=rings[rings.length-1]; if(last&&Math.abs(last.g-r.g)<=tol*r.g){ last.m++; last.hkls.push(r.h); if(r.I!=null) last.I+=r.I; } else rings.push({g:r.g,d:r.d,m:1,hkls:[r.h],I:r.I!=null?r.I:null}); }
  return rings; };

/* ---------------- CIF ---------------- */
/* 숫자 '1.234(5)' → 1.234 */
const cifNum=v=>{ if(v==null) return NaN; const s=String(v).replace(/\(.*?\)/g,"").trim(); if(s==="."||s==="?") return NaN; return parseFloat(s); };
function cifTokens(text){ const toks=[]; const lines=String(text).replace(/\r/g,"").split("\n");
  for(let i=0;i<lines.length;i++){ let ln=lines[i];
    if(ln.startsWith(";")){ let s=ln.slice(1); i++; while(i<lines.length&&!lines[i].startsWith(";")){ s+="\n"+lines[i]; i++; } toks.push({v:s.trim(),q:true}); continue; }
    let j=0; while(j<ln.length){ const ch=ln[j]; if(ch===" "||ch==="\t"){ j++; continue; } if(ch==="#") break;
      if(ch==="'"||ch==='"'){ let k=j+1; while(k<ln.length&&!(ln[k]===ch&&(k+1>=ln.length||/\s/.test(ln[k+1])))) k++; toks.push({v:ln.slice(j+1,k),q:true}); j=k+1; continue; }
      let k=j; while(k<ln.length&&!/\s/.test(ln[k])) k++; toks.push({v:ln.slice(j,k),q:false}); j=k; } }
  return toks; }
/* CIF 글 → 블록들 [{name, items:{key:value}, loops:[{keys:[...], rows:[[...]]}]}] */
DC.cifBlocks=function(text){ const T=cifTokens(text); const blocks=[]; let cur=null; let i=0;
  while(i<T.length){ const t=T[i];
    if(!t.q&&/^data_/i.test(t.v)){ cur={name:t.v.slice(5),items:{},loops:[]}; blocks.push(cur); i++; continue; }
    if(!cur){ i++; continue; }
    if(!t.q&&/^loop_$/i.test(t.v)){ i++; const keys=[]; while(i<T.length&&!T[i].q&&T[i].v[0]==="_"){ keys.push(T[i].v.toLowerCase()); i++; }
      const vals=[]; while(i<T.length&&!(!T[i].q&&(T[i].v[0]==="_"||/^loop_$/i.test(T[i].v)||/^data_/i.test(T[i].v)))){ vals.push(T[i].v); i++; }
      const rows=[]; for(let k=0;k+keys.length<=vals.length;k+=keys.length) rows.push(vals.slice(k,k+keys.length)); cur.loops.push({keys,rows}); continue; }
    if(!t.q&&t.v[0]==="_"){ const key=t.v.toLowerCase(); const v=T[i+1]; if(v&&!(!v.q&&(v.v[0]==="_"||/^loop_$/i.test(v.v)||/^data_/i.test(v.v)))){ cur.items[key]=v.v; i+=2; } else { cur.items[key]=""; i++; } continue; }
    i++; }
  return blocks; };
const loopOf=(b,key)=>b.loops.find(l=>l.keys.includes(key));
/* 블록 → 상 (DB 항목 모양). 원자가 없으면 공간군만 */
DC.phaseFromCIF=function(b){ const it=b.items; const g=k=>it[k];
  const a=cifNum(g("_cell_length_a")), bb=cifNum(g("_cell_length_b")), c=cifNum(g("_cell_length_c")), al=cifNum(g("_cell_angle_alpha")), be=cifNum(g("_cell_angle_beta")), ga=cifNum(g("_cell_angle_gamma"));
  if(!(a>0&&bb>0&&c>0)) throw new Error(`CIF '${b.name}': 격자상수가 없습니다`);
  const cell=DC.cell(a/10,bb/10,c/10,al||90,be||90,ga||90);   // Å → nm
  let ops=null; const L=loopOf(b,"_space_group_symop_operation_xyz")||loopOf(b,"_symmetry_equiv_pos_as_xyz");
  if(L){ const k=L.keys.indexOf("_space_group_symop_operation_xyz")>=0?L.keys.indexOf("_space_group_symop_operation_xyz"):L.keys.indexOf("_symmetry_equiv_pos_as_xyz"); ops=L.rows.map(r=>DC.parseSymop(r[k])); }
  const num=parseInt(g("_space_group_it_number")||g("_symmetry_int_tables_number")||"",10);
  const hm=(g("_space_group_name_h-m_alt")||g("_symmetry_space_group_name_h-m")||"").trim();
  let sgEntry=null; if(!ops){ sgEntry=(hm&&DC.sgBySymbol(hm))||(num>0?DC.sgByNumber(num):null); if(!sgEntry) throw new Error(`CIF '${b.name}': 대칭 연산·공간군이 없습니다`); ops=DC.sgOps(sgEntry); }
  const atoms=[]; const AL=loopOf(b,"_atom_site_fract_x");
  if(AL){ const ix=n=>AL.keys.indexOf(n); const iT=ix("_atom_site_type_symbol"), iL=ix("_atom_site_label"), ixx=ix("_atom_site_fract_x"), iy=ix("_atom_site_fract_y"), iz=ix("_atom_site_fract_z"), io=ix("_atom_site_occupancy"), iu=ix("_atom_site_u_iso_or_equiv"), ib=ix("_atom_site_b_iso_or_equiv");
    for(const r of AL.rows){ const el=DC.elementOf(iT>=0?r[iT]:r[iL])||DC.elementOf(iL>=0?r[iL]:""); if(!el) throw new Error(`CIF '${b.name}': 원소를 모르는 자리 '${iL>=0?r[iL]:"?"}'`);
      const U=iu>=0?cifNum(r[iu]):NaN, B=ib>=0?cifNum(r[ib]):NaN; const occ=io>=0&&isFinite(cifNum(r[io]))?cifNum(r[io]):1;
      atoms.push({el,x:cifNum(r[ixx]),y:cifNum(r[iy]),z:cifNum(r[iz]),occ,B:isFinite(B)?B:isFinite(U)?8*Math.PI*Math.PI*U:0,label:iL>=0?r[iL]:el}); } }
  const formula=(g("_chemical_formula_sum")||g("_chemical_formula_structural")||"").replace(/['"]/g,"").trim();
  const name=(g("_chemical_name_mineral")||g("_chemical_name_common")||g("_chemical_name_systematic")||formula||b.name).replace(/['"]/g,"").trim();
  const codId=g("_cod_database_code")||(/^\d{7}$/.test(b.name)?b.name:"");
  const doi=g("_journal_paper_doi")||"";
  const number=num>0?num:(sgEntry?sgEntry.number:null);
  return {name,formula,cell,ops,atoms,sg:{number,symbol:hm||(sgEntry&&sgEntry.short)||"",system:number?DC.crystalSystem(number):""},source:{codId,doi,block:b.name}}; };
DC.parseCIF=function(text){ const bl=DC.cifBlocks(text).filter(b=>b.items["_cell_length_a"]!=null); if(!bl.length) throw new Error("CIF 에 결정 구조 블록이 없습니다"); return bl.map(DC.phaseFromCIF); };
/* 원자 없이 (DB 수동 입력): 공간군 번호/기호 + 격자상수 */
DC.phaseFromParams=function(p){ const e=p.sgSymbol?DC.sgBySymbol(p.sgSymbol):null; const s=e||DC.sgByNumber(+p.sgNumber); const cell=DC.cell(+p.a,+(p.b||p.a),+(p.c||p.a),+(p.al||90),+(p.be||90),+(p.ga||90));
  return {name:p.name||"",formula:p.formula||"",cell,ops:DC.sgOps(s),atoms:p.atoms||[],sg:{number:s.number,symbol:s.short,system:DC.crystalSystem(s.number),hall:s.hall},source:p.source||{}}; };

/* ---------------- 저장 형식 (PhaseRec v1) ----------------
   { id, name, formula, cell:{a,b,c,al,be,ga} (nm·도), sg:{number, symbol, hall?}, symops?:["x,y,z",…], atoms?:[{label,el,x,y,z,occ,B}],
     source:{kind:"COD"|"DOI"|"paper"|"measured"|"manual", codId?, doi?, citation?, note?}, verified?, note?, (공유 DB 필드: createdBy… version, hidden) }
   symops 가 있으면 그것(CIF 설정 그대로), 없으면 sg.hall → sg.symbol → sg.number 순으로 표의 설정 */
const fmtOp=op=>{ const ax=["x","y","z"]; return [0,1,2].map(i=>{ let s=""; for(let j=0;j<3;j++){ const r=op.R[i*3+j]; if(!r) continue; s+=(r>0?(s?"+":""):"-")+ax[j]; } const t=op.t[i]; if(t){ const n=Math.round(t*24), g=gcd(n,24); s+=(s?"+":"")+(n/g)+"/"+(24/g); } return s||"0"; }).join(","); };
DC.formatSymop=fmtOp;
DC.phaseFromRecord=function(r){ if(!r||!r.cell) throw new Error("결정상 기록에 격자상수가 없습니다");
  const c=r.cell; const cell=DC.cell(+c.a,+(c.b||c.a),+(c.c||c.a),+(c.al||90),+(c.be||90),+(c.ga||90)); let ops=null, e=null;
  if(Array.isArray(r.symops)&&r.symops.length) ops=r.symops.map(DC.parseSymop);
  else { const sg=r.sg||{}; if(sg.hall){ e=sgList().find(x=>x.hallSymbol===sg.hall)||null; } if(!e&&sg.symbol){ try{ e=DC.sgBySymbol(sg.symbol); }catch(_){ e=null; } } if(!e&&sg.number) e=DC.sgByNumber(+sg.number); if(!e) throw new Error(`'${r.name||r.id}': 공간군이 없습니다`);
    // 설정을 명시하지 않았으면 격자 모양으로 고른다: 능면체 R 군 — 능면체 축(a=b=c, α=β=γ≠90) 이면 :R, 단사정 — 90° 가 아닌 각의 축(α→a, β→b, γ→c) 이 고유축
    if(!sg.hall&&!/:/.test(String(sg.symbol||""))) e=DC.settingForCell(e,c)||e; ops=DC.sgOps(e); }
  const number=(r.sg&&+r.sg.number)||(e&&e.number)||null;
  return {id:r.id,name:r.name||r.formula||r.id,formula:r.formula||"",cell,ops,atoms:(r.atoms||[]).map(a=>({label:a.label||a.el,el:a.el,x:+a.x,y:+a.y,z:+a.z,occ:a.occ==null?1:+a.occ,B:a.B==null?0:+a.B})),
    sg:{number,symbol:(r.sg&&r.sg.symbol)||(e&&e.short)||"",system:number?DC.crystalSystem(number):""},source:r.source||{},rec:r}; };
DC.settingForCell=function(e,c){ const n=e.number; const al=+(c.al||90), be=+(c.be||90), ga=+(c.ga||90); const a=+c.a, b=+(c.b||c.a), cc=+(c.c||c.a); const n90=x=>Math.abs(x-90)>1e-3; const L=sgList().filter(x=>x.number===n); if(L.length<2) return e;
  if([146,148,155,160,161,166,167].includes(n)){ const rh=Math.abs(a-b)<=1e-4*a&&Math.abs(b-cc)<=1e-4*a&&Math.abs(al-be)<1e-3&&Math.abs(be-ga)<1e-3&&n90(al); return L.find(x=>x.choice===(rh?"R":"H"))||e; }
  if(n>=3&&n<=15){ const ax=n90(ga)&&!n90(al)&&!n90(be)?"c":n90(al)&&!n90(be)&&!n90(ga)?"a":"b"; if(String(e.choice||"").replace("-","")[0]===ax) return e; return L.find(x=>String(x.choice||"")===ax+"1")||L.find(x=>String(x.choice||"").startsWith(ax))||e; }
  return e; };
DC.recordFromPhase=function(p,extra){ const c=p.cell; const r={name:p.name||"",formula:p.formula||"",cell:{a:c.a,b:c.b,c:c.c,al:c.al,be:c.be,ga:c.ga},sg:{number:p.sg&&p.sg.number||null,symbol:p.sg&&p.sg.symbol||""},
    symops:p.ops?p.ops.map(fmtOp):undefined,atoms:(p.atoms||[]).map(a=>({label:a.label,el:a.el,x:a.x,y:a.y,z:a.z,occ:a.occ,B:a.B})),source:Object.assign({},p.source||{})};
  if(r.source&&r.source.codId&&!r.source.kind) r.source.kind="COD"; return Object.assign(r,extra||{}); };
/* 입력 검사: errors(저장 불가) · warnings(확인 필요) */
DC.validateRecord=function(r){ const E=[], W=[]; const c=r&&r.cell||{}; const a=+c.a, b=+(c.b||c.a), cc=+(c.c||c.a), al=+(c.al||90), be=+(c.be||90), ga=+(c.ga||90);
  if(!r||!String(r.name||r.formula||"").trim()) E.push("이름 또는 화학식이 필요합니다");
  if(!(a>0&&b>0&&cc>0)) E.push("격자상수 a·b·c 는 양수여야 합니다");
  for(const [k,v] of [["a",a],["b",b],["c",cc]]) if(v>0&&(v<0.15||v>2.5)) W.push(`${k} = ${v} nm — 흔한 범위(0.15–2.5 nm) 밖입니다. Å 로 넣지 않았는지 확인하세요 (nm 단위)`);
  for(const [k,v] of [["α",al],["β",be],["γ",ga]]) if(!(v>0&&v<180)) E.push(`${k} 는 0–180° 사이여야 합니다`);
  let n=r&&r.sg&&+r.sg.number; if(!n&&r&&r.sg&&r.sg.symbol){ try{ n=DC.sgBySymbol(r.sg.symbol).number; }catch(_){ E.push(`공간군 기호 '${r.sg.symbol}' 를 모릅니다`); } }
  if(!n&&!(r&&r.symops&&r.symops.length)) E.push("공간군(번호 또는 기호)이 필요합니다");
  if(n&&(n<1||n>230)) E.push("공간군 번호는 1–230 입니다");
  if(n&&n>=1&&n<=230&&a>0&&b>0&&cc>0){ const sys=DC.crystalSystem(n), eq=(x,y)=>Math.abs(x-y)<=1e-4*Math.max(x,y), r90=x=>Math.abs(x-90)<1e-3, r120=x=>Math.abs(x-120)<1e-3;
    const rhombAxes=r.sg&&/:r$/i.test(String(r.sg.symbol||"")); let bad="";
    if(sys==="Cubic"&&!(eq(a,b)&&eq(b,cc)&&r90(al)&&r90(be)&&r90(ga))) bad="입방정: a = b = c, α = β = γ = 90°";
    else if(sys==="Tetragonal"&&!(eq(a,b)&&r90(al)&&r90(be)&&r90(ga))) bad="정방정: a = b, α = β = γ = 90°";
    else if(sys==="Orthorhombic"&&!(r90(al)&&r90(be)&&r90(ga))) bad="사방정: α = β = γ = 90°";
    else if(sys==="Hexagonal"&&!(eq(a,b)&&r90(al)&&r90(be)&&r120(ga))) bad="육방정: a = b, α = β = 90°, γ = 120°";
    else if(sys==="Trigonal"&&!rhombAxes&&!(eq(a,b)&&r90(al)&&r90(be)&&r120(ga))&&!(eq(a,b)&&eq(b,cc)&&eq(al,be)&&eq(be,ga))) bad="삼방정: 육방 축(a = b, γ = 120°) 또는 능면체 축(a = b = c, α = β = γ)";
    else if(sys==="Monoclinic"&&[r90(al),r90(be),r90(ga)].filter(Boolean).length<2) bad="단사정: 세 각 중 둘이 90°";
    if(bad) E.push(`공간군 #${n} (${sys}) 과 격자상수가 맞지 않습니다 — ${bad}`); }
  for(const at of (r&&r.atoms)||[]){ if(!DC.elementOf(at.el||at.label)) E.push(`원자 '${at.label||at.el}': 원소를 모릅니다`); if(at.occ!=null&&!(+at.occ>0&&+at.occ<=1.0001)) E.push(`원자 '${at.label||at.el}': 점유율은 0–1 입니다`); for(const k of ["x","y","z"]) if(!isFinite(+at[k])) E.push(`원자 '${at.label||at.el}': ${k} 좌표가 숫자가 아닙니다`); }
  if(r&&(!r.atoms||!r.atoms.length)) W.push("원자 좌표가 없습니다 — 소멸칙은 공간군으로만, 세기(구조인자)는 계산하지 않습니다");
  if(r&&(!r.source||!(r.source.codId||r.source.doi||r.source.citation||r.source.note))) W.push("출처(COD id·DOI·논문)가 없습니다");
  if(!E.length&&r){ try{ DC.phaseFromRecord(r); }catch(e){ E.push(String(e.message||e)); } }
  return {errors:E,warnings:W,ok:!E.length}; };
/* 중복 후보: 같은 COD id, 또는 같은 화학식 원소 + 같은 공간군 번호 + 격자상수 ±0.5 % */
DC.sameRecord=function(p,q){ if(p.source&&q.source&&p.source.codId&&p.source.codId===q.source.codId) return true;
  const fp=DC.parseFormula(p.formula), fq=DC.parseFormula(q.formula); const ka=Object.keys(fp).sort().join(), kb=Object.keys(fq).sort().join(); if(!ka||ka!==kb) return false;
  if(!(p.sg&&q.sg&&+p.sg.number&&+p.sg.number===+q.sg.number)) return false; const rel=(x,y)=>Math.abs(x-y)/Math.max(x,y)<=0.005;
  const cp=p.cell, cq=q.cell; return rel(+cp.a,+cq.a)&&rel(+(cp.b||cp.a),+(cq.b||cq.a))&&rel(+(cp.c||cp.a),+(cq.c||cq.a)); };

/* ---------------- 화학식 ---------------- */
DC.parseFormula=function(f){ const out={}; const s=String(f||"").replace(/\s+/g,""); const re=/([A-Z][a-z]?)(\d*\.?\d*)/g; let m; while((m=re.exec(s))){ if(!m[1]) continue; out[m[1]]=(out[m[1]]||0)+(m[2]?+m[2]:1); } return out; };

/* ---------------- 정대축의 2D 역격자 ---------------- */
/* [uvw] 에 수직인 정수 hkl 평면 격자의 원시 기저 p, q (cross(p,q) = ±[uvw]/gcd) */
DC.planeBasis=function(uvw){ const z=DC.reduceInt(uvw); const M=Math.max(1,...z.map(Math.abs)); const L=[];
  for(let h=-M;h<=M;h++) for(let k=-M;k<=M;k++) for(let l=-M;l<=M;l++){ if(!h&&!k&&!l) continue; if(h*z[0]+k*z[1]+l*z[2]) continue; L.push([h,k,l]); }
  L.sort((a,b)=>dot(a,a)-dot(b,b)); const zz=dot(z,z);
  for(let i=0;i<L.length;i++) for(let j=i+1;j<L.length;j++){ const c=cross(L[i],L[j]); if(dot(c,c)===zz) return [L[i],L[j]]; }
  return null; };
/* 정수 2D 점들이 만드는 격자의 기저 (Hermite 꼴): (a,b), (0,d) */
function hnf2(pts){ let a=0,b=0,d=0; for(const v of pts){ let x=v[0], y=v[1];
    while(x){ const q=Math.floor(a/x); a-=q*x; b-=q*y; let t=a; a=x; x=t; t=b; b=y; y=t; }
    d=gcd(d,y); if(a<0){ a=-a; b=-b; } if(d) b=((b%d)+d)%d; }
  return a&&d?[[a,b],[0,d]]:null; }
/* 정대축 [uvw] 의 0 차 Laue 반사가 만드는 2D 격자를 REDI 규약으로: |g1| ≤ |g2|, g1·g2 ≤ 0 (둔각), g3 = g1+g2 (짧은 대각선).
   반사 = 공간군 소멸이 아니고 |F|² > Fmin·최대 인 것 (구조상 0 인 Si 222 같은 것 제외). 이중 회절 자리는 이 격자의 점이므로 저절로 들어간다.
   opt.refl: 미리 구한 DC.reflections(…,{F:true}) 를 넘기면 다시 계산하지 않는다 */
DC.zoneBasis=function(phase,uvw,opt){ opt=opt||{}; const z=DC.reduceInt(uvw); const pb=DC.planeBasis(z); if(!pb) return null; const [p,q]=pb;
  const refl=opt.refl||DC.reflections(phase,opt.gmax||12,{F:!!(phase.atoms&&phase.atoms.length)}); let Im=opt.Imax; if(Im==null){ Im=0; for(const r of refl) if(r.I!=null&&r.I>Im) Im=r.I; }
  const thr=(opt.Fmin!=null?opt.Fmin:1e-4)*Im; const pts=[]; const pp=dot(p,p), qq=dot(q,q), pq=dot(p,q), det=pp*qq-pq*pq;
  for(const r of refl){ const h=r.h; if(h[0]*z[0]+h[1]*z[1]+h[2]*z[2]) continue; if(r.I!=null&&Im>0&&r.I<=thr) continue; const hp=dot(h,p), hq=dot(h,q); const m=Math.round((hp*qq-hq*pq)/det), n=Math.round((hq*pp-hp*pq)/det); pts.push([m,n]); }
  if(pts.length<2) return null; const H=hnf2(pts); if(!H) return null;
  let h1=[0,1,2].map(i=>H[0][0]*p[i]+H[0][1]*q[i]), h2=[0,1,2].map(i=>H[1][1]*q[i]);
  const cell=phase.cell; const gv=h=>DC.gvec(cell,h); const add=(a,b,k)=>a.map((x,i)=>x+k*b[i]);
  for(let it=0;it<100;it++){ let g1=gv(h1), g2=gv(h2); if(dot(g1,g1)>dot(g2,g2)){ const t=h1; h1=h2; h2=t; const tg=g1; g1=g2; g2=tg; } const mu=Math.round(dot(g1,g2)/dot(g1,g1)); if(!mu) break; h2=add(h2,h1,-mu); }
  if(dot(gv(h1),gv(h2))>0) h2=h2.map(x=>-x);
  const h3=add(h1,h2,1); const g1=DC.gLen(cell,h1), g2=DC.gLen(cell,h2), g3=DC.gLen(cell,h3);
  return {zone:z,h1,h2,h3,g:[g1,g2,g3],d:[1/g1,1/g2,1/g3],ang12:DC.angleG(cell,h1,h2),ang13:DC.angleG(cell,h1,h3),ang23:DC.angleG(cell,h2,h3),ratio21:g2/g1,ratio31:g3/g1,nSpots:pts.length}; };
/* 정대축 표: |u|,|v|,|w| ≤ maxIdx 인 방향을 Laue 대칭으로 묶어 대표 하나씩 (같은 무리의 방향은 같은 회절도) */
DC.zoneTable=function(phase,opt){ opt=opt||{}; const N=opt.maxIdx||5; const gmax=opt.gmax||12; const refl=DC.reflections(phase,gmax,{F:!!(phase.atoms&&phase.atoms.length)}); let Im=0; for(const r of refl) if(r.I!=null&&r.I>Im) Im=r.I;
  const rots=phase.ops?DC.lauePoints(phase.ops):[[1,0,0,0,1,0,0,0,1]]; const seen=new Set(); const out=[];
  const key=u=>u.join(","); const canon=u=>{ let best=null; for(const R of rots){ const v=[R[0]*u[0]+R[1]*u[1]+R[2]*u[2],R[3]*u[0]+R[4]*u[1]+R[5]*u[2],R[6]*u[0]+R[7]*u[1]+R[8]*u[2]]; for(const s of [1,-1]){ const w=v.map(x=>x*s); if(!best||cmpU(w,best)>0) best=w; } } return best; };
  for(let u=-N;u<=N;u++) for(let v=-N;v<=N;v++) for(let w=-N;w<=N;w++){ if(!u&&!v&&!w) continue; const z=DC.reduceInt([u,v,w]); if(z.some((x,i)=>x!==[u,v,w][i])) continue; const c=canon(z); const k=key(c); if(seen.has(k)) continue; seen.add(k);
    const zb=DC.zoneBasis(phase,c,{refl,Imax:Im,Fmin:opt.Fmin}); if(zb) out.push(zb); }
  out.sort((a,b)=>Math.max(...a.zone.map(Math.abs))-Math.max(...b.zone.map(Math.abs))||a.zone.reduce((s,x)=>s+x*x,0)-b.zone.reduce((s,x)=>s+x*x,0));
  return out; };
const cmpU=(a,b)=>{ // 보기 좋은 대표: 양수 많고, 큰 성분 먼저
  const pa=a.filter(x=>x>=0).length, pb=b.filter(x=>x>=0).length; if(pa!==pb) return pa-pb; for(let i=0;i<3;i++) if(a[i]!==b[i]) return a[i]-b[i]; return 0; };

if(typeof module!=="undefined"&&module.exports) module.exports=DC;
G.DIFF_CORE=DC;
})();
