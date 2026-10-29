# Proposed xps_core.js changes from the CasaXPS Round-4/5/6 readings (SPEC_casa_v4.md)

**Applied in v1.9.2** (2026-10-29, `casa_v4_apply.py --live`, all anchors matched; follow-ups in SPEC_casa_v4.md §7). Originally: not applied (another agent was editing `src/xps_core.js`, LM stall rule). Anchors are exact text of the 2026-10-29 snapshot (after casa_v3_patch).
`test/ref/casa_v4_apply.py <copy>` applies §1–§6 (§4/§6 = Rounds 5–6 ST shape and position) to a copy by these anchors (refuses the live file without `--live`); `node test/ref/check_casa_v4.js <copy>` prints the result.
Verified on a full module copy (all of §1–§6): `test_shapes2.js` 73/73; `test_casa_v3.js` (12) / `test_xps_casa.js` (3) fail only on the expectations listed in §7.

---------------------------------------------------------------------------------------------------
## 1. `casaWidthAsym0` — one rule for the CasaXPS fwhm number and asymmetry index  [PINNED, 14 shapes]
Replace the comment block starting `/* CasaXPS 'fwhm' 숫자와 비대칭 지수 (SPEC_casa_v3 §3` up to the constants line, the line
`const CASA_DELTA_LA=0.0025, CASA_DELTA_TAIL=0.011, CASA_DELTA_ST=0.0054, CWA=new Map();` → `const CASA_GRID=1/40, CASA_DS_T=600, CWA=new Map();`
and the whole `function casaWidthAsym0(c){…}` by (keep `refineHalfMax` and the cached wrapper `casaWidthAsym`):
```js
/* CasaXPS 'fwhm' number and asymmetry index (SPEC_casa_v4 §3): x̂ = parabola vertex on nodes k·d (d = F/40, DS wL/40; ST: d = D = base FWHM/40)
   around the largest node; half level = vertex value / 2; fwhm = 2·min(l, h), asym = 1 − l/h. node = k·d (CasaXPS ST position, §4). */
function casaWidthAsym0(c){
  const md=compMod(c), sp=shapeParams(c,1), Fb=c.shape==="V"?sp.wL:sp.F;
  // ST: CasaXPS tabulates S(μ)F on nodes k·D, D = (base FWHM)/40, with the step lagged by D/2, and puts its position on the largest node (SPEC_casa_v4 §4)
  const D=md==="ST"?(c.stCasaK>0?c.stCasaK:casaStK(c))*Fb/40:0; if(md==="ST") sp.stLag=D/2;
  const w0=trueFWHM(c.shape,sp);
  if(w0.fwhm==null||w0.lo==null||w0.hi==null) return {fwhm:null,asym:null,exact:false,lo:null,hi:null,xmax:null,node:null,side:null,note:"no half-maximum points (CasaXPS fwhm not defined here)"};
  const W=w0.fwhm, F=(c.shape==="DS"?sp.wL:sp.F)||W, d=md==="ST"?D:CASA_GRID*F, ctx={lo:-60*W,hi:60*W};   // 같은 ctx → 창 정규화 상수는 상쇄
  const ev=xs=>shapeUnit(c.shape,Float64Array.from(xs),sp,ctx);
  let k=Math.round(w0.xmax/d), v=ev([(k-1)*d,k*d,(k+1)*d]);
  for(let it=0;it<16;it++){ if(v[2]>v[1]) k++; else if(v[0]>v[1]) k--; else break; v=ev([(k-1)*d,k*d,(k+1)*d]); }
  const den=v[0]-2*v[1]+v[2], xs=den<0?(k+0.5*(v[0]-v[2])/den)*d:k*d, half=(den<0?v[1]-0.125*(v[0]-v[2])*(v[0]-v[2])/den:v[1])/2;
  const cross=(u,w)=>{ for(let pass=0;pass<3;pass++){ const K=64, g=new Float64Array(K+1); for(let i=0;i<=K;i++) g[i]=u+(w-u)*i/K; const f=ev(g); let j=-1;
      for(let i=0;i<K;i++) if((f[i]-half)*(f[i+1]-half)<=0){ j=i; break; } if(j<0) return null; if(pass===2) return g[j]+(f[j]-half)/(f[j]-f[j+1])*(g[j+1]-g[j]); u=g[j]; w=g[j+1]; } return null; };
  const lo=cross(w0.lo-0.03*W,Math.min(w0.lo+0.03*W,xs)), hi=cross(Math.max(w0.hi-0.03*W,xs),w0.hi+0.03*W);
  if(lo==null||hi==null) return {fwhm:null,asym:null,exact:false,lo:null,hi:null,xmax:null,node:null,side:null,note:"half-maximum points not found"};
  const l=xs-lo, h=hi-xs, sym=Math.abs(h-l)<=1e-7*(h+l);
  const stConv=md==="ST"&&(sp.wG>1e-6||sp.stG>1e-9);
  const note=md==="ST"?"ST: CasaXPS shape = S(μ)F with the step lagged by (base FWHM)/80; position = INNO origin + largest node k·(base FWHM)/40 (measured, 16 readings ≤ 0.2 meV); fwhm = 2·min(l,h), asym = 1 − l/h"+(stConv?" (Gaussian-broadened base: −0.5…+0.3 % — approximate)":" (measured ≤ 4e-4 / 1.1e-3)")
    :sym?"CasaXPS fwhm = true FWHM"
    :"CasaXPS fwhm = 2·min(l, h), asym = 1 − l/h; maximum = parabola vertex on nodes k·F/40 (measured: 14 shapes, fwhm ≤ 5.3e-4, asym ≤ 6e-4)";
  return {fwhm:sym?w0.fwhm:2*Math.min(l,h),asym:sym?0:1-l/h,side:l<=h?"low":"high",exact:!stConv,lo,hi,xmax:xs,node:k*d,note};
}
/* ST: base FWHM / base width (F, or wL for V) — scale-free, fixed by the base's shape parameters (m, α, β, n …) */
function casaStK(c){ const b=clone(c); b.mod=null; delete b.stCasaK; const sk=b.shape==="V"?"wL":"F"; b[sk]=1; const w=trueFWHM(b.shape,shapeParams(b,1)).fwhm; return w>0?w:1; }
```
- Cost: ≤ ~22 `shapeUnit` calls of ≤ 65 points (cached per shape in `CWA`); the check script's 17 fits run in ~1.5 s.
- `xmax` = CasaXPS's x̂, `node` = k·d (ST position, §4), both relative to the INNO origin in eV; `casaStK` is §6. `side` is new ("low"/"high": which half width is the fwhm).
- Results (`check_casa_v4.js`, INNO fit → compToCasa vs CasaXPS): LA(2.43,1.53,0) 0.866770/0.866617, asym −0.2950/−0.2944; LA(1.53,2.43,0) F1·F2 +1.3e-4·+5.2e-4, 0.2278/0.2275;
  LA(2,3,0) +5.3e-4, 0.1939/0.1936; LA(2,3,30) +2.2e-4, 0.1302/0.1296; GL(30)T(2) +1.2e-4, 0.2137/0.2138; GL(30)T(1)·GL(0)T(2) asym 0.4151·0.1972 (0.4151·0.1973);
  TLA ±3.8e-4; DS −3.7e-4…+1.3e-4; ST(1.7)LA(1,1,0) +3.3e-4, 0.1529/0.1529; ST(1.7)GL(30) +1.1e-3, −0.0112/−0.0131; ST(3)LA(1,1,200) +2.3e-3, 0.0412/0.0440; ST(1.7)LA(1,1,200) −5.7e-3, 0.0184/0.0071.

## 2. New `dsCasaAreaFactor(c, emin, emax)` (after `casaWidthAsym0`)  [PINNED]
```js
/* CasaXPS DS area = amplitude × ∫ P over E0 ± 600·wL (SPEC_casa_v4 §2); INNO DS area = area inside [emin, emax] → CasaXPS area = INNO area × factor.
   Raw-DS integrals (closed form); the n-Gaussian changes the window integral by ≤ 5e-5. Cut variant (dsl): null (not read). */
function dsCasaAreaFactor(c,emin,emax){ if(!c||c.shape!=="DS"||c.dsl>0||!(emax>emin)) return null; const sp=shapeParams(c,1), g=sp.wL/2, al=c.alpha, T=CASA_DS_T*sp.wL;
  const N600=dsPrim(T,al,g)-dsPrim(-T,al,g), Nw=dsPrim(emax-c.be,al,g)-dsPrim(emin-c.be,al,g); return Nw>0?N600/Nw:null; }
```
(emin/emax = absolute BE of the INNO fit window, i.e. `prob.E[0]`, `prob.E[n−1]` — the same window `compEval` normalises over.)
Check: INNO fit × factor vs CasaXPS Area — DS(0,100) +3.3e-4, DS(0.2,100) +2.3e-4, DS(0.1,0) −7e-5, DS(0.1,100) +2.7e-4, DS(0.1,100) on 280.03–290.08 +1.9e-4 (factor 1.03306 · 1.96612 · 1.50361 · 1.45068 · 1.59584).

## 3. DS in `compToCasa` / `casaToComp`
- `function compToCasa(c){` → `function compToCasa(c,o){\n  o=o||{};` (optional `o.emin`, `o.emax` = fit window).
- DS branch head: append `, ak=o.emin!=null&&o.emax!=null?dsCasaAreaFactor(c,+o.emin,+o.emax):null;` to the `const n=…, cw=…, sh=…` line.
- DS `base={…}`: `posShift:sh,exact:false,note:…` → `posShift:sh,areaFactor:ak,area:ak?c.area*ak:null,exact:!!ak&&n.ok,note:"DS(α,n…), Cookbook p. 66: Lorentzian FWHM = wL; shape, tail side and n rule measured (α 0, 0.1, 0.2; n 0, 100); "+cw.note+"; CasaXPS position = INNO origin + wL·|cot(π/(2−α))| (measured ≤ 1.6e-4 eV)"+…+"; CasaXPS area = amplitude × ∫ over ±600·wL"+(ak?" = INNO window area × "+ak.toFixed(5):" — pass the fit window {emin, emax} to convert …")`.
- `casaToComp`, DS parse: note → `"DS(α,n): shape, tail side and n rule (Gaussian FWHM = n·wL/144.87) measured (α 0–0.2, n 0, 100)"` (drop 'area not pinned … ≈ 1.43 ×').
- `casaToComp`, DS position note: `(measured, one reading DS(0.1,100))` → `(measured α 0, 0.1, 0.2; ≤ 1.6e-4 eV)`.
- `casaToComp`, before the final `delete r.keepF; return r; }`:
```js
  // DS area: CasaXPS Area (amplitude × ∫ over ±600·wL) → INNO area inside the fit window [o.emin, o.emax] (SPEC_casa_v4 §2)
  if(c.shape==="DS"&&!o.innoUnits){ const k=o.emin!=null&&o.emax!=null?dsCasaAreaFactor(c,+o.emin,+o.emax):null;
    if(k){ c.area=area/k; r.areaFactor=k; r.exact=Number.isInteger(+c.casaN)&&!(c.dsl>0); r.note+="; area: INNO window area = CasaXPS area / "+k.toFixed(5)+" (CasaXPS DS area = amplitude × ∫ over ±600·wL)"; }
    else r.note+="; area kept as given — pass the fit window {emin, emax} to convert the CasaXPS DS area (CasaXPS area = amplitude × ∫ over ±600·wL)"; }
```
  (must run after the DS position shift, since the factor depends on `c.be`.) Check: `casaToComp('DS(0.2,100)',{be:285.0927,fwhm:0.795054,area:27119.4,emin:280.03,emax:297.13})` → area 13792.3 vs INNO fit 13796.6 (−3.1e-4); round trip compToCasa(…,{emin,emax}).area = 27119.4.
- `casaModelToINNO`: DS component areas are CasaXPS areas — convert with `dsCasaAreaFactor(comp, BE(region KE end), BE(region KE start))` per region (same for area bounds;
  area links between a DS and a non-DS component then carry the factor). Not in the apply script (needs the region→window mapping there).

## 4. ST position  [PINNED — Round 6, 16 readings ≤ 0.16 meV]
CasaXPS pos = INNO (base) origin + D·k*, D = K·F/40 with K = base-shape FWHM / F (`casaStK`), k* = largest node of CasaXPS's ST table (nodes k·D aligned to the base origin, no sub-node refinement).
`casaWidthAsym0` (§1 code, as applied by the script) returns `node:k*D` for every shape; for ST it uses D as its grid and sets `sp.stLag=D/2` (§6).
- `casaToComp`: `  else r.posShift=mod&&mod.n==="ST"?null:0;` →
```js
  else if(mod&&mod.n==="ST"&&!o.innoUnits){ const cw=casaWidthAsym(c); if(cw.node!=null){ c.be=be-cw.node; r.posShift=cw.node; r.note+="; position: INNO origin = CasaXPS position − "+cw.node.toFixed(5)+" eV (largest node of CasaXPS's ST table, step (base FWHM)/40)"; } else r.posShift=null; }
  else r.posShift=0;
```
- `compToCasa` ST wrapper: `fwhm:cw.fwhm,asym:cw.asym,fwhmExact:false,posShift:null,exact:false,` → `fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,posShift:cw.node,exact:cw.exact&&base.exact,`;
  note tail `; CasaXPS position offset not pinned (…+0.0984 eV); "+base.note` → `; "+base.note` (cw.note carries the rule). `cw.exact` is false only for a Gaussian-broadened base (LA/GL n > 0 or ST γ > 0).
- casaToComp ST note: `CasaXPS position and fwhm number not pinned (…) — refit position and width; ` → `CasaXPS draws it with the step lagged by (base FWHM)/80 and positions it on the largest node of a (base FWHM)/40 grid (SPEC_casa_v4 §4); `.
- k* jumps by one node when the width crosses a half-node: this is CasaXPS's own (discontinuous) behaviour; INNO's fitted origin stays continuous, only the reported CasaXPS position steps.
- Effect (CasaXPS reading → INNO origin vs INNO fit): 0.07–0.32 eV → ≤ 1.5e-4 eV for 14 n = 0 readings, 2e-4 / 6.2e-4 eV for the two LA(1,1,200) bases.

## 6. ST shape: CasaXPS's step lags by D/2 = (base FWHM)/80  [PINNED for bases without ⊗G — Rounds 5–6]
CasaXPS evaluates S(μ)F = F·[Φ + A/μ] with Φ an exclusive rectangle sum on the D grid, i.e. Φ(x − D/2); area 1 = the numerical area of that shape.
(12 GL-based files + LA(1,1,0): CasaXPS RMS reproduced (0.95–1.07×), area ≤ 7.6e-5, fwhm ≤ 3.3e-4, asym ≤ 1.0e-3; LA(1,1,200) bases −4.9e-3/+2.9e-3, area −1.1e-3/−1.7e-3 stay approximate.)
- `s2Eval`: `const muI=mod==="ST"?1/Math.max(1e-9,p.stMu):0, stN=1/(0.5+muI);` →
  `const muI=mod==="ST"?1/Math.max(1e-9,p.stMu):0, L=mod==="ST"&&p.stLag>0?p.stLag:0; let stN=1/(0.5+muI);`
- `s2Eval`, end of the `if(mod){ … }` Φ block (after the `Phi[i0]=B.left; …` line's closing `}`, before the block's own `}`):
```js
    if(L){ const P0=Phi.slice(); let dz=0; for(let i=0;i<M;i++){ Phi[i]=B.cdf?B.cdf(gx[i]-L):lag6(P0,k0,h,gx[i]-L); dz+=v[i]*(Phi[i]-P0[i]); } stN=1/(0.5+muI+dz*h); }
```
  (v = unconvolved base on the fine grid there; dz·h = ∫F(Φ(x−L) − Φ(x)) keeps area 1.)
- `s2Eval`, direct path: `const ph=B.cdf?B.cdf(x[i]):lag6(Pf,kf,hf,x[i]);` → `const ph=B.cdf?B.cdf(x[i]-L):lag6(Pf,kf,hf,x[i]);` (Pf is already lagged).
- `shapeParams`: after the `p.stG=…` assignment add `if(md==="ST"&&c.stCasaK>0) p.stLag=c.stCasaK*(c.shape==="V"?p.wL:p.F)/80;`.
- New helper (after `casaWidthAsym0`): `function casaStK(c){ const b=clone(c); b.mod=null; delete b.stCasaK; const sk=b.shape==="V"?"wL":"F"; b[sk]=1; const w=trueFWHM(b.shape,shapeParams(b,1)).fwhm; return w>0?w:1; }`
  (K = 1 for GL(0) and LA(1,1,0), 0.96136 for GL(30), 1.99129 for LA(1,1,200)); exported.
- `casaToComp`, ST modifier: `c.stG=0; }` → `c.stG=0; if(!o.innoUnits) c.stCasaK=casaStK(c); }` — an imported CasaXPS ST component is fitted with CasaXPS's own shape.
  New comp field `stCasaK` (number, not a fit parameter; fixed by the base's m/α/β/n, so it stays valid when F changes). INNO-native ST components (no stCasaK) keep the Cookbook shape;
  `casaWidthAsym0` always evaluates CasaXPS's shape for the CasaXPS fwhm/asym/position.
- Result (`check_casa_v4.js` on a patched copy, INNO fit with the casaToComp component, 16 ST rows): ST(μ)GL(0) ×9 area ≤ 4.7e-5, fwhm ≤ 2.6e-4, asym ≤ 1e-4, pos ≤ 0.09 meV;
  ST(1.7)GL(30) F 0.6/1.2/2.4 −2.6…−3.3e-4, asym −0.0141/−0.0131, pos ≤ 0.09 meV; ST(1.7)LA(1,1,0) +2.8e-4, pos +0.10 meV.

## 5. `compToCasa0` / GL(m,n) tail flags — widths now measured
- LA asymmetric: `fwhmExact:cw.exact,n,exact:false,note:"asymmetric LA: shape and tail side measured; "+cw.note+(ni>0?" (n > 0 asymmetric LA not read yet)":"")` →
  `fwhmExact:cw.exact,n,exact:Math.abs(n-ni)<0.05,note:"asymmetric LA: shape, tail side and n rule measured (LA(1.53,2.43,0), LA(2.43,1.53,0), LA(2,3,0), LA(2,3,30)); "+cw.note`.
- AGT GL(p)T(k): `fwhmExact:false,exact:false,` → `fwhmExact:cw.exact,exact:true,`.
- `compToCasa` GLN/SGLN: `fwhmExact:!t,exact:!t&&n.ok,` → `fwhmExact:true,exact:n.ok,`.
- `casaToComp0` asymmetric LA note `(n > 0 asymmetric LA not read yet)` may stay or be dropped (now read: LA(2,3,30)); not touched by the apply script.
- Exports: `CASA_DELTA_LA,CASA_DELTA_TAIL,CASA_DELTA_ST,` → `CASA_GRID,CASA_DS_T,dsCasaAreaFactor,` (no other file uses the old names).

## 6. UI (`src/xps_ui.js` l. 679, component card 'CasaXPS:' line)
`ps=q.posShift>1e-6?`, pos INNO + ${f(q.posShift,4)} eV`:…` → prefix `≈ ` when `q.posShiftExact===false` (ST); keep ', pos ≈ (offset not pinned)' for `posShift===null`.
When the card knows the fit window, call `compToCasa(c,{emin,emax})` and show `area ${f(q.area,1)}` for DS (CasaXPS Area).

## 7. Test expectation changes (`test/test_casa_v3.js`, `test/test_xps_casa.js`)
Failures on a patched copy, all expected:
- test_casa_v3 §1(b): `KIND` laKink/tail are now exact → `exOK` must expect `fwhmExact===true` for laKink/tail (only ST false); tolerances FW_TOL {smooth 6e-4, laKink 6e-4, tail 6e-4, ST 6e-4, ST with ⊗G base 6e-3}, ASYM_TOL {smooth 1.2e-3, laKink 1e-3, tail 1e-3, ST 1.2e-3, ST ⊗G 1.2e-2}
  (ST(1.7)LA(1,1,200) now 0.0164 vs 0.0071 and area −1.1e-3 — widen its area tolerance to 2e-3); ST pos now checked at 2e-4 eV (posShift = cw.node).
- test_casa_v3 §1b shape_only: `fwhmExact===false` → `true`; tolerance 2e-3 → 6e-4 (now 0.4151/0.4151, 0.1972/0.1973, 0.2137/0.2138).
- test_casa_v3 §2: ST row 'Casa→INNO F ≤ 0.005' → 6e-3 and add origin ≤ 1e-3 eV (6.2e-4 measured; n = 0 ST rows ≤ 2e-4); `exact` expectation for laKink/tail rows → true (casaToComp LA/AGT now exact when n integer — casaToComp0's `exact:sym&&…` for asymmetric LA stays false unless also changed; decide together).
- test_casa_v3: `casaToComp('DS(0.1,100)', …)` — replace `/not pinned/` by checking `r.posShift` and, with `{area:20811.2, emin:280.03, emax:297.13}`, `r.comp.area ≈ INNO fitted area (≤ 5e-4)`; `ST: … 'not pinned'` → expect `a.posShift === casaWidthAsym(a.comp).node` (k*·D, e.g. 4·1.99129·F/40 for ST(1.7)LA(1,1,200)) and `a.comp.stCasaK ≈ 1.99129`; the 'true FWHM = given width' part only with `innoUnits`.
- test_xps_casa §5: V2_LA_153_243_F2, Shape_LA153_243, V2_GL30T2_F1.2 — the '≈ empirical δ' expectations → exact (fwhm ≤ 6e-4, asym ≤ 1e-3).
- test_casa_v3 §1(a): K_ST17_LA11200_F1 'INNO fit = scipy reference' now differs because casaToComp sets `stCasaK` (lagged step) — compare against `ref_casa_v4.json` `st5` rows with sh = F/80 instead of ref_casa_v3.
- New rows to add (Round 4/5, from `ref_casa_v4.json` `fits`, `st5`, `ds5`): the Round-4 K files, the 7 Round-5 and 5 Round-6 files (`st6` rows; ST via stCasaK: area ≤ 2e-4, fwhm ≤ 5e-4, asym ≤ 1.2e-3, pos ≤ 2e-4 eV), DS(0.3,0) (area ≤ 5e-4 via factor), the 280.03–290.08 (i1 = 220) refit; DS area via `dsCasaAreaFactor` ≤ 5e-4; DS pos ≤ 2e-4 eV; LA(2.43,1.53,0) fwhm side 'high' and asym sign.
