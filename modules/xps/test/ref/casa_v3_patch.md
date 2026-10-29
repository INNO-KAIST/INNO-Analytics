# Proposed xps_core.js changes from the CasaXPS Round-3 readings (SPEC_casa_v3.md)

Not applied (another agent was editing `src/xps_core.js`). Line numbers are from the 2026-10-29 snapshot and will drift;
each item names the function and quotes the current text to anchor it. Reference numbers: `test/ref/ref_casa_v3.py` → `ref_casa_v3.json`.

---------------------------------------------------------------------------------------------------
## 1. `fitStats` (≈ l. 913) — CasaXPS DoF counts only channels with y ≠ background  [PINNED]

Current:
```js
  let casa=null; if(prob.weight==="casa"){ let s2=0; for(let i=1;i<n-1;i++) s2+=(y[i]-mu[i])*(y[i]-mu[i]); casa={rms:Math.sqrt(s2),dof:Math.max(1,n-2-p),n:n-2}; }
```
Proposed (B = the static region background that `makeEvaluator` caches on `prob._Bst`; active backgrounds fall back to all channels):
```js
  // CasaXPS 'RMS'·DoF (SPEC_casa_v2 §1, SPEC_casa_v3 §2): √Σ(y − μ)² over channels i0+1 … i1−1; DoF = #{y_i ≠ B_i} − p
  // (GL(30,100) 154−3 = 151, GL(30) 159−3 = 156, GL(30)T(2) 340−3 = 337 — channels where the data equal the background exactly are not counted)
  let casa=null; if(prob.weight==="casa"){ const Bc=prob._Bst||null; let s2=0, ne=0;
    for(let i=1;i<n-1;i++){ const d=y[i]-mu[i]; s2+=d*d; if(!Bc||Math.abs(y[i]-Bc[i])>1e-9*Math.max(1,Math.abs(y[i]))) ne++; }
    casa={rms:Math.sqrt(s2),dof:Math.max(1,ne-p),n:ne}; }
```
Check: the K files give 151 / 156 / 337 / 338 exactly with this rule (`casa_dof()` in ref_casa_v3.py).
Leave `covariance()` (l. 893, `n-k-(casa?2:0)`) as is; CasaXPS's error-bar DoF has not been read.
Test: `test_xps_casa.js` — fit K_GL30n100_F1_2 with GLN(m .3, casaN 100), weight "casa", bg linear conv casa aw 1 → `casa.dof === 151`;
Shape_GLP30 GL(30) → 156; V2_GL30T2_F1.2 → 337; every other K file → 338.

---------------------------------------------------------------------------------------------------
## 2. New helper: CasaXPS width number and asymmetry index  [PINNED for smooth-top shapes; empirical δ otherwise]

### 2a. Return the half-maximum points from the FWHM routines
- `trueFWHM` (l. 194): last line → `return {fwhm:a!=null&&b!=null?b-a:null,xmax:x[im],ymax:v[im],lo:a,hi:b};`
- `s2FWHM` (l. 487): `const res={fwhm:A!=null&&Bx!=null?Bx-A:null,xmax:xm,ymax:ym,lo:A,hi:Bx};`
  (`lo` = low-BE half-max point, `hi` = high-BE point, both relative to the component origin; the cache key does not change.)

### 2b. Helper (place in §16 next to `casaWidthFactor`, ≈ l. 1757)
```js
/* CasaXPS 'fwhm' number and asymmetry index (SPEC_casa_v3 §3): x̂ = peak maximum, a < x̂ < b half-maximum points (a on the low-BE side),
   L = x̂ − a, H = b − x̂ → fwhm = 2(L + δ), asym = 1 − (L + δ)/(H − δ).  Symmetric shapes: = true FWHM, 0.
   δ = 0 when the top is smooth (TLA, any shape with n > 0, DS; matches CasaXPS to ≤ 4.4e-4 in fwhm, ≤ 5e-4 in asym).
   Kinked / flat tops: CasaXPS's maximum sits δ towards the tail — empirical, from the readings, not a derived rule:
   LA(α≠β, n = 0) δ = 0.0025·F (LA(2,3,0), LA(1.53,2.43,0)); GL/SGL with T(k)/PHI and AGT tu "fwhm" δ = 0.011·F (GL(30)T(1), GL(30)T(2), GL(0)T(2)).
   ST(μ): not pinned — returns the true FWHM of the ST shape (−0.37 % vs CasaXPS for ST(1.7)LA(1,1,200)) and asym from δ = 0.0054·F (one reading). */
const CASA_DELTA_LA=0.0025, CASA_DELTA_TAIL=0.011, CASA_DELTA_ST=0.0054;
function casaWidthAsym(c){
  const sp=shapeParams(c,1), w=trueFWHM(c.shape,sp); if(w.fwhm==null||w.lo==null||w.hi==null) return {fwhm:null,asym:null,exact:false,note:"no half-maximum points"};
  const F=c.shape==="DS"?sp.wL:sp.F, L=w.xmax-w.lo, H=w.hi-w.xmax, md=compMod(c);
  const flatTail=!!compTail(c)||(c.shape==="AGT"&&c.tu==="fwhm"), laKink=c.shape==="LA"&&Math.abs(c.alpha-c.beta)>1e-9&&!(sp.wG>1e-6);
  const d=md==="ST"?CASA_DELTA_ST*F:flatTail?CASA_DELTA_TAIL*F:laKink?CASA_DELTA_LA*F:0;
  const sym=Math.abs(H-L)<=1e-9*(H+L);
  const fwhm=md==="ST"?w.fwhm:2*(L+d), asym=sym?0:1-(L+d)/(H-d);
  return {fwhm,asym,exact:!(md||flatTail||laKink),
    note:md==="ST"?"ST: CasaXPS fwhm ≈ true FWHM (−0.4 %), position offset not pinned":flatTail||laKink?"CasaXPS fwhm = 2(x̂−a) with an empirical peak offset (≈1–2 %)":"CasaXPS fwhm = 2 × (maximum − low-BE half-maximum point)"};
}
```
Expected values (unit F; from ref_casa_v3.json, CasaXPS number / F in brackets):
TLA(1.2,3,0) 0.795678 (0.795642) · TLA(1.2,3,100) 1.180252 (1.180775) · DS(0.1,100) 1.422235 per wL (1.422827) ·
LA(2,3,0) 0.509941+0.005 = 0.514941 (0.514509) · LA(1.53,2.43,0) 0.579666 (0.580096 / 0.579871) · GL(30)T(2) 0.962115+0.022 = 0.984115 (0.983758) ·
asym: TLA 0.1529 / 0.0822, DS 0.1180, LA(2,3) 0.1947, LA(1.53,2.43) 0.2283, GL(30)T(2) 0.2141, GL(30)T(1) ≈ 0.4151, GL(0)T(2) ≈ 0.1973.
Side caveat: all tested shapes have the tail on the high-BE side, so 'low-BE side' and 'steep side' coincide — reading R3 (LA(2.43,1.53,0)) decides; until then use the low-BE side.

---------------------------------------------------------------------------------------------------
## 3. `compToCasa` (l. 1823) / `compToCasa0` (l. 1759)

Add `asym` and `posShift` (= CasaXPS position − INNO `be`) to every returned object (default `asym:0, posShift:0`).
- `compToCasa0`, LA branch, `if(!sym) return {…fwhm:null…}` → `const cw=casaWidthAsym(c); return {notation:not,fwhm:cw.fwhm,asym:cw.asym,n,exact:false,note:"asymmetric LA: "+cw.note};`
  (with n > 0 the top is smooth → `exact:cw.exact` once an n > 0 asymmetric LA is read, R3.)
- `compToCasa0`, AGT (GL(p)T(k)) branch: `fwhm:casaWidthAsym(c).fwhm` (≈ 0.984·F for p 30, k 2), keep `exact:false`.
- `compToCasa`, GLN/SGLN branch: no tail → `exact:n.ok` (n rule measured: GL(30,100) fwhm −1.2e-4, SGL(30,100) −9e-5, area 3e-5 / 5.7e-4);
  tail → `fwhm:casaWidthAsym(c).fwhm, asym:…, exact:false`. Change note "Gaussian FWHM = n·F/144.87 assumed (LA rule)" → "… (measured: LA(1,1,n), GL(30,100), SGL(30,100))".
- TLA branch: `const cw=casaWidthAsym(c); base={notation:`TLA(${…},${…},${n.t})`,fwhm:cw.fwhm,asym:cw.asym,exact:n.ok,note:"TLA(α,μ,n), Cookbook p. 66 — notation, μ, tail side and n rule measured (TLA(1.2,3,0), TLA(1.2,3,100)); fwhm = 2·(max − low-BE half point)"};`
- DS branch: `const cw=casaWidthAsym(c), sh=c.alpha>0?sp.wL*Math.abs(1/Math.tan(PI/(2-c.alpha))):0;`
  `base={…,fwhm:cw.fwhm,asym:cw.asym,posShift:sh,exact:false,note:"DS(α,n): shape, tail side and n rule measured; CasaXPS position = INNO origin + wL·|cot(π/(2−α))| (one reading, DS(0.1,100)); CasaXPS area normalisation not pinned (INNO area = area inside the window)"};`
  (`dsl > 0` variant: same shift assumed, note it.)
- ST wrapper (`if(md==="ST")`): `fwhm: casaWidthAsym(c).fwhm` (true FWHM of the modified shape), `posShift:null`, note "position offset not pinned (ST(1.7)LA(1,1,200): CasaXPS pos − INNO origin = +0.0984 eV)".
  SB wrapper unchanged.

---------------------------------------------------------------------------------------------------
## 4. `casaToComp` (l. 1792)

(a) Width mapping. Replace the block
```js
  if(!r.keepF){ const u=unitFWHM(c); if(c.shape==="DS") c.wL=u>0?fw/u:fw; else c.F=u>0?fw/u:fw;
    const asym=(c.shape==="GLN"||c.shape==="SGLN")&&compTail(c)||c.shape==="TLA"||c.shape==="DS"; if(asym) r.note+="; CasaXPS fwhm taken as the true FWHM (not verified for asymmetric shapes) — refit the width"; }
```
by (the modifier is attached first so that ST uses its own width; move the `if(mod){…}` block above this one):
```js
  if(!r.keepF){ const sk=c.shape==="DS"?"wL":"F", u0=clone(c); u0[sk]=1; const cw=casaWidthAsym(u0), u=cw.fwhm>0?cw.fwhm:unitFWHM(c);
    c[sk]=u>0?fw/u:fw; if(!cw.exact) r.note+="; "+cw.note+" — refit the width"; }
```
- GL/SGL tail path: delete `if(tl){ r.comp.F=fw; r.keepF=true; }` so the tail goes through (a) (currently F = CasaXPS number: −1.6 % for GL(30)T(2); with (a) −0.04 %).
- Today TLA(1.2,3,0) gets F 8 % too small (true FWHM used: 0.911717/0.867476 = 1.051 vs 1.146) and ST(1.7)LA(1,1,200) 4.5 % too small (base width used, modifier dropped by `unitFWHM`).
- `casaToComp0` LA asymmetric branch: `c.F=fw/casaWidthFactor(...)` → `const u0=clone(c); u0.F=1; c.F=fw/casaWidthAsym(u0).fwhm;` (note stays exact:false).

(b) DS position. After the width is set (DS branch):
```js
  if(c.shape==="DS"&&c.alpha>0){ c.be=be+c.wL/Math.tan(PI/(2-c.alpha)); r.note+="; position: INNO origin = CasaXPS position − wL·|cot(π/(2−α))| (measured DS(0.1,100))"; }   // tan < 0 → shift towards lower BE
```
(c) ST position: no change possible yet; add to the note "CasaXPS ST position is 0.0984 eV (0.197·F) above the INNO origin for ST(1.7)LA(1,1,200) — offset rule not pinned; refit the position".

(d) Notes: TLA "notation assumed by analogy" → "TLA(α,μ,n) notation accepted by CasaXPS 2.3.27; n rule measured"; GLN/SGLN "assumed (LA rule, measured for LA(1,1,n) only)" → "measured for GL(30,100), SGL(30,100)"; DS "Gaussian FWHM = n·wL/144.87 assumed" → "measured (DS(0.1,100))".
    Set `exact:true` for GL(m,n)/SGL(m,n) without tail and for TLA(α,μ,n) (width via casaWidthAsym), keep false for DS (area), tails, ST, asymmetric LA n = 0.

---------------------------------------------------------------------------------------------------
## 5. `casaModelToINNO` (l. 1855)
- l. 1870: `const k=c.fwhm&&c.fwhm.v>0&&r.comp.F>0?r.comp.F/c.fwhm.v:1;` → use `const W=r.comp.shape==="DS"?r.comp.wL:r.comp.F; const k=c.fwhm&&c.fwhm.v>0&&W>0?W/c.fwhm.v:1;`
  (DS has no F, so a DS fwhm link is currently scaled by 1).
- Position links/offsets between components whose `posShift` differs (DS with different α or wL, ST): the CasaXPS offset applies to CasaXPS positions;
  convert with `be_INNO = pos_Casa − posShift` per component and add a note when a linked partner has a different shift.

---------------------------------------------------------------------------------------------------
## 6. Tests to add (`test/test_xps_casa.js`, new §6 "Round 3")
- For each Round-3 row of ref_casa_v3.json: INNO fit of the K file (weight casa, Linear conv casa aw 1, casaToComp shape) reproduces sim `area`, `pos`, `rms` (≤ 1e-4 rel., this is INNO vs scipy) and CasaXPS area within 1e-3, pos 2e-4 eV (DS: after posShift), RMS within 3.5 %.
- `compToCasa(fitted comp).fwhm` vs CasaXPS: TLA ≤ 5e-4, DS ≤ 5e-4, GLN/SGLN n ≤ 2e-4, asymmetric LA / GL T ≤ 1e-3 (empirical δ), ST ≤ 5e-3; `asym` ≤ 1.2e-3 (ST 3e-4).
- `casaToComp('DS(0.1,100)',{be:285.0523,fwhm:0.907537}).comp.be` ≈ 284.9996 (± 2e-4).
- DoF as in §1.
