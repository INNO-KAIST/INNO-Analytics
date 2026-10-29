#!/usr/bin/env python3
"""Apply casa_v4_patch.md (§1–§5) to a COPY of xps_core.js, by exact-text anchors (fails loudly if an anchor moved).
   python3 -I test/ref/casa_v4_apply.py /path/to/copy/xps_core.js            (refuses .../src/xps_core.js unless --live)
   Applies casa_v4_patch.md §1–§6 (Rounds 4–6).
   Then: node test/ref/check_casa_v4.js /path/to/copy/xps_core.js
Written against the 2026-10-29 snapshot (after casa_v3_patch)."""
import sys, os

args = [a for a in sys.argv[1:] if not a.startswith('--')]
if len(args) != 1: sys.exit(__doc__)
P = args[0]
if os.path.abspath(P).endswith(os.path.join('src', 'xps_core.js')) and '--live' not in sys.argv:
    sys.exit('refusing to edit the live src/xps_core.js without --live (another agent may be editing it)')
s = open(P, encoding='utf8').read()


def rep(old, new):
    global s
    n = s.count(old)
    if n != 1: sys.exit(f'anchor found {n}×: {old[:90]!r}')
    s = s.replace(old, new)


# §1 width/asym rule: replace the δ comment, constants and casaWidthAsym0 (refineHalfMax and casaWidthAsym stay)
i = s.index("/* CasaXPS 'fwhm' 숫자와 비대칭 지수 (SPEC_casa_v3 §3")
j = s.index('const CASA_DELTA_LA=0.0025, CASA_DELTA_TAIL=0.011, CASA_DELTA_ST=0.0054, CWA=new Map();')
s = s[:i] + r'''/* CasaXPS 'fwhm' 숫자와 비대칭 지수 (SPEC_casa_v4 §3, CasaXPS 2.3.27 Round 1–4 실측 14 모양): CasaXPS 의 최댓값 x̂ = 성분 원점에 맞춘 격자 k·d (d = F/40,
   DS 는 wL/40) 에서 가장 큰 마디와 이웃 두 마디를 지나는 포물선의 꼭짓점, 반높이 = 꼭짓점 값/2, a < x̂ < b 반높이 점 (a 낮은 BE 쪽),
   l = x̂ − a, h = b − x̂ → fwhm = 2·min(l, h), asym = 1 − l/h (낮은 BE 쪽이 넓으면 음수). 꺾인 꼭대기(LA α≠β n 0, GL/SGL T·PHI)의 'δ' 는
   이 포물선 꼭짓점이 꼬리 쪽으로 밀리는 것(LA: δ ≈ d(β−α)/(2(α+β))). 대칭: 실제 FWHM, 0. ST: 같은 규칙, 위치·폭 미정 (exact:false). */
''' + s[j:]
rep('const CASA_DELTA_LA=0.0025, CASA_DELTA_TAIL=0.011, CASA_DELTA_ST=0.0054, CWA=new Map();', 'const CASA_GRID=1/40, CASA_DS_T=600, CWA=new Map();')
i = s.index('function casaWidthAsym0(c){'); j = s.index('/* INNO 성분 → CasaXPS 표기', i)
s = s[:i] + r'''function casaWidthAsym0(c){
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
/* CasaXPS DS area = amplitude × ∫ P over E0 ± 600·wL (SPEC_casa_v4 §2, α 0 / 0.1 / 0.2, n 0 / 100, two windows: ≤ 4e-4); INNO DS area = area inside [emin, emax]
   → CasaXPS area = INNO area × factor. Raw-DS integrals (closed form); the n-Gaussian changes the window integral by ≤ 5e-5. Cut variant (dsl): null (not read). */
function dsCasaAreaFactor(c,emin,emax){ if(!c||c.shape!=="DS"||c.dsl>0||!(emax>emin)) return null; const sp=shapeParams(c,1), g=sp.wL/2, al=c.alpha, T=CASA_DS_T*sp.wL;
  const N600=dsPrim(T,al,g)-dsPrim(-T,al,g), Nw=dsPrim(emax-c.be,al,g)-dsPrim(emin-c.be,al,g); return Nw>0?N600/Nw:null; }
''' + s[j:]

# §2 compToCasa0: asymmetric LA and GL(p)T(k) widths are now measured rules
rep('if(!sym){ const cw=casaWidthAsym(c); return {notation:not,fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,n,exact:false,note:"asymmetric LA: shape and tail side measured; "+cw.note+(ni>0?" (n > 0 asymmetric LA not read yet)":"")}; }',
    'if(!sym){ const cw=casaWidthAsym(c); return {notation:not,fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,n,exact:Math.abs(n-ni)<0.05,note:"asymmetric LA: shape, tail side and n rule measured (LA(1.53,2.43,0), LA(2.43,1.53,0), LA(2,3,0), LA(2,3,30)); "+cw.note}; }')
rep('fwhm:cw.fwhm,asym:cw.asym,fwhmExact:false,exact:false,note:"same shape as CasaXPS GL(p)T(k) (area/position measured); "+cw.note+" (GL(30)T(2), F 1.2: 1.18051 read)"}; }',
    'fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,exact:true,note:"same shape as CasaXPS GL(p)T(k) (area/position measured); "+cw.note+" (GL(30)T(2), F 1.2: 1.18051 read)"}; }')
# GL/SGL(m,n) with tail: width via the measured rule
rep('fwhm:t?cw.fwhm:s2FWHM(c.shape,sp).fwhm,asym:t?cw.asym:0,fwhmExact:!t,exact:!t&&n.ok,',
    'fwhm:t?cw.fwhm:s2FWHM(c.shape,sp).fwhm,asym:t?cw.asym:0,fwhmExact:true,exact:n.ok,')

# §3 DS: compToCasa(c, o) gets the CasaXPS area when the fit window is known; casaToComp converts the area back
rep('function compToCasa(c){', 'function compToCasa(c,o){\n  o=o||{};')
rep('fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,posShift:sh,exact:false,note:"DS(α,n"+(c.dsl>0?",l":"")+"), Cookbook p. 66: Lorentzian FWHM = wL; shape, tail side and n rule measured (DS(0.1,100)); "+cw.note+"; CasaXPS position = INNO origin + wL·|cot(π/(2−α))| (measured, one reading)"+(c.dsl>0?" — assumed for the cut variant":"")+"; CasaXPS area normalisation not pinned (INNO area = area inside the fit window)"}; }',
    'fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,posShift:sh,areaFactor:ak,area:ak?c.area*ak:null,exact:!!ak&&n.ok,note:"DS(α,n"+(c.dsl>0?",l":"")+"), Cookbook p. 66: Lorentzian FWHM = wL; shape, tail side and n rule measured (α 0, 0.1, 0.2; n 0, 100); "+cw.note+"; CasaXPS position = INNO origin + wL·|cot(π/(2−α))| (measured ≤ 1.6e-4 eV)"+(c.dsl>0?" — assumed for the cut variant":"")+"; CasaXPS area = amplitude × ∫ over ±600·wL"+(ak?" = INNO window area × "+ak.toFixed(5):" — pass the fit window {emin, emax} to convert (INNO area = area inside the fit window)")}; }')
rep('else if(c.shape==="DS"){ const n=nTxt(nOf(sp.wG,sp.wL)), cw=casaWidthAsym(c), sh=c.alpha>0?sp.wL*Math.abs(1/Math.tan(PI/(2-c.alpha))):0;',
    'else if(c.shape==="DS"){ const n=nTxt(nOf(sp.wG,sp.wL)), cw=casaWidthAsym(c), sh=c.alpha>0?sp.wL*Math.abs(1/Math.tan(PI/(2-c.alpha))):0, ak=o.emin!=null&&o.emax!=null?dsCasaAreaFactor(c,+o.emin,+o.emax):null;')
rep('note:"DS(α,n): shape, tail side and n rule (Gaussian FWHM = n·wL/144.87) measured (DS(0.1,100)); CasaXPS area normalisation not pinned — INNO area = area inside the fit window (CasaXPS area ≈ 1.43 × for DS(0.1,100)) — refit the area"',
    'note:"DS(α,n): shape, tail side and n rule (Gaussian FWHM = n·wL/144.87) measured (α 0–0.2, n 0, 100)"')
rep('r.note+="; position: INNO origin = CasaXPS position − wL·|cot(π/(2−α))| (measured, one reading DS(0.1,100))"',
    'r.note+="; position: INNO origin = CasaXPS position − wL·|cot(π/(2−α))| (measured α 0, 0.1, 0.2; ≤ 1.6e-4 eV)"')
rep('  delete r.keepF; return r; }\nfunction compToCasa(c,o){',
    '''  // DS area: CasaXPS Area (amplitude × ∫ over ±600·wL) → INNO area inside the fit window [o.emin, o.emax] (SPEC_casa_v4 §2)
  if(c.shape==="DS"&&!o.innoUnits){ const k=o.emin!=null&&o.emax!=null?dsCasaAreaFactor(c,+o.emin,+o.emax):null;
    if(k){ c.area=area/k; r.areaFactor=k; r.exact=Number.isInteger(+c.casaN)&&!(c.dsl>0); r.note+="; area: INNO window area = CasaXPS area / "+k.toFixed(5)+" (CasaXPS DS area = amplitude × ∫ over ±600·wL)"; }
    else r.note+="; area kept as given — pass the fit window {emin, emax} to convert the CasaXPS DS area (CasaXPS area = amplitude × ∫ over ±600·wL)"; }
  delete r.keepF; return r; }
function compToCasa(c,o){''')

# §4 ST position: CasaXPS pos = INNO origin + largest node k·(base FWHM)/40 of CasaXPS's ST shape  [PINNED, 16 readings ≤ 0.2 meV]
rep('  else r.posShift=mod&&mod.n==="ST"?null:0;',
    '  else if(mod&&mod.n==="ST"&&!o.innoUnits){ const cw=casaWidthAsym(c); if(cw.node!=null){ c.be=be-cw.node; r.posShift=cw.node; r.note+="; position: INNO origin = CasaXPS position − "+cw.node.toFixed(5)+" eV (largest node of CasaXPS\'s ST table, step (base FWHM)/40)"; } else r.posShift=null; }\n  else r.posShift=0;')
rep('fwhm:cw.fwhm,asym:cw.asym,fwhmExact:false,posShift:null,exact:false,',
    'fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,posShift:cw.node,exact:cw.exact&&base.exact,')
rep('; CasaXPS position offset not pinned (ST(1.7)LA(1,1,200): CasaXPS pos − INNO origin = +0.0984 eV); "+base.note}; }',
    '; "+base.note}; }')
rep('CasaXPS position and fwhm number not pinned (ST(1.7)LA(1,1,200): CasaXPS pos − INNO origin = +0.0984 eV) — refit position and width; ',
    'CasaXPS draws it with the step lagged by (base FWHM)/80 and positions it on the largest node of a (base FWHM)/40 grid (SPEC_casa_v4 §4); ')

# §6 ST shape: CasaXPS's step = exclusive rectangle cumulative on the (base FWHM)/40 grid ≈ Φ(x − base FWHM/80)  (SPEC_casa_v4 §4, Rounds 5–6)
#   s2Eval: optional lag p.stLag (eV) on Φ for mod ST + numerical normalisation of the lagged area; shapeParams sets it from c.stCasaK (= base FWHM / F);
#   casaWidthAsym0 always evaluates CasaXPS's ST; casaToComp stores stCasaK on imported ST components.
rep('const muI=mod==="ST"?1/Math.max(1e-9,p.stMu):0, stN=1/(0.5+muI);',
    'const muI=mod==="ST"?1/Math.max(1e-9,p.stMu):0, L=mod==="ST"&&p.stLag>0?p.stLag:0; let stN=1/(0.5+muI);   // L: CasaXPS step lag (base FWHM)/80 (SPEC_casa_v4 §4)')
rep('      Phi[i0]=B.left; for(let i=i0+1;i<M;i++) Phi[i]=Phi[i-1]+cell(gx[i-1],gx[i]); for(let i=i0-1;i>=0;i--) Phi[i]=Phi[i+1]-cell(gx[i],gx[i+1]); } }',
    '      Phi[i0]=B.left; for(let i=i0+1;i<M;i++) Phi[i]=Phi[i-1]+cell(gx[i-1],gx[i]); for(let i=i0-1;i>=0;i--) Phi[i]=Phi[i+1]-cell(gx[i],gx[i+1]); }\n'
    '    if(L){ const P0=Phi.slice(); let dz=0; for(let i=0;i<M;i++){ Phi[i]=B.cdf?B.cdf(gx[i]-L):lag6(P0,k0,h,gx[i]-L); dz+=v[i]*(Phi[i]-P0[i]); } stN=1/(0.5+muI+dz*h); } }')
rep('const ph=B.cdf?B.cdf(x[i]):lag6(Pf,kf,hf,x[i]);', 'const ph=B.cdf?B.cdf(x[i]-L):lag6(Pf,kf,hf,x[i]);')
rep('const md=compMod(c); if(md){ p.mod=md; if(md==="ST"&&c.stN!=null) p.stG=Math.max(0,+c.stN)*CASA_LA_RHO*(c.shape==="V"?p.wL:p.F);',
    'const md=compMod(c); if(md){ p.mod=md; if(md==="ST"&&c.stN!=null) p.stG=Math.max(0,+c.stN)*CASA_LA_RHO*(c.shape==="V"?p.wL:p.F); if(md==="ST"&&c.stCasaK>0) p.stLag=c.stCasaK*(c.shape==="V"?p.wL:p.F)/80;')
rep('c.mod=mod.n; if(mod.n==="ST"){ c.stMu=mod.a[0]; c.stN=mod.a[1]!=null?mod.a[1]:0; c.stG=0; }',
    'c.mod=mod.n; if(mod.n==="ST"){ c.stMu=mod.a[0]; c.stN=mod.a[1]!=null?mod.a[1]:0; c.stG=0; if(!o.innoUnits) c.stCasaK=casaStK(c); }')

# §5 exports
rep('CASA_DELTA_LA,CASA_DELTA_TAIL,CASA_DELTA_ST,', 'CASA_GRID,CASA_DS_T,dsCasaAreaFactor,casaStK,')
open(P, 'w', encoding='utf8').write(s)
print('patched', P)
