/* ============================================================
   INNO Analytics — XPS fit guidance (v1.9.3; v1.9.6 transition-metal 2p3/2 envelope recipes; v1.9.7 Ba · Y · lanthanide recipes, overlap rules): fit check (diagnosis), region recipes, staged fit
   DOM-free: browser main thread, Web Worker (load after xps_core.js) and node.
   Definitions, thresholds, citations and validation: modules/xps/SPEC.md §17 (if they disagree, fix SPEC first).
   General rules only — nothing here is tuned to one data set; the Ce 3d recipe (SPEC §17.3), the Biesinger envelopes (SPEC §17.7) and the Ba/Y/lanthanide recipes (SPEC §17.8) are literature data.
   ============================================================ */
(function(root){
"use strict";
function getCore(){ if(root&&root.XPS_CORE) return root.XPS_CORE; if(typeof require==="function"){ try{ return require("./xps_core.js"); }catch(e){} } return null; }
let _C=null; const CORE=()=>_C||(_C=getCore());

/* ------------------------------------------------------------
   0) Thresholds (SPEC §17.1) — conservative: a finding must be clear in the data, not a matter of taste
   ------------------------------------------------------------ */
const GUIDE_TH={
  boundRel:1e-3,     // parameter within 1e-3·(hi − lo) of a bound (same rule as the stall-rule 'bound' flag, SPEC §3.3b)
  corr:0.95,         // |r| between two component parameters
  corrBg:0.995,      // |r| when one parameter is the in-fit background (c0·k·s) — area↔k is high by construction on wide windows
  corrMax:3,         // report at most this many pairs
  negShare:0.01,     // area < 1 % of the summed component area → negligible
  relSig:1.0,        // σ(area)/area > 100 % → ill-determined
  merge:0.5,         // |ΔBE| < 0.5·mean FWHM between two free main lines → merging
  edgeComp:0.05,     // fitted component sum at a window end > 5 % of its maximum → truncated
  edgeRise:0.10,     // smoothed data rise towards the outside within edgeLen > 10 % of the net peak height → a peak continues outside
  edgeLen:1.0,       // eV
  bgAbove:0.20,      // background above the smoothed data by more than max(1 σ noise, bgAboveAmp · net peak height) on > 20 % of the window
  bgAboveAmp:0.03,
  blobZ:5,           // residual blob |z| ≥ 5 (noise σ and AR(1) ρ estimated from the data's second differences, independent of the model)
  blobHalf:0.6,      // eV, half width of the residual blob window
  blobAmp:0.03,      // and the blob's mean residual ≥ 3 % of the net peak height (a static-background endpoint error of ~1 % is not a missing line)
};
const SEV_ORDER={bad:0,warn:1,info:2};

/* ------------------------------------------------------------
   1) Complex regions (multiplet / satellite-rich) — element·orbital rule, independent of data
   ------------------------------------------------------------ */
const SAT_TM2P=["Cr","Mn","Fe","Co","Ni","Cu"], SAT_RE3D=["Ce","Pr","Nd"];
function complexInfo(el,orb){
  const C=CORE(); if(!el||!orb) return null; const rc=recipesFor(el,orb).map(r=>r.id); const m=C.multipletInfo(el,orb);
  if(m) return {kind:"multiplet",why:`${el} ${orb}: exchange (multiplet) splitting with the open ${m.kind} shell is comparable to or larger than the spin–orbit splitting — the envelope is many lines, not one doublet per chemical state`+(el==="Eu"?"; Eu²⁺ and Eu³⁺ each give their own multiplet (Phys. Rev. B 51, 10146 (1995))":""),closed:el==="Yb"?"Yb²⁺ (4f¹⁴) — a plain doublet (Auto suggest)":(m.closed||null),recipes:rc,note:m.kind==="4f"?RE_NOTE(el,orb,rc):null};
  if(orb==="2p"&&SAT_TM2P.includes(el)) return {kind:"satellite",why:`${el} 2p: multiplet splitting and shake-up satellites of the open 3d shell (closed-shell states such as Cu⁺/Cu⁰ excepted) — one doublet per chemical state does not describe the envelope`,closed:el==="Cu"?"Cu⁺·Cu⁰ (3d¹⁰)":null,recipes:rc};
  if(orb==="3d"&&SAT_RE3D.includes(el)) return {kind:"satellite",why:`${el} 3d: 4f final-state (screening) satellites — several peaks per oxidation state, overlapping across the two spin–orbit halves`,closed:null,recipes:rc,note:rc.length&&recipeById(rc[0]).kind==="emp"?RE_LCF:null};
  if(el==="La"&&(orb==="3d"||orb==="4d")) return {kind:"satellite",why:`La ${orb}: charge-transfer (c4f¹L) satellites — a bonding and an antibonding satellite next to each c4f⁰ spin–orbit line (Sunding et al. 2011); their share differs between La(OH)₃ and La₂O₃`,closed:null,recipes:rc};
  if(orb==="3d"&&RE3D_MORE.includes(el)) return {kind:"satellite",why:`${el} 3d: open-4f final-state multiplet and satellite structure — not one doublet per chemical state`+(el==="Eu"?"; Eu²⁺ and Eu³⁺ each give their own 3d multiplet (mixed valence is a multiplet superposition, not shake-up — Phys. Rev. B 51, 10146 (1995))":el==="Yb"?" (Yb³⁺, 4f¹³; Yb²⁺ is 4f¹⁴)":""),closed:el==="Yb"?"Yb²⁺ (4f¹⁴)":null,recipes:rc,note:RE_NOTE(el,orb,rc)};
  if(rc.length) return {kind:"recipe",why:`${el} ${orb}: literature starting model available`,closed:null,recipes:rc};
  return null;
}
const RE3D_MORE=["Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"];
/* no tabulated parametric source (SPEC §17.8): Sm, Eu, Ho, Er, Tm, Yb → reference-spectrum LCF; Haasch envelopes (Pr, Nd, Gd, Tb, Dy) → also LCF for chemistry */
function RE_NOTE(el,orb,rc){ if(rc.length) return RE_LCF;
  return `No tabulated parametric source is built in for ${el} ${orb} (no recipe). `+RE_LCF+(el==="Eu"?" For Eu²⁺/Eu³⁺ use one reference per valence (e.g. EuO or EuS for Eu²⁺, Eu₂O₃ for Eu³⁺) — the valence ratio is the LCF weight ratio.":"")+(el==="Yb"?" Yb²⁺ (4f¹⁴) 4d — like Lu 4d — is a plain doublet: Auto suggest handles it.":""); }

/* ------------------------------------------------------------
   2) Region recipes (SPEC §17.3) — data, not code. Positions are charge-corrected BE (C 1s C–C/C–H 284.8 eV).
      line = one spin–orbit pair (main = 3d5/2 'v' line, partner = 3d3/2 'u' line at +ΔE, area ratio r, equal FWHM).
   ------------------------------------------------------------ */
const RECIPES=[{
  id:"ce3d", el:"Ce", orb:"3d", name:"Ce 3d (Burroughs v/u)",
  window:[875,925], bg:{type:"active",n:5,c0:null,k:null,s:0,cons:{s:{mode:"free"}}}, shape:"GLS", m:{v:0.3,lo:0,hi:1},
  dE:{v:18.6,lo:18.3,hi:18.9}, dEs:{"Ce3+":{v:18.2,lo:17.9,hi:18.5}}, r:2/3, beTol:0.6, align:3,
  states:{"Ce4+":"Ce⁴⁺","Ce3+":"Ce³⁺"},
  lines:[
    {key:"v0", name:"v0/u0", state:"Ce3+", be:880.5, F:2.4, Flo:1.2, Fhi:3.5},
    {key:"v",  name:"v/u",   state:"Ce4+", be:882.4, F:2.3, Flo:1.2, Fhi:3.5},
    {key:"v1", name:"v′/u′", state:"Ce3+", be:885.0, F:3.1, Flo:1.5, Fhi:4.5},
    {key:"v2", name:"v″/u″", state:"Ce4+", be:888.9, F:5.0, Flo:2.5, Fhi:6.5},
    {key:"v3", name:"v‴/u‴", state:"Ce4+", be:898.2, F:2.3, Flo:1.2, Fhi:3.5}],
  options:[
    {id:"ce4", label:"CeO₂ only (6)", keys:["v","v2","v3"]},
    {id:"ce3", label:"Ce³⁺ only (4)", keys:["v0","v1"]},
    {id:"mix", label:"Mixed (10)", keys:["v0","v","v1","v2","v3"]}],
  derived:[
    {id:"ce3", label:"Ce³⁺", state:"Ce3+", tip:"Ce³⁺ fraction = (v0 + v′ + u0 + u′) / (all ten lines), areas from the fit (each pair = v + u). ± = covariance σ (statistics only). The fraction depends on the background and on the line shapes; beam reduction under X-rays raises it during the measurement (Isaacs et al. 2023; Morgan 2023)."},
    {id:"u3", label:"u‴", key:"v3", part:"partner", tip:"u‴ fraction = u‴ area / (all ten lines). u‴ (≈ 916.7 eV) has no Ce³⁺ counterpart, so it falls as Ce⁴⁺ is reduced; for fully oxidised CeO₂ Shyu et al. (J. Catal. 114, 23, 1988, Fig. 2; nonlinear Shirley) found ≈13–14 % of the Ce 3d region (read from the graph; ±10 % relative), but the value depends on background and fit model, so compare it only between spectra treated the same way."}],
  cite:"Burroughs et al., J. Chem. Soc. Dalton Trans. 1976, 1686 (v/u notation); Romeo et al., Surf. Interface Anal. 20, 508 (1993); Bêche et al., Surf. Interface Anal. 40, 264 (2008); Morgan, Surf. Interface Anal. 55 (2023), doi:10.1002/sia.7254; Isaacs et al., Appl. Surf. Sci. Adv. 18, 100469 (2023); Paparazzo, Mater. Res. Bull. 46, 323 (2011); Shyu et al., J. Catal. 114, 23 (1988) (u‴ fraction)",
  pitfalls:"Five spin–orbit pairs overlap (v‴ ≈ 898 eV lies under u0 ≈ 899 eV, u ≈ 901 eV under u′), so single-line areas are strongly correlated — report the Ce³⁺ fraction, not single peaks. Symmetric lines cannot follow the asymmetric Ce⁴⁺ peaks; the leftover intensity near 885–888 eV is then taken by v′ (Ce³⁺), so 10-line fits tend to overestimate Ce³⁺ — cross-check with the u‴ fraction, which has no Ce³⁺ counterpart. The Ce³⁺ fraction depends on the background (here in-fit Shirley + slope over 875–925 eV; a static Shirley on a sloping baseline biased it by +5 % on synthetic spectra), the u‴ fraction hardly does. X-rays and UHV reduce CeO₂ during acquisition (measure Ce 3d first, low flux). Positions are for C 1s = 284.8 eV and may sit off by charging; u‴ ≈ 916.7 eV is an internal reference. Lines start within ±0.6 eV; ΔE free per state (Ce⁴⁺ 18.6, 18.3–18.9 eV; Ce³⁺ 18.2, 17.9–18.5 eV; Bêche 2008), u/v area = 2/3 (Bêche: I3d5/2/I3d3/2 = 1.5), equal widths within a pair."
}];
/* ------------------------------------------------------------
   2b) Transition-metal 2p3/2 envelope recipes (v1.9.6, SPEC §17.7) — Biesinger et al., Appl. Surf. Sci. 257, 2717 (2011) Tables 1, 3–6;
       Biesinger, Surf. Interface Anal. 49, 1325 (2017) Table 4. Numbers copied from test/ref/biesinger/*.txt (journal tables).
       A species = one or more 'sets'; a set = peaks fitted as ONE envelope (published spacings, area %, FWHM fixed);
       free per set: one position (± beTol around the charge-corrected literature BE) and one area. Optional common FWHM scale.
       Peak rows: first [BE, %, FWHM PE10, FWHM PE20], then [spacing to the previous peak, %, F10, F20] (abs:true → every row starts with its BE).
       Shapes in CasaXPS notation (casaToComp): GL(30) unless stated (Biesinger).
   ------------------------------------------------------------ */
const ENV_COMMON={kind:"env", orb:"2p", bg:{type:"active",n:5,c0:null,k:null,s:0,cons:{s:{mode:"fixed"}}}, beTol:0.25, wScale:[0.6,1.6]};   // Shirley c0 + k·∫P fitted with the envelopes (no slope) — SPEC §17.7: a static Shirley biased fractions by up to 11–12 % where satellite tails reach the window end
const TM_PIT_COMMON="Each species is one envelope: peak spacings, area % and FWHMs are fixed at the published values; only its position (± tolerance around the charge-corrected literature BE, C 1s 284.8 eV) and its total area are free (option: one common FWHM scale for other instruments / pass energies). Report species fractions, not single peaks. Biesinger fits the 2p3/2 only, with a standard Shirley background; here the Shirley (c0 + k·∫P, no slope) is fitted together with the envelopes, because satellite tails reach the window ends (a static Shirley then biased fractions by up to 12 % on synthetic spectra). Differential charging (insulating oxide on a conducting metal) shifts the oxide sets against the metal — then the shared C 1s correction is wrong for one of them and positions hit their window (Fit check). Species not in the list (mixed oxides, carbonates, salts) are taken up by the nearest envelope.";
const RECIPES_ENV=[
 Object.assign({},ENV_COMMON,{id:"ni2p3", el:"Ni", name:"Ni 2p3/2 (Biesinger 2011)", window:[848,869],
  species:[
   {key:"Ni0", label:"Ni metal", short:"Ni", def:true, casa:["LA(1.1,2.2,10)"], peaks:[[852.6,81.2,0.94,0.95],[3.65,6.3,2.70,2.70],[2.38,12.5,2.70,2.70]]},
   {key:"NiO", label:"NiO", short:"NiO", def:true, peaks:[[853.7,14.3,0.98,1.02],[1.71,44.2,2.04,3.25],[5.44,34.0,3.85,3.76],[3.10,3.6,1.97,2.04],[2.38,3.9,2.60,2.44]]},
   {key:"NiOH2", label:"Ni(OH)₂", short:"Ni(OH)₂", def:true, peaks:[[854.9,7.4,1.12,1.16],[0.77,45.3,1.5,2.29],[2.02,3.0,1.59,1.59],[2.79,1.4,1.06,1.06],[1.00,39.2,4.64,4.65],[4.96,3.7,3.08,3.01]]},
   {key:"NiOOH", label:"γ-NiOOH", short:"γ-NiOOH", def:false, peaks:[[854.6,13.8,1.40,null],[0.70,12.4,6.50,null],[0.36,9.7,1.40,null],[0.78,20.7,1.40,null],[1.33,8.7,1.90,null],[3.20,23.3,4.00,null],[3.38,11.4,4.40,null]]}],
  sums:[{id:"Ni2",label:"Ni²⁺",keys:["NiO","NiOH2"],tip:"Ni²⁺ = NiO + Ni(OH)₂ (sum of the two envelopes). The sum is better determined than the split: the two envelopes overlap strongly (main lines 853.7/855.4 vs 854.9/855.7 eV)."}],
  cite:"M.C. Biesinger, B.P. Payne, A.P. Grosvenor, L.W.M. Lau, A.R. Gerson, R.St.C. Smart, Appl. Surf. Sci. 257, 2717 (2011), Table 6 (Ni metal: new line shape LA(1.1,2.2,10))",
  pitfalls:"Ni metal main line is asymmetric (LA(1.1,2.2,10)); its plasmon/satellite lines at +3.65 and +6.0 eV overlap the NiO/Ni(OH)₂ main lines. NiO and Ni(OH)₂ envelopes overlap — their sum (Ni²⁺) is robust, the split is not; γ-NiOOH has a 6.5 eV wide component and only PE 10 widths. Ni metal 2p1/2 starts at ≈869 eV (window end). BEs significant to 0.1 eV, spacings more accurate (Biesinger)."}),
 Object.assign({},ENV_COMMON,{id:"co2p3", el:"Co", name:"Co 2p3/2 (Biesinger 2011)", window:[775,792],
  species:[
   {key:"Co0", label:"Co metal", short:"Co", def:true, peaks:[[778.1,81.0,0.70,0.75],[3.00,11.0,3.00,3.00],[2.00,8.0,3.00,3.00]]},
   {key:"CoO", label:"CoO", short:"CoO", def:false, peaks:[[780.0,46.6,2.23,2.24],[2.10,25.7,2.59,2.66],[3.40,1.6,2.42,2.29],[1.00,26.1,5.28,4.98]]},
   {key:"CoOH2", label:"Co(OH)₂", short:"Co(OH)₂", def:true, peaks:[[780.4,38.1,2.01,2.04],[1.80,26.6,2.60,2.55],[3.79,33.0,4.47,4.47],[4.40,2.4,2.33,2.33]]},
   {key:"CoOOH", label:"CoOOH", short:"CoOOH", def:false, peaks:[[780.1,61.4,null,1.48],[1.32,24.5,null,1.48],[1.68,5.2,null,1.48],[7.07,8.9,null,3.30]]},
   {key:"Co3O4", label:"Co₃O₄", short:"Co₃O₄", def:true, peaks:[[779.6,40.5,1.38,1.39],[1.30,29.1,1.55,1.62],[1.30,15.2,1.94,2.18],[3.00,8.1,4.28,4.44],[4.30,7.2,3.15,3.29]]}],
  sums:[{id:"Co2",label:"Co²⁺",keys:["CoO","CoOH2"],tip:"Co²⁺ = CoO + Co(OH)₂ (sum of the two envelopes; their main lines are 0.4 eV apart, so the split is weakly determined)."}],
  cite:"M.C. Biesinger et al., Appl. Surf. Sci. 257, 2717 (2011), Table 5",
  pitfalls:"Co metal is fitted with GL(30) here (the transcribed table gives no asymmetric shape for Co). CoOOH, Co₃O₄ and CoO/Co(OH)₂ main lines lie within 0.8 eV — the satellites (≈786 eV for Co²⁺, weak for Co³⁺) carry the distinction, so keep them inside the window. CoOOH has PE 20 widths only. Co 2p1/2 starts at ≈793 eV (window end)."}),
 Object.assign({},ENV_COMMON,{id:"fe2p3", el:"Fe", name:"Fe 2p3/2 (Biesinger 2011)", window:[703,721], abs:true, stoich:true,
  species:[
   {key:"Fe0", label:"Fe metal", short:"Fe", def:true, casa:["LA(1.2,4.8,3)"], peaks:[[706.6,100.0,0.88,null]]},
   {key:"FeO", label:"FeO", short:"FeO", def:false, peaks:[[708.4,24.2,1.4,null],[709.7,30.1,1.6,null],[710.9,14.5,1.6,null],[712.1,25.6,2.9,null],[715.4,5.6,2.5,null]]},
   {key:"Fe2O3", label:"Fe₂O₃", short:"Fe₂O₃", def:true, peaks:[[709.8,26.8,1.1,null],[710.8,24.7,1.3,null],[711.6,18.9,1.3,null],[712.7,10.1,1.4,null],[713.7,10.0,2.0,null],[719.3,9.3,2.6,null]]},
   {key:"FeOOH", label:"FeOOH", short:"FeOOH", def:true, peaks:[[710.3,27.0,1.4,null],[711.3,26.5,1.3,null],[712.2,20.6,1.4,null],[713.3,11.3,1.4,null],[714.4,6.3,1.8,null],[719.7,8.3,2.9,null]]},
   {key:"Fe3O4", label:"Fe₃O₄", short:"Fe₃O₄", def:false, sets:[
     {key:"Fe3O4_2", label:"Fe₃O₄ Fe²⁺", peaks:[[708.4,16.6,1.2,null],[709.2,14.8,1.2,null]]},
     {key:"Fe3O4_3", label:"Fe₃O₄ Fe³⁺", peaks:[[710.2,23.7,1.4,null],[711.2,17.8,1.4,null],[712.3,12.2,1.4,null],[713.4,5.7,1.4,null],[714.5,9.1,3.3,null]]}]}],
  sums:[{id:"Fe3",label:"Fe³⁺",sets:["Fe2O3","FeOOH","Fe3O4_3"],tip:"Fe³⁺ = Fe₂O₃ + FeOOH (+ the Fe³⁺ set of Fe₃O₄). Fe₂O₃ and FeOOH envelopes are 0.5 eV apart with the same multiplet pattern, so their split is weakly determined — this sum is not."},
        {id:"Fe2",label:"Fe²⁺",sets:["FeO","Fe3O4_2"],tip:"Fe²⁺ = FeO (+ the Fe²⁺ set of Fe₃O₄)."}],
  cite:"M.C. Biesinger et al., Appl. Surf. Sci. 257, 2717 (2011), Table 4 (printed averages for Fe₂O₃ and FeOOH; Fe metal LA(1.2,4.8,3))",
  pitfalls:"Fe₃O₄ = an Fe²⁺ set (2 peaks) + an Fe³⁺ set (5 peaks); its Fe³⁺ satellite (≈719 eV) is buried under the Fe²⁺ 2p1/2 and is not in the 2p3/2 table — option 1:2 ties the Fe³⁺ set to the Fe²⁺ set (area 2:1, fixed separation 1.8 eV; the table itself gives 68.5 : 31.4). Fe₂O₃ and FeOOH differ by 0.5 eV: report Fe³⁺. Fe metal 2p1/2 (≈719.8 eV) lies inside this window — with metal present the Fit check will show missing intensity there (add a line or end the window at ≈717.5 eV, which cuts the Fe³⁺ satellite). Widths are PE 10 values only. Table typos (α-Fe₂O₃ P4, γ-Fe₂O₃ P3, α-FeOOH P4) are avoided by using the printed averages."}),
 Object.assign({},ENV_COMMON,{id:"mn2p3", el:"Mn", name:"Mn 2p3/2 (Biesinger 2011)", window:[635.5,648.5],
  species:[
   {key:"Mn0", label:"Mn metal", short:"Mn", def:true, peaks:[[638.6,87.0,0.74,0.79],[1.00,13.0,0.74,0.79]]},
   {key:"MnO", label:"MnO", short:"MnO", def:true, peaks:[[640.2,24.0,1.21,1.23],[0.97,27.8,1.21,1.23],[0.93,22.1,1.21,1.23],[0.95,12.5,1.21,1.23],[1.14,4.7,1.21,1.23],[1.75,9.1,1.21,1.23]]},
   {key:"Mn2O3", label:"Mn₂O₃", short:"Mn₂O₃", def:true, peaks:[[640.8,18.9,1.65,1.75],[1.10,44.5,1.65,1.75],[1.27,25.3,1.65,1.75],[1.50,8.5,1.65,1.75],[1.62,3.1,1.65,1.75]]},
   {key:"MnOOH", label:"MnOOH", short:"MnOOH", def:false, peaks:[[641.0,24.0,1.34,1.35],[0.70,24.0,1.34,1.35],[0.81,27.8,1.34,1.35],[1.02,17.5,1.34,1.35],[1.37,6.7,1.34,1.35]]},
   {key:"MnO2", label:"MnO₂", short:"MnO₂", def:true, peaks:[[641.9,41.7,0.84,0.91],[0.86,26.5,0.84,0.91],[0.70,15.5,0.84,0.91],[0.75,9.1,0.84,0.91],[0.85,4.9,0.84,0.91],[1.00,2.5,0.84,0.91]]}],
  sums:[{id:"Mn3",label:"Mn³⁺",keys:["Mn2O3","MnOOH"],tip:"Mn³⁺ = Mn₂O₃ + MnOOH (sum of the two envelopes)."}],
  cite:"M.C. Biesinger et al., Appl. Surf. Sci. 257, 2717 (2011), Table 3 (MnO₂ row, not pyrolusite)",
  pitfalls:"All Mn oxide envelopes span 640–646 eV with spacings ≈ 1 eV and one FWHM per compound — the species are distinguished by small shape differences, so mixtures of three or more oxides are poorly determined (check the σ). MnO has a shake-up at ≈646 eV (in the table, 9.1 %). Mn metal is fitted with GL(30) here (no asymmetric shape in the transcribed table). Mn 2p1/2 starts at ≈649.5 eV (window end)."}),
 Object.assign({},ENV_COMMON,{id:"cr2p3", el:"Cr", name:"Cr 2p3/2 (Biesinger 2011)", window:[571,582.5],
  species:[
   {key:"Cr0", label:"Cr metal", short:"Cr", def:true, peaks:[[574.2,100,0.80,0.90]]},
   {key:"Cr2O3", label:"Cr₂O₃", short:"Cr₂O₃", def:true, peaks:[[575.7,36,0.88,0.94],[1.01,35,0.88,0.94],[0.78,19,0.88,0.94],[1.00,8,0.88,0.94],[0.41,5,0.88,0.94]]},
   {key:"CrOH3", label:"Cr(OH)₃", short:"Cr(OH)₃", def:true, peaks:[[577.3,100,2.58,2.60]]},
   {key:"Cr6", label:"Cr(VI) (CrO₃)", short:"Cr(VI)", def:false, peaks:[[579.6,100,1.28,1.38]]}],
  sums:[{id:"Cr3",label:"Cr³⁺",keys:["Cr2O3","CrOH3"],tip:"Cr³⁺ = Cr₂O₃ + Cr(OH)₃."}],
  cite:"M.C. Biesinger et al., Appl. Surf. Sci. 257, 2717 (2011), Table 1 (Cr(VI): CrO₃ row; the 'mixed species' row is 579.5 eV, FWHM 1.40/1.50)",
  pitfalls:"Cr₂O₃ multiplet (5 peaks, 575.7–578.9 eV) and the broad Cr(OH)₃ line (577.3 eV, FWHM 2.6) overlap — the Cr³⁺ sum is robust. Cr(VI) is reduced by X-rays (measure Cr 2p first). Cr metal is fitted with GL(30) here (no asymmetric shape in the transcribed table). Cr₂O₃ peak % in the table sum to 103 (used as printed — only ratios matter). Cr 2p1/2 starts at ≈583 eV (window end)."}),
 Object.assign({},ENV_COMMON,{id:"cu2p3", el:"Cu", name:"Cu 2p3/2 (Biesinger 2017)", window:[927,948], beTol:0.2,
  species:[
   {key:"Cu0", label:"Cu(0)", short:"Cu(0)", def:true, casa:["GL(90)"], peaks:[[932.63,100,null,0.83]]},
   {key:"Cu2O", label:"Cu₂O", short:"Cu₂O", def:true, casa:["GL(80)"], peaks:[[932.18,100,null,1.00]]},
   {key:"CuO", label:"CuO", short:"CuO", def:true, sat:[2,3,4], shake:1.89, peaks:[[933.11,31,null,2.07],[1.37,33,null,3.05],[6.04,3,null,1.03],[1.13,28,null,3.55],[2.05,6,null,1.17]]},
   {key:"CuOH2", label:"Cu(OH)₂", short:"Cu(OH)₂", def:true, sat:[1,2,3], shake:1.57, peaks:[[934.67,60,null,2.85],[4.63,6,null,2.80],[2.90,28,null,3.66],[1.92,7,null,1.76]]}],
  sums:[{id:"Cu01",label:"Cu(0)+Cu(I)",keys:["Cu0","Cu2O"],tip:"Cu(0) + Cu(I) = Cu(0) + Cu₂O envelopes. Their 2p3/2 lines are 0.45 eV apart with no satellites, so only the sum is reliable from 2p3/2 — separate them with the modified Auger parameter α′ = BE(Cu 2p3/2) + KE(Cu L3M4,5M4,5): Cu(0) 1851.24, Cu₂O 1849.17, CuO 1851.33, Cu(OH)₂ 1850.92 eV (Biesinger 2017 Table 2)."},
        {id:"Cu2",label:"Cu(II)",keys:["CuO","CuOH2"],tip:"Cu(II) = CuO + Cu(OH)₂ envelopes (main lines + shake-up)."},
        {id:"CuShk",label:"Cu(II) shake-up",shake:true,tip:"Cu(II) from the main/shake-up areas (Biesinger 2017 Eq. 2): %Cu(II) = B·(1 + A1s/Bs)/(A + B), %Cu(0)+Cu(I) = (A − A1s·B/Bs)/(A + B) (Eq. 1). A = all main-line area, B = shake-up area (here: the fitted satellite peaks, ≥ 938 eV), A1s/Bs = main/shake-up ratio of pure Cu(II): CuO 1.89, Cu(OH)₂ 1.57 (B-weighted when both are fitted). Independent of how Cu(0) and Cu(I) are split. Note: the Table 4 peak % give main/shake-up = 1.73 (CuO) and 1.46 (Cu(OH)₂), not 1.89/1.57, so on a pure Cu(II) envelope this estimate reads ≈ 6 % (CuO) / 4 % (Cu(OH)₂) above the species fraction."}],
  cite:"M.C. Biesinger, Surf. Interface Anal. 49, 1325 (2017), Table 4 (PE 20 eV, Kratos; Cu(0) GL(90), Cu₂O GL(80)), Table 2 (modified Auger parameters), A1s/Bs ratios, Eq. 1–2",
  pitfalls:"Cu(0) and Cu₂O 2p3/2 lines are 0.45 eV apart (932.63 vs 932.18 eV) with no shake-up: their split from 2p3/2 alone is not reliable — use the modified Auger parameter α′ (Cu(0) 1851.24 vs Cu₂O 1849.17 eV; Table 2) and report Cu(0)+Cu(I). Cu(II) is reduced by X-rays and in UHV (acquire Cu 2p first and fast; watch the shake-up shrink over repeated scans). Differential charging shifts CuO/Cu(OH)₂ against Cu(0). Position tolerance ±0.2 eV (Biesinger 2017 for C 1s-referenced data) keeps Cu(0) and Cu₂O from swapping. Widths are PE 20 values only. Table 4 % give main/shake-up 1.73 (CuO) and 1.46 (Cu(OH)₂) vs the quoted A1s/Bs 1.89 and 1.57."})
];
/* ------------------------------------------------------------
   2c) v1.9.7 Ba · Y · lanthanide recipes (SPEC §17.8) — numbers copied from test/ref/lit2/*.txt (journal tables; test_guide §9 compares them).
       kind "par": parametric spin–orbit pairs with literature links (La 3d/4d Sunding 2011, Y 3d Cole 2020, Ba 4d Karslıoğlu 2018, Ba 3d generic).
       kind "emp": empirical envelope (Haasch 2014, single-crystal REScO₃): published BE/FWHM components, one common position, areas free.
       Species lines: be = charge-corrected start (on C 1s 284.8 after 'c1s' or 'off'); sep:[key, Δ] = position linked to another line of the
       species; aeq: area linked 1:1 (Sunding constraint 4). F/Flo/Fhi: start and free range (Sunding and Cole give no width limits — the
       ranges are fit ranges, not literature values; Cole's FWHMs are the starts of Y).
   ------------------------------------------------------------ */
const PAR_BG={type:"active",n:5,c0:null,k:null,s:0,cons:{s:{mode:"free"}}};   // in-fit Shirley + slope (as the Ce 3d recipe)
const LA_HOW="Sunding et al. (2011) model: per compound a c4f⁰ main line and two charge-transfer c4f¹L satellites (bonding, antibonding; equal areas), each a spin–orbit pair with the 3d (4d) area ratio 2/3 and equal FWHM in the pair; satellite separations fixed at the published values and the same for both spin–orbit components (pairs share ΔE); Voigt ≈ GL sum with one shared Lorentzian fraction. One compound: main and satellite areas free → 'c4f¹L share' (compare Table 5). Two compounds (or 'Fix ratio'): main:satellite ratio fixed at Table 5 and the compounds held at their published relative position → La(OH)₃ fraction. 'Plasmon' (one compound): one plasmon pair at the published separation from the main line, area and width free (constraints 5–6 read as one plasmon per orbital).";
const RE_LCF="Measured reference spectra (Store as reference → LCF, SPEC §12) carry the multiplet and satellite structure; a component fit of this region is an envelope description, not a chemical-state analysis.";
const RECIPES_PAR=[
 {kind:"par", ct:true, id:"la3d", el:"La", orb:"3d", name:"La 3d (Sunding 2011)", window:[826,872], bg:PAR_BG, shape:"GLS", m:{v:0.3,lo:0,hi:1},
  dE:{v:16.8,lo:16.6,hi:17.0}, r:2/3, beTol:0.6, align:3, relShift:true,
  flags:[{id:"fixR",label:"Fix ratio",def:false,tip:"Fix the c4f⁰ : c4f¹L area ratio at Sunding Table 5 (La(OH)₃ 43:57, La₂O₃ 10:90, ±5 %). Always on when both compounds are ticked (otherwise the two compounds cannot be told apart)."},
         {id:"plasmon",label:"Plasmon",def:false,tip:"One plasmon pair per orbital at the published separation from the c4f⁰ line (La(OH)₃ 12.3, La₂O₃ 12.5 eV; Sunding Table 4), area and width free; 2/3 area ratio and shared ΔE as all pairs. One compound only. Under Al Kα its 3d3/2 half (≈ 864 eV) coincides with the La MNN Auger."},
         {id:"auger",label:"MNN Auger",def:true,al:true,tip:"La M45N45N45 Auger as one free line at BE = hν − KE (KE 622.4 La(OH)₃ / 623.1 La₂O₃ eV, Sunding Table 3 → ≈ 864 eV under Al Kα, inside the 3d3/2 satellite region); position ± 1 eV, width free. Only added for Al Kα (other sources move it: the Auger KE does not depend on hν)."}],
  species:[
   {key:"LaOH3", label:"La(OH)₃", short:"La(OH)₃", def:true, ratio:[43,57], pl:12.3, ke:622.4, lines:[
     {key:"m", name:"c4f⁰", be:835.1, F:1.6, Flo:0.6, Fhi:4.5},
     {key:"b", name:"c4f¹L b", sep:["m",3.9], F:2.0, Flo:0.6, Fhi:4.5, role:"sat"},
     {key:"a", name:"c4f¹L ab", sep:["b",2.0], aeq:"b", F:2.0, Flo:0.6, Fhi:4.5, role:"sat"}]},
   {key:"La2O3", label:"La₂O₃", short:"La₂O₃", def:false, ratio:[10,90], pl:12.5, ke:623.1, lines:[
     {key:"m", name:"c4f⁰", be:834.9, F:1.6, Flo:0.6, Fhi:4.5},
     {key:"b", name:"c4f¹L b", sep:["m",4.9], F:2.0, Flo:0.6, Fhi:4.5, role:"sat"},
     {key:"a", name:"c4f¹L ab", sep:["b",3.7], aeq:"b", F:2.0, Flo:0.6, Fhi:4.5, role:"sat"}]}],
  how:LA_HOW,
  cite:"A. Sunding, K. Hadidi, S. Diplas, O.M. Løvvik, T.E. Norby, A.E. Gunnæs, J. Electron Spectrosc. Relat. Phenom. 184, 399–409 (2011), Sect. 2.4.2 (method, six constraints), Tables 3–5",
  pitfalls:"Positions are Sunding's (charge reference tailored to Au, ±0.2 eV), not C 1s 284.8 — the lines are aligned to the data within ±3 eV; ΔE = 851.9 − 835.1 = 16.8 eV (Table 3), free 16.6–17.0. La M45N45N45 Auger: BE ≈ 1486.6 − 622.4 ≈ 864 eV under Al Kα, inside the 3d3/2 satellite region — Sunding measured it (Table 3) and treat the overlap; here 'MNN Auger' adds one free line for it (or use Mg Kα). La 3d3/2 (≈ 852 eV) overlaps Ni 2p3/2 (852.7 eV; LaNiO₃) — Fit check warns when Ni is present. La₂O₃ takes up water and CO₂ in air (La(OH)₃, carbonate): a pure-oxide ratio (10:90) is rarely met on air-exposed samples. Table 5 ratios are ±5 % and depend on Sunding's background and line shapes; the c4f¹L share is the observable to compare. FWHM ranges are fit ranges (the transcribed tables give no widths)."},
 {kind:"par", ct:true, id:"la4d", el:"La", orb:"4d", name:"La 4d (Sunding 2011)", window:[96,117], bg:PAR_BG, shape:"GLS", m:{v:0.3,lo:0,hi:1},
  dE:{v:3.2,lo:3.0,hi:3.4}, r:2/3, beTol:0.6, align:3, relShift:true,
  flags:[{id:"fixR",label:"Fix ratio",def:false,tip:"Fix the c4f⁰ : c4f¹L area ratio at Sunding Table 5 (La 4d: La(OH)₃ 71:29, La₂O₃ 42:58, ±5 %). Always on when both compounds are ticked."},
         {id:"plasmon",label:"Plasmon",def:false,win:[96,129],tip:"One plasmon pair at the published separation from the c4f⁰ line (12.3 / 12.5 eV, Sunding Table 4), area and width free; the window is extended to 129 eV. One compound only."}],
  species:[
   {key:"LaOH3", label:"La(OH)₃", short:"La(OH)₃", def:true, ratio:[71,29], pl:12.3, lines:[
     {key:"m", name:"c4f⁰", be:102.8, F:1.4, Flo:0.5, Fhi:4.0},
     {key:"b", name:"c4f¹L b", sep:["m",3.5], F:1.8, Flo:0.5, Fhi:4.0, role:"sat"},
     {key:"a", name:"c4f¹L ab", sep:["b",1.6], aeq:"b", F:1.8, Flo:0.5, Fhi:4.0, role:"sat"}]},
   {key:"La2O3", label:"La₂O₃", short:"La₂O₃", def:false, ratio:[42,58], pl:12.5, lines:[
     {key:"m", name:"c4f⁰", be:102.6, F:1.4, Flo:0.5, Fhi:4.0},
     {key:"b", name:"c4f¹L b", sep:["m",4.1], F:1.8, Flo:0.5, Fhi:4.0, role:"sat"},
     {key:"a", name:"c4f¹L ab", sep:["b",2.6], aeq:"b", F:1.8, Flo:0.5, Fhi:4.0, role:"sat"}]}],
  how:LA_HOW,
  cite:"A. Sunding et al., J. Electron Spectrosc. Relat. Phenom. 184, 399–409 (2011), Sect. 2.4.2, Tables 3–5",
  pitfalls:"ΔE = 106.0 − 102.8 = 3.2 eV (Table 3), free 3.0–3.4; the 4d5/2 bonding satellite (+3.5 / +4.1 eV) falls on the 4d3/2 main line, so single-line areas are strongly correlated — report the c4f¹L share or the compound fraction. Positions are Sunding's (Au reference), aligned to the data within ±3 eV. FWHM ranges are fit ranges (no widths in the transcribed tables)."},
 {kind:"par", id:"y3d", el:"Y", orb:"3d", name:"Y 3d (Cole 2020)", window:[151,166], bg:PAR_BG, shape:"GLS", m:{v:0.3,lo:0,hi:1},
  dE:{v:2.0,lo:1.9,hi:2.2}, r:0.689, beTol:0.3, align:0, c1s:284.8,
  species:[
   {key:"Y0", label:"Y metal", short:"Y", def:true, lines:[{key:"m", name:"3d", be:156.4, F:0.86, Flo:0.52, Fhi:1.38}]},
   {key:"Y2O3", label:"Y₂O₃", short:"Y₂O₃", def:true, lines:[{key:"m", name:"3d", be:157.5, F:2.38, Flo:1.43, Fhi:3.81}]},
   {key:"YOH3", label:"Y(OH)₃", short:"Y(OH)₃", def:false, lines:[{key:"m", name:"3d", be:157.9, F:1.52, Flo:0.91, Fhi:2.43}]},
   {key:"YCl3", label:"YCl₃", short:"YCl₃", def:false, lines:[{key:"m", name:"3d", be:158.8, F:2.30, Flo:1.38, Fhi:3.68}]}],
  sums:[{id:"Y3",label:"Y³⁺",keys:["Y2O3","YOH3","YCl3"],tip:"Y³⁺ = Y₂O₃ + Y(OH)₃ + YCl₃ (sum of the fitted Y³⁺ species). Better determined than the split: Y₂O₃ and Y(OH)₃ are only 0.4 eV apart (Cole 2020)."}],
  how:"One spin–orbit pair per species (Cole, Kirk, Thorpe 2020): area ratio 3d3/2 : 3d5/2 = 0.689 (their areas, = Avantage SF ratio 2.742/3.978), shared ΔE (2.0–2.1 eV in their tables; free 1.9–2.2), equal FWHM in the pair; start BE and FWHM = their values (C 1s 284.8 eV, as here), position ± 0.3 eV, FWHM 0.6–1.6 × theirs; one shared Lorentzian fraction.",
  cite:"S.L. Cole, K.A. Kirk, R.E. Thorpe, Surf. Sci. Spectra 27, 024010 (2020), doi:10.1116/6.0000475 (Y metal, Y₂O₃); 024007 (2020), doi:10.1116/6.0000317 (Y(OH)₃, YCl₃)",
  pitfalls:"Y₂O₃ (157.5) and Y(OH)₃ (157.9 eV) are 0.4 eV apart with different widths only — their split is weakly determined (check σ; report Y³⁺ = their sum). Y metal is fitted symmetric here (no asymmetric shape in the transcribed data). The ΔE (2.0–2.1 eV) is unresolved at their oxide widths (2.4 eV): keep it inside its range. Positions are for C 1s = 284.8 eV (K-Alpha, Avantage) and are not aligned — a charge-correction error shows as a position at its limit (Fit check)."},
 {kind:"par", id:"ba4d", el:"Ba", orb:"4d", name:"Ba 4d (Karslıoğlu 2018)", window:[84,96.5], bg:PAR_BG, shape:"GLS", m:{v:0.3,lo:0,hi:1},
  dE:{v:2.6,lo:2.5,hi:2.7}, r:2/3, beTol:0.3, align:1.0, off:-0.6, relShift:true,
  species:[
   {key:"BaO", label:"BaO", short:"BaO", def:true, lines:[{key:"m", name:"4d", be:88.8, F:1.6, Flo:0.96, Fhi:2.56}]},
   {key:"BaOH2", label:"Ba(OH)₂", short:"Ba(OH)₂", def:false, lines:[{key:"m", name:"4d", be:90.0, F:1.3, Flo:0.78, Fhi:2.08}]},
   {key:"BaCO3", label:"BaCO₃", short:"BaCO₃", def:true, lines:[{key:"m", name:"4d", be:90.1, F:1.7, Flo:1.02, Fhi:2.72}]},
   {key:"BaNO32", label:"Ba(NO₃)₂", short:"Ba(NO₃)₂", def:false, lines:[{key:"m", name:"4d", be:91.0, F:1.4, Flo:0.84, Fhi:2.24}]}],
  sums:[{id:"BaOHCO",label:"Ba(OH)₂+BaCO₃",keys:["BaOH2","BaCO3"],tip:"Ba(OH)₂ + BaCO₃: their Ba 4d lines are 0.1 eV apart (90.0 / 90.1 eV, Karslıoğlu Table 1) — only the sum is determined; use C 1s (carbonate 289.9 eV) and O 1s (530.8 vs 531.7 eV, same scale) to tell them apart."}],
  how:"One spin–orbit pair per species, ΔE 2.6 eV (all four compounds; free 2.5–2.7), 4d area ratio 2/3, equal FWHM in the pair. Karslıoğlu's BEs are referenced to the Pt Fermi edge, so only their relative shifts are used: the species are held at the published separations (BaO → Ba(OH)₂ +1.2, BaCO₃ +1.3, Ba(NO₃)₂ +2.2 eV) and the whole set is aligned to the data (±1 eV), then one position is free (± 0.3 eV). Start on the C 1s scale = their BE − 0.6 eV (their BaCO₃ carbonate C 1s 289.9 eV vs Schmitz 2001 BaCO₃ 289.1 eV at adventitious 284.6 eV, i.e. 289.3 eV at 284.8) — a one-point estimate; the alignment decides. FWHM start = theirs, free 0.6–1.6 ×.",
  cite:"O. Karslıoğlu, L. Trotochaud, I. Zegkinoglou, H. Bluhm, J. Electron Spectrosc. Relat. Phenom. 225 (2018), doi:10.1016/j.elspec.2018.03.008, Table 1; P.J. Schmitz, Surf. Sci. Spectra 8, 190 (2001) (BaCO₃ C 1s, scale estimate)",
  pitfalls:"BEs are Fermi-edge referenced (in-situ films on Pt, AP-XPS) — absolute values differ from C 1s-referenced data; the alignment (±1 eV) cannot tell a global offset from a change of species (BaO → Ba(OH)₂ is +1.2 eV), so confirm with O 1s (BaO 527.8, Ba(OH)₂ 530.8, BaCO₃ 531.7, Ba(NO₃)₂ 533.6 eV) and C 1s/N 1s on the same scale. Ba(OH)₂ and BaCO₃ are 0.1 eV apart (sum only). Loss features (BaO ≈ 3.2 eV wide, Ba(OH)₂ ≈ 6 eV wide; positions not given) are not modelled. Ba 4d is the Ba line to use when Co is present (Ba 3d5/2 overlaps Co 2p3/2)."},
 {kind:"par", id:"ba3d", el:"Ba", orb:"3d", name:"Ba 3d (generic doublet)", window:[772,801], bg:PAR_BG, shape:"GLS", m:{v:0.3,lo:0,hi:1},
  dE:{v:15.2,lo:15.2,hi:15.2,fix:true}, r:2/3, beTol:0.6, align:3, c1s:284.6,
  species:[{key:"Ba", label:"Ba 3d", short:"Ba", def:true, lines:[{key:"m", name:"3d", be:779.4, F:1.56, Flo:0.8, Fhi:3.5}]}],
  sums:[],
  how:"One generic spin–orbit pair: ΔE fixed at the X-ray Data Booklet value 795.7 − 780.5 = 15.2 eV, area ratio 2/3, equal FWHM; start = Schmitz 2001 BaCO₃ Ba 3d5/2 779.4 eV (FWHM 1.56) at C 1s 284.6 → 779.6 eV at 284.8, aligned to the data within ±3 eV. Surface / lattice Ba components are not built in (no tabulated source).",
  cite:"X-ray Data Booklet (Ba 3d5/2 780.5, 3d3/2 795.7 eV; library BE table); P.J. Schmitz, Surf. Sci. Spectra 8, 190 (2001) (BaCO₃ Ba 3d5/2 779.4 eV, FWHM 1.56)",
  pitfalls:"Ba 3d5/2 (≈ 780 eV) lies on Co 2p3/2 (778.1 eV) and Ba 3d3/2 (795.7) on Co 2p1/2 (793.2 eV): in Ba–Co oxides (e.g. BaCo₀.₄Fe₀.₄Zr₀.₁Y₀.₁O₃, BCFZY) a Ba 3d fit takes up Co intensity and vice versa — use Ba 4d (recipe) for Ba, or fit Co 2p and Ba 3d in one model. The Fe L3VV Auger (≈ 784 eV under Al Kα, library KE 703 eV) also falls here when Fe is present. Two Ba 3d components (surface/lattice) are reported for perovskites, but no tabulated parameters are built in — add a line only if the residual (Fit check) asks for it."}
];
/* Haasch et al. (2014) empirical envelopes — rows [BE, FWHM] as printed (C 1s 285.0 eV) */
const HAASCH={
  pr3d:{el:"Pr",orb:"3d",doi:"10.1116/11.20140906",pg:131,rows:[[928.0,2.37],[932.8,4.70],[949.2,3.24],[953.5,3.56],[957.4,3.57]]},
  pr4d:{el:"Pr",orb:"4d",doi:"10.1116/11.20140906",pg:131,rows:[[109.6,1.00],[112.8,3.25],[116.4,3.66],[119.9,3.00],[122.6,2.63]]},
  nd3d:{el:"Nd",orb:"3d",pg:140,rows:[[977.8,3.92],[982.4,4.34],[1000.4,3.92],[1005.0,4.34]]},
  nd4d:{el:"Nd",orb:"4d",pg:140,rows:[[121.3,6.00],[125.9,5.60]]},
  gd3d:{el:"Gd",orb:"3d",pg:149,rows:[[1186.9,5.34],[1196.6,3.67],[1219.0,5.76]]},
  gd4d:{el:"Gd",orb:"4d",pg:149,rows:[[142.0,3.99],[147.2,3.99]]},
  tb3d:{el:"Tb",orb:"3d",pg:157,rows:[[1243.6,7.58],[1247.1,4.15],[1281.1,6.09]]},
  tb4d:{el:"Tb",orb:"4d",pg:157,rows:[[145.7,2.03],[148.4,3.48],[151.3,3.80],[155.0,3.38]]},
  dy3d:{el:"Dy",orb:"3d",pg:165,rows:[[1298.2,7.66],[1333.4,4.04],[1338.3,6.96]]},
  dy4d:{el:"Dy",orb:"4d",pg:165,rows:[[153.5,3.33],[156.1,1.65],[157.8,2.84],[161.2,3.52],[164.4,1.85],[170.0,3.93],[178.1,5.80],[183.2,6.12]]}};
const HAASCH_CMP={Pr:"PrScO₃",Nd:"NdScO₃",Gd:"GdScO₃",Tb:"TbScO₃",Dy:"DyScO₃"};
Object.keys(HAASCH).forEach(id=>{ const h=HAASCH[id]; const lo=Math.floor(Math.min(...h.rows.map(r=>r[0]-2.5*r[1]))), hi=Math.ceil(Math.max(...h.rows.map(r=>r[0]+2.5*r[1])));   // window = components ± 2.5 FWHM
  const cmp=HAASCH_CMP[h.el];
  RECIPES_PAR.push({kind:"emp", id, el:h.el, orb:h.orb, name:`${h.el} ${h.orb} (Haasch 2014, empirical)`, window:[lo,hi], bg:PAR_BG, casa:"GL(30)", dE:null, beTol:0.3, align:3, c1s:285.0, relShift:true, wScale:[0.6,1.6],
    flags:[{id:"wscale",label:"FWHM×",def:false,tip:"One common FWHM scale (0.6–1.6×) for all components — Haasch's widths are for a Kratos Axis Ultra on single crystals; off: widths fixed at their values."}],
    species:[{key:"env", label:`${cmp} envelope`, short:cmp, def:true, lines:h.rows.map((r,j)=>({key:"c"+(j+1), name:`${h.orb} c${j+1}`, be:r[0], F:r[1], Ffix:true}))}], sums:[],
    how:`Empirical envelope (single crystal ${cmp}(110), areas free): Haasch's ${h.rows.length} components with their BE and FWHM (CasaXPS, Tougaard background; shape not given — GL(30), the CasaXPS default, is used), held at the published spacings; the set is aligned to the data (±3 eV, positions converted from C 1s 285.0 to 284.8 eV), then one common position (± 0.3 eV) and every component area are free. The components are a description of the envelope of one compound, not chemical states — for chemical-state or valence analysis use reference-spectrum LCF.`,
    cite:`R.T. Haasch, E. Breckenfeld, L.W. Martin, Surf. Sci. Spectra 21, ${h.pg} (2014)${h.doi?", doi:"+h.doi:""} (${cmp}(110), Kratos Axis Ultra, mono Al Kα, C 1s 285.0 eV)`,
    pitfalls:`Empirical envelope (single crystal, areas free): no component areas were published, so no area ratio or doublet constraint is imposed and the areas carry no chemical meaning — use the summed area for quantification and measured reference spectra (Store as reference → LCF) for chemical states. Haasch used a Tougaard background; here the in-fit Shirley + slope is the default (change it in the background row if the window ends are not on background).`+(h.el==="Nd"&&h.orb==="3d"?" Nd 3d overlaps the O KLL Auger under Al Kα (O KLL KE ≈ 511 / 492 / 477 eV → BE ≈ 976 / 995 / 1010 eV): with oxygen present part of this envelope is O Auger intensity — use Nd 4d, or Mg Kα.":"")+(h.el==="Pr"&&h.orb==="3d"?" The Pr 3d5/2 satellite (932.8 eV) coincides with Cu 2p3/2 (932.7 eV).":"")}); });
/* spin–orbit Δ per oxidation state (Bêche et al. 2008: ≈18.6 eV for Ce⁴⁺/CeO₂, ≈18.2 eV for Ce³⁺/CePO₄) */
function dEof(R,state){ return (R.dEs&&R.dEs[state])||R.dE; }
RECIPES_ENV.forEach(R=>{ R.options=[]; R.derived=[]; RECIPES.push(R); });
RECIPES_PAR.forEach(R=>{ R.options=[]; R.derived=[]; R.flags=R.flags||[]; RECIPES.push(R); });
/* ---- parametric / empirical recipes: species lines → absolute (charge-corrected, C 1s 284.8 eV) positions and links ---- */
function parDefaults(R){ return R.species.filter(s=>s.def).map(s=>s.key); }
function parFlags(R,o){ const f={}; (R.flags||[]).forEach(q=>{ f[q.id]=o&&o.flags&&o.flags[q.id]!=null?!!o.flags[q.id]:!!q.def; }); return f; }
const isAl=hv=>hv==null||Math.abs(hv-1486.6)<2;
function parLines(R,keys,fl,hv){
  const off=R.c1s!=null?+(284.8-R.c1s).toFixed(4):(R.off||0), multi=keys.length>1, out=[];
  keys.forEach(k=>{ const s=R.species.find(x=>x.key===k); if(!s) return; const be={};
    s.lines.forEach(l=>{ const b=l.sep?be[l.sep[0]]+l.sep[1]:l.be+off; be[l.key]=b;
      out.push({sp:k,key:k+"_"+l.key,lk:l.key,name:(R.species.length>1?(s.short||s.label)+" ":"")+l.name,be:+b.toFixed(4),F:l.F,Flo:l.Flo,Fhi:l.Fhi,Ffix:!!l.Ffix,dbl:!!R.dE,
        sep:l.sep?k+"_"+l.sep[0]:null,sepB:l.sep?l.sep[1]:0,aRef:l.aeq?k+"_"+l.aeq:null,aK:1,role:l.role||"main"}); });
    if(R.ct){ if((fl.fixR||multi)&&s.ratio){ const lb=out.find(x=>x.key===k+"_b"); lb.aRef=k+"_m"; lb.aK=(s.ratio[1]/2)/s.ratio[0]; }
      if(fl.plasmon&&!multi&&s.pl) out.push({sp:k,key:k+"_pl",lk:"pl",name:"plasmon",be:+(be.m+s.pl).toFixed(4),F:3,Flo:1,Fhi:8,Ffix:false,dbl:true,sep:k+"_m",sepB:s.pl,aRef:null,aK:1,role:"pl",rx:true}); } });
  if(R.ct&&fl.auger&&R.flags.some(q=>q.id==="auger")&&isAl(hv)){ const s=R.species.find(x=>x.key===keys[0]); const h=hv||1486.6;
    if(s&&s.ke) out.push({sp:keys[0],key:"aug",lk:"aug",name:R.el+" MNN",be:+(h-s.ke).toFixed(4),F:4,Flo:1.5,Fhi:10,Ffix:false,dbl:false,sep:null,sepB:0,aRef:null,aK:1,role:"aug",rx:true,beTol:1.0}); }
  return out; }
/* Model (SPEC §17.8). o = {shift, be, y, species:[keys], flags:{id:bool}, hv, align (false: no search)} */
function parModel(R,o){
  const C=CORE(); const sh=o.shift||0; let keys=(o.species||[]).filter(k=>R.species.some(s=>s.key===k)); if(!keys.length) keys=parDefaults(R);
  const fl=parFlags(R,o); const lines=parLines(R,keys,fl,o.hv);
  const win=(fl.plasmon&&keys.length===1&&(R.flags.find(q=>q.id==="plasmon")||{}).win)||R.window;
  let lo=win[0]-sh, hi=win[1]-sh; const clipped=[];
  if(o.be&&o.be.length){ const a=o.be[0], b=o.be[o.be.length-1]; if(lo<a){ clipped.push("low"); lo=a; } if(hi>b){ clipped.push("high"); hi=b; } }
  const tag=Math.random().toString(36).slice(2,6);
  const mk=ds=>lines.map(l=>{ const id="r"+R.id+"_"+l.key+"_"+tag; let c;
    if(R.casa){ c=C.casaToComp(R.casa,{be:l.be-sh+ds,fwhm:l.F,area:1,id}).comp; } else c=C.newComp(R.shape,l.be-sh+ds,{F:l.F,m:R.m.v,el:R.el,line:R.orb,id});
    c.name=l.name; c.el=R.el; c.line=R.orb; if(l.dbl) c.doublet={dE:R.dE.v,r:R.r,q:1};
    c.recipe=R.id; c.rk=l.key; c.rsp=l.sp; c.rrole=l.role; if(l.rx) c.rx=1; c.area=1; return c; });
  const rootOf=l=>{ let x=l, k=1, g=0; while(x.aRef&&g++<10){ k*=x.aK; x=lines.find(q=>q.key===x.aRef); } return {key:x.key,k}; };
  let ds=0, score=null, areas=lines.map(l=>rootOf(l).k);
  if(o.be&&o.y){ const reg=C.sliceRegion(o.be,o.y,lo,hi); if(reg&&reg.E.length>=8){
    const E=Float64Array.from(reg.E), y=Float64Array.from(reg.y), n=E.length, ctx={emin:E[0],emax:E[n-1]};
    const B=C.bgStatic(E,y,{type:"smart",n:5}); const z=Float64Array.from(y,(v,i)=>v-B[i]);
    const ys=C.smoothQ(E,y,C.smoothHW(E)); const w=Float64Array.from(ys,v=>1/Math.sqrt(Math.max(v,1e-9)));
    // columns follow the area links (a linked line adds to its root's column), as the fit sees them
    const roots=[...new Set(lines.map(l=>rootOf(l).key))];
    const trial=d=>{ const cs=mk(d); const cols=roots.map(()=>new Float64Array(n));
      cs.forEach((c,i)=>{ const rt=rootOf(lines[i]); const col=cols[roots.indexOf(rt.key)]; const u=C.compEval(c,E,ctx); for(let j=0;j<n;j++) col[j]+=rt.k*u[j]*w[j]; });
      const rhs=Float64Array.from(z,(v,i)=>v*w[i]); const a=C.nnls(cols,rhs); let ss=0; for(let i=0;i<n;i++){ let m=0; for(let j=0;j<cols.length;j++) m+=a[j]*cols[j][i]; const r=rhs[i]-m; ss+=r*r; } return {ss,a}; };
    let best=null; const st=0.05, K=o.align===false?0:Math.round((R.align||0)/st);
    for(let k=-K;k<=K;k++){ const d=k*st, t=trial(d); if(!best||t.ss<best.ss-1e-12*Math.abs(t.ss)) best={d,ss:t.ss,a:t.a}; }
    ds=best.d; score=best.ss; const amax=Math.max(...best.a,0)||1;
    areas=lines.map(l=>{ const rt=rootOf(l); return rt.k*Math.max(best.a[roots.indexOf(rt.key)],0.02*amax,1e-6); }); } }
  const comps=mk(ds); comps.forEach((c,i)=>{ c.area=areas[i]; });
  const byKey={}; lines.forEach((l,i)=>byKey[l.key]=comps[i]);
  const cons={}; comps.forEach(c=>{ cons[c.id]={}; });
  if(R.casa) comps.forEach(c=>{ C.paramNames(c).forEach(k=>{ if(k!=="be"&&k!=="area"&&!C.derivedParam(c,k)) cons[c.id][k]={mode:"fixed"}; }); });
  // positions: sep → linked; position roots (no sep) free ± tolerance, or (relShift/emp) linked to the largest root at the published separation
  const rootsP=lines.map((l,i)=>i).filter(i=>!lines[i].sep), grpA=i=>comps.reduce((s,c,j)=>{ let x=lines[j]; let g=0; while(x.sep&&g++<10) x=lines.find(q=>q.key===x.sep); return x===lines[i]?s+c.area:s; },0);
  const rigid=rootsP.filter(i=>lines[i].role!=="aug"); const ancP=R.relShift&&rigid.length>1?rigid.reduce((a,i)=>a==null||grpA(i)>grpA(a)?i:a,null):null;
  lines.forEach((l,i)=>{ const c=comps[i], q=cons[c.id];
    if(l.sep){ q.be={mode:"link",ref:byKey[l.sep].id,k:1,b:+(c.be-byKey[l.sep].be).toFixed(4)}; return; }
    if(ancP!=null&&i!==ancP&&rigid.includes(i)){ q.be={mode:"link",ref:comps[ancP].id,k:1,b:+(c.be-comps[ancP].be).toFixed(4)}; return; }
    const tol=l.beTol||R.beTol; q.be={mode:"free",lo:+(c.be-tol).toFixed(4),hi:+(c.be+tol).toFixed(4)}; });
  // areas
  lines.forEach((l,i)=>{ const c=comps[i], q=cons[c.id]; if(l.aRef){ q.area={mode:"link",ref:byKey[l.aRef].id,k:+l.aK.toFixed(10),b:0}; c.aref=byKey[l.aRef].id; c.amul=l.aK; } else { q.area={mode:"free"}; c.afloor=0.02; } });
  // widths
  const wscale=!!(fl.wscale&&R.wScale);
  if(wscale){ const anc=comps.reduce((a,c)=>!a||c.area>a.area?c:a,null); cons[anc.id].F={mode:"free",lo:+(anc.F*R.wScale[0]).toFixed(5),hi:+(anc.F*R.wScale[1]).toFixed(5)};
    comps.forEach(c=>{ if(c!==anc) cons[c.id].F={mode:"link",ref:anc.id,k:+(c.F/anc.F).toFixed(8),b:0}; }); }
  else lines.forEach((l,i)=>{ const c=comps[i]; cons[c.id].F=l.Ffix?{mode:"fixed"}:{mode:"free",lo:l.Flo,hi:l.Fhi}; });
  // doublets: shared ΔE (largest pair free in its range, or fixed), r and q fixed; GLS: one shared m
  const dbl=comps.filter(c=>c.doublet); if(dbl.length){ const a=dbl.reduce((x,c)=>!x||c.area>x.area?c:x,null);
    dbl.forEach(c=>{ const q=cons[c.id]; q.r={mode:"fixed"}; q.q={mode:"fixed"}; q.dE=R.dE.fix?{mode:"fixed"}:c===a?{mode:"free",lo:R.dE.lo,hi:R.dE.hi}:{mode:"link",ref:a.id,k:1,b:0}; }); }
  if(!R.casa&&R.shape==="GLS"){ const a=comps.reduce((x,c)=>!x||c.area>x.area?c:x,null); comps.forEach(c=>{ cons[c.id].m=c===a?{mode:"free",lo:R.m.lo,hi:R.m.hi}:{mode:"link",ref:a.id,k:1,b:0}; }); }
  const lab=R.species.filter(s=>keys.includes(s.key)).map(s=>s.short||s.label).join(" + ")+(R.ct&&(fl.fixR||keys.length>1)?" · ratio fixed":"")+(lines.some(l=>l.role==="pl")?" · plasmon":"")+(lines.some(l=>l.role==="aug")?" · Auger":"")+(wscale?" · FWHM×":"");
  return {lo,hi,bg:C.clone(R.bg),comps,cons,recipe:{id:R.id,opt:"par",kind:R.kind,species:keys,flags:fl,label:lab,align:ds,score,clipped,name:R.name}};
}
/* Derived outputs of parametric recipes (SPEC §17.8) → [{id,label,w(c),tip}]; plasmon and Auger lines (c.rx) are outside the total */
function parDerived(R,cs){
  const out=[]; const sp=[...new Set(cs.map(c=>c.rsp))];
  if(R.ct&&sp.length===1){ const s=R.species.find(x=>x.key===sp[0]); const t=s&&s.ratio?s.ratio[1]:null;
    out.push({id:"ct",label:"c4f¹L",w:c=>c.rrole==="sat"?1:0,tip:`c4f¹L share = (bonding + antibonding satellite areas) / (main + satellites) of ${R.el} ${R.orb}, plasmon and Auger lines excluded. Sunding Table 5: ${R.species.map(x=>x.label+" "+x.ratio[1]+" %").join(", ")} (±5 %; their background and shapes)${t!=null?"":""} — a higher share points to the oxide, a lower one to the hydroxide. ± = covariance σ (statistics only).`}); }
  if(sp.length>1) R.species.forEach(s=>{ if(!sp.includes(s.key)) return; out.push({id:"sp_"+s.key,label:s.short||s.label,w:c=>c.rsp===s.key?1:0,
    tip:`${s.label} fraction = its area / all ${R.el} ${R.orb} lines of the recipe`+(R.ct?" (main + satellites; ratios fixed at Sunding Table 5)":"")+(R.kind==="emp"?"":"")+". ± = covariance σ (delta method, statistics only; the literature separations and ratios are held fixed, so model error is not included)."}); });
  (R.sums||[]).forEach(d=>{ const inS=c=>d.keys.includes(c.rsp); if(new Set(cs.filter(inS).map(c=>c.rsp)).size<2||cs.every(inS)) return; out.push({id:d.id,label:d.label,tip:d.tip,w:c=>inS(c)?1:0}); });
  return out; }
/* ---- envelope recipes: data → sets with absolute (charge-corrected) peak positions ---- */
function envF(F10,F20,pe){ const a=pe<=15?F10:F20, b=pe<=15?F20:F10; return a!=null?{F:a,pe:pe<=15?10:20}:{F:b,pe:pe<=15?20:10}; }
function envSpecies(R,s){ return s.sets?s.sets:[{key:s.key,label:s.label,peaks:s.peaks}]; }
function envSets(R,keys,pe){
  const out=[]; R.species.forEach(s=>{ if(!keys.includes(s.key)) return;
    envSpecies(R,s).forEach(st=>{ let be=null; const pk=st.peaks.map((r,j)=>{ be=(R.abs||j===0)?r[0]:+(be+r[0]).toFixed(4); const f=envF(r[2],r[3],pe);
        return {be,pct:r[1],F:f.F,pe:f.pe,casa:(s.casa&&s.casa[j])||"GL(30)",sat:!!(s.sat&&s.sat.includes(j))}; });
      out.push({sp:s.key,key:st.key,label:st.label,spLabel:s.label,shake:s.shake||null,peaks:pk}); }); });
  return out; }
function envDefaults(R){ return R.species.filter(s=>s.def).map(s=>s.key); }
/* Envelope model (SPEC §17.7). o = {shift, be, y, species:[keys], pe (10|20), wscale (bool), stoich (bool, Fe₃O₄ 1:2)} */
function envModel(R,o){
  const C=CORE(); const sh=o.shift||0, pe=o.pe===10||o.pe===20?o.pe:20;
  let keys=(o.species||[]).filter(k=>R.species.some(s=>s.key===k)); if(!keys.length) keys=envDefaults(R);
  const sets=envSets(R,keys,pe);
  let lo=R.window[0]-sh, hi=R.window[1]-sh; const clipped=[];
  if(o.be&&o.be.length){ const a=o.be[0], b=o.be[o.be.length-1]; if(lo<a){ clipped.push("low"); lo=a; } if(hi>b){ clipped.push("high"); hi=b; } }
  const tag=Math.random().toString(36).slice(2,6); const comps=[], pesUsed=new Set();
  sets.forEach(st=>{ const S=st.peaks.reduce((s,p)=>s+p.pct,0); st.jA=st.peaks.reduce((a,p,j)=>p.pct>st.peaks[a].pct?j:a,0); st.S=S; st.comps=[];
    st.peaks.forEach((p,j)=>{ pesUsed.add(p.pe); const id="r"+R.id+"_"+st.key+"_"+j+"_"+tag;
      const cv=C.casaToComp(p.casa,{be:p.be-sh,fwhm:p.F,area:p.pct/S,id}); const c=cv.comp; c.name=st.label+" "+(j+1); c.el=R.el; c.line=R.orb;
      c.recipe=R.id; c.renv=1; c.rsp=st.sp; c.rset=st.key; c.rlab=st.label; c.rk=st.key+"_"+j; c.rpct=p.pct; if(p.sat) c.rsat=1; if(j===st.jA){ c.ranc=1; c.afloor=0.02; }
      st.comps.push(c); comps.push(c); }); });
  // starting areas: NNLS of the set envelopes on the Smart-subtracted data (weights 1/√y), as for the Ce recipe
  let score=null; if(o.be&&o.y){ const reg=C.sliceRegion(o.be,o.y,lo,hi); if(reg&&reg.E.length>=8){
    const E=Float64Array.from(reg.E), y=Float64Array.from(reg.y), n=E.length, ctx={emin:E[0],emax:E[n-1]};
    const B=C.bgStatic(E,y,{type:"smart",n:5}); const ys=C.smoothQ(E,y,C.smoothHW(E)); const w=Float64Array.from(ys,v=>1/Math.sqrt(Math.max(v,1e-9)));
    const cols=sets.map(st=>{ const v=new Float64Array(n); st.comps.forEach(c=>{ const u=C.compEval(c,E,ctx); for(let i=0;i<n;i++) v[i]+=u[i]*w[i]; }); return v; });
    const rhs=Float64Array.from(y,(v,i)=>(v-B[i])*w[i]); const a=C.nnls(cols,rhs); let ss=0; for(let i=0;i<n;i++){ let m=0; for(let j=0;j<cols.length;j++) m+=a[j]*cols[j][i]; ss+=(rhs[i]-m)*(rhs[i]-m); } score=ss;
    const amax=Math.max(...a,0)||1; sets.forEach((st,j)=>{ const A=Math.max(a[j],0.02*amax,1e-6); st.comps.forEach(c=>{ c.area=A*c.area; }); }); } }
  // constraints: per set one position (± beTol) and one area; spacings/area ratios fixed by links; widths fixed or one common scale
  const cons={}; comps.forEach(c=>{ const q=cons[c.id]={}; C.paramNames(c).forEach(k=>{ if(k!=="be"&&k!=="area"&&!C.derivedParam(c,k)) q[k]={mode:"fixed"}; }); });
  sets.forEach(st=>{ const A=st.comps[st.jA];
    cons[A.id].be={mode:"free",lo:+(A.be-R.beTol).toFixed(4),hi:+(A.be+R.beTol).toFixed(4)}; cons[A.id].area={mode:"free"};
    st.comps.forEach((c,j)=>{ c.aref=A.id; c.amul=st.peaks[j].pct/st.peaks[st.jA].pct; if(c===A) return;
      cons[c.id].be={mode:"link",ref:A.id,k:1,b:+(c.be-A.be).toFixed(4)}; cons[c.id].area={mode:"link",ref:A.id,k:c.amul,b:0}; }); });
  const stoich=!!(R.stoich&&o.stoich); const s2=sets.find(s=>s.key==="Fe3O4_2"), s3=sets.find(s=>s.key==="Fe3O4_3");
  if(stoich&&s2&&s3){ const A2=s2.comps[s2.jA], A3=s3.comps[s3.jA]; const kA=2*s2.S/s3.S*(s3.peaks[s3.jA].pct/s2.peaks[s2.jA].pct);   // total Fe³⁺ set = 2 × total Fe²⁺ set
    cons[A3.id].be={mode:"link",ref:A2.id,k:1,b:+(A3.be-A2.be).toFixed(4)}; cons[A3.id].area={mode:"link",ref:A2.id,k:kA,b:0}; A3.area=kA*A2.area;
    s3.comps.forEach((c,j)=>{ if(c!==A3) c.area=A3.area*c.amul; c.aref=A2.id; c.amul=kA*s3.peaks[j].pct/s3.peaks[s3.jA].pct; }); }
  if(o.wscale&&comps.length){ const anc=comps.reduce((a,c)=>!a||c.area>a.area?c:a,null);
    cons[anc.id].F={mode:"free",lo:+(anc.F*R.wScale[0]).toFixed(5),hi:+(anc.F*R.wScale[1]).toFixed(5)};
    comps.forEach(c=>{ if(c!==anc) cons[c.id].F={mode:"link",ref:anc.id,k:+(c.F/anc.F).toFixed(8),b:0}; }); }
  const lab=R.species.filter(s=>keys.includes(s.key)).map(s=>s.short||s.label).join(" + ");
  return {lo,hi,bg:C.clone(R.bg),comps,cons,recipe:{id:R.id,opt:"env",kind:"env",species:keys,pe:[...pesUsed].sort().join("/"),peReq:pe,wscale:!!o.wscale,stoich:stoich&&!!s2,label:lab,align:0,score,clipped,name:R.name,nSets:sets.length}};
}
/* Derived outputs of an envelope recipe: species fractions, sums, Cu(II) shake-up (SPEC §17.7) → [{id,label,w(c),tip}] */
function envDerived(R,cs){
  const out=[]; const present=new Set(cs.map(c=>c.rsp));
  R.species.forEach(s=>{ if(!present.has(s.key)) return; out.push({id:"sp_"+s.key,label:s.short||s.label,w:c=>c.rsp===s.key?1:0,
    tip:`${s.label} fraction = ${s.label} envelope area / all fitted envelopes (${R.el} 2p3/2 only). ± = covariance σ (delta method, statistics only; the envelope shapes are fixed at the literature values, so model error is not included).`}); });
  (R.sums||[]).forEach(d=>{ if(d.shake){ if(!cs.some(c=>c.rsat)) return; const rOf={}; R.species.forEach(s=>{ if(s.shake) rOf[s.key]=s.shake; });
      out.push({id:d.id,label:d.label,tip:d.tip,w:null,shake:rOf}); return; }
    const inS=c=>d.keys?d.keys.includes(c.rsp):d.sets.includes(c.rset); if(cs.filter(inS).length===0) return;
    const nSp=new Set(cs.filter(inS).map(c=>d.keys?c.rsp:c.rset)).size; if(nSp<2) return;   // a sum of one species = that species
    if(cs.every(inS)) return;   // the sum of everything fitted = 100 %
    out.push({id:d.id,label:d.label,tip:d.tip,w:c=>inS(c)?1:0}); });
  return out; }
function recipesFor(el,orb){ return RECIPES.filter(r=>r.el===el&&r.orb===orb); }
function recipeById(id){ return RECIPES.find(r=>r.id===id)||null; }
/* Build a recipe model on raw (uncorrected) BE: shift = charge correction Δ (BE_corr = BE_raw + Δ).
   o = {shift, be, y (full spectrum, raw BE ascending), opt (option id), align (true: search ±recipe.align eV for the best NNLS match, SPEC §17.3)}
   → {lo, hi, bg, comps, cons, recipe:{id,opt,align,...}} */
function recipeModel(id,o){
  const C=CORE(); const R=recipeById(id); if(!R) throw new Error("unknown recipe "+id);
  if(R.kind==="env") return envModel(R,o||{});
  if(R.kind==="par"||R.kind==="emp") return parModel(R,o||{});
  o=o||{}; const sh=o.shift||0; const op=R.options.find(x=>x.id===o.opt)||R.options[R.options.length-1];
  const lines=R.lines.filter(l=>op.keys.includes(l.key));
  let lo=R.window[0]-sh, hi=R.window[1]-sh; const clipped=[];
  if(o.be&&o.be.length){ const a=o.be[0], b=o.be[o.be.length-1]; if(lo<a){ clipped.push("low"); lo=a; } if(hi>b){ clipped.push("high"); hi=b; } }
  const mk=ds=>lines.map((l,i)=>{ const id="r"+R.id+"_"+l.key+"_"+Math.random().toString(36).slice(2,6);
    const c=C.newComp(R.shape,l.be-sh+ds,{F:l.F,m:R.m.v,el:R.el,line:R.orb,id,name:l.name});
    c.doublet={dE:dEof(R,l.state).v,r:R.r,q:1}; c.rk=l.key; c.rstate=l.state; c.recipe=R.id; c.area=1; return c; });
  let ds=0, score=null, areas=null;
  if(o.be&&o.y&&o.align!==false){
    const reg=C.sliceRegion(o.be,o.y,lo,hi); if(reg&&reg.E.length>=8){
      const E=Float64Array.from(reg.E), y=Float64Array.from(reg.y), n=E.length, ctx={emin:E[0],emax:E[n-1]};
      const B=C.bgStatic(E,y,{type:"smart",n:5}); const z=Float64Array.from(y,(v,i)=>v-B[i]);
      const ys=C.smoothQ(E,y,C.smoothHW(E)); const w=Float64Array.from(ys,v=>1/Math.sqrt(Math.max(v,1e-9)));
      const trial=d=>{ const cs=mk(d); const cols=cs.map(c=>Float64Array.from(C.compEval(c,E,ctx),(v,i)=>v*w[i])); const rhs=Float64Array.from(z,(v,i)=>v*w[i]);
        const a=C.nnls(cols,rhs); let ss=0; for(let i=0;i<n;i++){ let m=0; for(let j=0;j<cols.length;j++) m+=a[j]*cols[j][i]; const r=rhs[i]-m; ss+=r*r; } return {ss,a}; };
      let best=null; const st=0.05, K=Math.round((R.align||0)/st);
      for(let k=-K;k<=K;k++){ const d=k*st, t=trial(d); if(!best||t.ss<best.ss) best={d,ss:t.ss,a:t.a}; }
      ds=best.d; score=best.ss; areas=best.a; } }
  const comps=mk(ds); if(areas) comps.forEach((c,i)=>{ c.area=Math.max(areas[i],1e-6); });
  const anc=comps.reduce((a,c)=>!a||c.area>a.area?c:a,null); const cons={};
  const ancS={}; comps.forEach(c=>{ const a=ancS[c.rstate]; if(!a||c.area>a.area) ancS[c.rstate]=c; });
  comps.forEach((c,i)=>{ const l=lines[i]; const q=cons[c.id]={};
    q.be={mode:"free",lo:+(c.be-R.beTol).toFixed(4),hi:+(c.be+R.beTol).toFixed(4)};
    q.F={mode:"free",lo:l.Flo,hi:l.Fhi};
    q.r={mode:"fixed"}; q.q={mode:"fixed"};
    const dS=dEof(R,l.state), aS=ancS[c.rstate]; q.dE=c===aS?{mode:"free",lo:dS.lo,hi:dS.hi}:{mode:"link",ref:aS.id,k:1,b:0};
    q.m=c===anc?{mode:"free",lo:R.m.lo,hi:R.m.hi}:{mode:"link",ref:anc.id,k:1,b:0}; });
  return {lo,hi,bg:C.clone(R.bg),comps,cons,recipe:{id:R.id,opt:op.id,label:op.label,align:ds,score,clipped,name:R.name}};
}
/* Derived recipe output (SPEC §17.3): fractions of summed areas with covariance σ (delta method on the free area parameters).
   fin = finalize() output; returns [{id,label,v,sig,tip}] in % */
function recipeDerived(fin,meta){
  const R=recipeById(meta&&meta.id); if(!R||!fin||!fin.comps) return [];
  const cs=fin.comps.filter(c=>c.recipe===R.id&&!c.rx); if(!cs.length) return [];   // v1.9.7: plasmon / Auger lines (rx) are outside the totals
  const T=cs.reduce((s,c)=>s+Math.max(0,c.area),0); if(!(T>0)) return [];
  const names=(fin.cov&&fin.cov.names)||[], np=(fin.cov&&fin.cov.np)||names.length, Cm=(fin.cov&&fin.cov.C)||[];
  const aIdx={}; names.forEach((q,j)=>{ if(q.k==="area") aIdx[q.id]=j; });
  /* each area → (free root area, factor): follow area links in meta.cons (linear links k·ref + b, derivative k); without cons, the
     build-time links of envelope recipes (c.aref, c.amul) — otherwise the area itself (Ce recipe: every area free) */
  const byId={}; fin.comps.forEach(c=>byId[c.id]=c); const cons=meta.cons||null;
  const root=c=>{ if(cons){ let id=c.id, mul=1, g=0; for(;;){ const q=cons[id]&&cons[id].area; if(q&&q.mode==="link"&&q.ref!=null&&byId[q.ref]&&(q.rk==null||q.rk==="area")&&g++<20){ mul*=q.k!=null?+q.k:1; id=q.ref; } else break; } return {id,mul}; }
    if(c.aref&&byId[c.aref]) return {id:c.aref,mul:c.amul!=null?c.amul:1}; return {id:c.id,mul:1}; };
  const roots=cs.map(root);
  const list=R.kind==="env"?envDerived(R,cs):(R.kind==="par"||R.kind==="emp")?parDerived(R,cs):R.derived;
  return list.map(d=>{
    const part=d.shake?(c=>c.rsat?1+(d.shake[c.rsp]||0):0):d.w?d.w:c=>{ if(d.state) return c.rstate===d.state?1:0; if(d.key) return c.rk===d.key?(d.part==="partner"?(c.doublet?c.doublet.r/(1+c.doublet.r):0):(c.doublet?1/(1+c.doublet.r):1)):0; return 0; };
    const num=cs.reduce((s,c)=>s+part(c)*Math.max(0,c.area),0), f=num/T;
    const g=new Float64Array(np); cs.forEach((c,i)=>{ const rt=roots[i]; const j=aIdx[rt.id]; if(j==null) return; g[j]+=rt.mul*(part(c)-f)/T; });   // 0 + 1·x = x: Ce results unchanged bit for bit
    let s=0; const ok=true; for(let a=0;a<np;a++) if(g[a]) for(let b=0;b<np;b++) if(g[b]){ const v=Cm[a*np+b]; if(isFinite(v)) s+=g[a]*v*g[b]; }   // areas at a bound count as fixed (SPEC §3.4)
    const present=cs.some(c=>part(c)>0);
    return {id:d.id,label:d.label,v:present?100*f:null,sig:present&&ok?100*Math.sqrt(Math.max(0,s)):null,tip:d.tip};
  });
}

/* ------------------------------------------------------------
   3) Noise statistics from the data alone (SPEC §17.1 R8). Second differences d_i = y_i − (y_{i−1} + y_{i+1})/2 remove the smooth
      signal (peaks are many channels wide). Noise model MA(1) — lag-1 correlation ρ, none beyond (measured on VGD flat regions:
      lag 1 ≈ 0.3–0.5, lag ≥ 2 ≈ 0; Poisson counting noise: ρ = 0). Then var(d) = σ²(1.5 − 2ρ), corr(d_i, d_{i+1}) = (1.75ρ − 1)/(1.5 − 2ρ)
      → ρ = (1.5c + 1)/(1.75 + 2c), clipped to [0, 0.5]; the variance of a sum of W noise points ≈ W·σ²(1 + 2ρ).
   ------------------------------------------------------------ */
function d2Corr(rho){ return (1.75*rho-1)/(1.5-2*rho); }
function noiseMA1(y,v){
  const n=y.length; if(n<8) return {rho:0,s2:1,lrv:1};
  const d=[]; for(let i=1;i<n-1;i++){ const sc=v?Math.sqrt(v[i]):1; d.push((y[i]-(y[i-1]+y[i+1])/2)/sc); }
  let a=0,b=0; for(let i=0;i<d.length;i++){ a+=d[i]*d[i]; if(i) b+=d[i]*d[i-1]; } const c=a>0?b/(a*(d.length-1)/d.length):-2/3;
  const rho=Math.max(0,Math.min(0.5,(1.5*c+1)/(1.75+2*c)));
  const q=d.map(x=>x*x).sort((x,y)=>x-y); const med=q[q.length>>1]/0.4549;   // robust variance of d: median of χ²₁ = 0.4549
  const s2=med/(1.5-2*rho); return {rho,s2,lrv:s2*(1+2*rho),c};
}
/* ------------------------------------------------------------
   4) Fit check (SPEC §17.1). ctx = {E, y, fin, el, orb, spec:{be0, be1} (raw spectrum extent), shift, recipe, dEref (literature Δ), anyRecipe}
      → [{id, sev:"bad"|"warn"|"info", label, tip, act|null, ids:[comp ids]}] sorted by severity
   ------------------------------------------------------------ */
function diagnose(ctx,th){
  const C=CORE(); th=Object.assign({},GUIDE_TH,th||{}); const out=[];
  const fin=ctx.fin; if(!fin||!fin.comps) return out;
  const E=ctx.E, y=ctx.y, n=E.length, sh=ctx.shift||0, step=n>1?(E[n-1]-E[0])/(n-1):0.05;
  const comps=fin.comps, byId={}; comps.forEach((c,i)=>byId[c.id]=i);
  const nm=id=>{ if(id==="bg") return "Background"; const c=comps[byId[id]]; return c?(c.name||("#"+(byId[id]+1))):"?"; };
  const lab=k=>(C.PDEF[k]||{lab:k}).lab;
  const fmt=(v,d)=>v==null||!isFinite(v)?"—":(+v).toFixed(d);
  const val=(id,k)=>{ if(id==="bg") return (fin.bg||{})[k]; const c=comps[byId[id]]; if(!c) return null; try{ return C.getP(c,k); }catch(e){ return null; } };
  const free=fin.free||[], cov=fin.cov||{}, np=cov.np||free.length, corr=cov.corr||[], atB=cov.atB||[];
  const summ=fin.summ||[]; const tot=comps.reduce((s,c)=>s+Math.max(0,c.area||0),0)||1;
  const fw=i=>(summ[i]&&summ[i].fwhm)||comps[i].F||1;
  const isFree=(id,k)=>free.some(f=>f.id===id&&f.k===k);
  const push=f=>{ if(!out.some(x=>x.key===f.key)) out.push(f); };
  const dropAct=id=>({type:"drop",id,label:"Remove "+nm(id)});
  const areaDead=new Set();
  // R3a area at 0 (lower bound) → handled as a component finding
  free.forEach((f,j)=>{ if(f.k!=="area"||f.ci<0) return; const p=val(f.id,"area"); if(p!=null&&(atB[j]||p<=Math.max(0,f.lo)+1e-6*tot)) areaDead.add(f.id); });
  // envelope sets that are absent within the statistics (< 5 % of the area, area at 0 or σ > 100 %): their position is not determined by the
  // data, so a position at its window end is a consequence — reported once, as the area finding (R3)
  const envWeak=new Set(); { const tS={}; comps.forEach(c=>{ if(c.renv) tS[c.rset]=(tS[c.rset]||0)+Math.max(0,c.area||0); });
    comps.forEach((c,i)=>{ if(!c.renv||!c.ranc) return; const s=summ[i]&&summ[i].sig?summ[i].sig.area:null; if((tS[c.rset]||0)<0.05*tot||areaDead.has(c.id)||(s!=null&&c.area>0&&s/c.area>th.relSig)) envWeak.add(c.id); }); }
  // R1 parameters at / near a bound
  free.forEach((f,j)=>{
    if(f.k==="area") return; if(f.k==="be"&&envWeak.has(f.id)) return; const p=val(f.id,f.k); if(p==null||!isFinite(p)) return;
    const L=isFinite(f.lo), H=isFinite(f.hi); if(!L&&!H) return;
    const dl=L?p-f.lo:Infinity, dh=H?f.hi-p:Infinity, ref=L&&H?f.hi-f.lo:Math.abs(p)+1;
    if(Math.min(dl,dh)>th.boundRel*ref) return;
    const side=dl<=dh?"lower":"upper", bv=side==="lower"?f.lo:f.hi, add=f.k==="be"?sh:0;
    const shapeP=["m","alpha","beta","T","tau","ET","TM","CT","mu","stMu","stG","tk","tj","ta","tb"].includes(f.k);
    let sev=shapeP||f.id==="bg"?"info":"warn", act=null, hint="";
    if(["F","wG","wL"].includes(f.k)){
      const ci=byId[f.id], others=comps.filter((c,i)=>i!==ci&&c[f.k]!=null&&isFree(c.id,f.k));
      const anc=others.reduce((a,c)=>!a||c.area>a.area?c:a,null);
      const rc=ctx.recipe?ctx.recipe.lines&&ctx.recipe.lines.find(l=>l.key===comps[ci].rk):null;
      if(rc&&f.k==="F"&&side==="upper"&&p>=rc.Fhi-1e-9){ act=null; hint=" The recipe width range is a typical span (not a literature limit) — a wider peak usually means a missing line, a wrong background or a charging/beam-damage broadening, not a wider component."; }
      else if(anc){ act={type:"link",id:f.id,k:f.k,ref:anc.id,label:`Link ${lab(f.k)} to ${nm(anc.id)}`}; hint=` Linking it to ${nm(anc.id)} (same ratio as now) removes one free width; if the line should be broader (e.g. a satellite), fix it at a literature value instead.`; }
      else if(side==="lower"){ act={type:"fix",id:f.id,k:f.k,v:p,label:`Fix ${lab(f.k)}`}; hint=" A width at its lower limit is narrower than the instrument resolution allows — fix it at the width of a well-resolved line in the same spectrum."; }
    } else if(f.k==="be"){
      const edge=Math.abs(bv-E[0])<=1.5*step||Math.abs(bv-E[n-1])<=1.5*step;
      const ce=comps[byId[f.id]];
      if(ce&&ce.renv){ act={type:"widen",id:f.id,k:"be",by:0.2,label:"Widen range ±0.2 eV"}; hint=` The ${ce.rlab} envelope reached the end of its position window (literature BE ± tolerance). Check the charge correction first (differential charging moves oxide sets against the metal); widen only if the reference is uncertain — otherwise this species is absent (its area is then small) or the compound is a different one.`; }
      else if(edge){ act=dropAct(f.id); hint=" The line drifted to the window edge — it is not supported by the data inside the window. Remove it, or widen the window if the line is really outside."; }
      else { act={type:"widen",id:f.id,k:"be",by:0.3,label:"Widen range ±0.3 eV"}; hint=" Widen the range only if the literature allows a larger shift (charge-reference error is ±0.2–0.3 eV); otherwise the model is missing a line or the line belongs to another state."; }
    } else if(f.k==="dE"){
      const ref=ctx.recipe&&ctx.recipe.dE?dEof(ctx.recipe,(comps[byId[f.id]]||{}).rstate).v:ctx.dEref;
      if(ref!=null){ act={type:"fix",id:f.id,k:"dE",v:ref,label:`Fix Δ at ${fmt(ref,2)} eV`}; hint=` Fix the spin–orbit splitting at the literature value (${fmt(ref,2)} eV) unless the doublet is well resolved and isolated.`; }
    } else if(f.k==="r"||f.k==="q"){
      const ref=f.k==="q"?1:null; if(ref!=null) act={type:"fix",id:f.id,k:f.k,v:ref,label:`Fix ${lab(f.k)} at 1`};
      hint=" A doublet ratio at its limit usually compensates for an overlapping line; fix it at the statistical value (2j+1 ratio) and look at the residual.";
    } else if(f.id==="bg"){ hint=" The in-fit background parameter reached its limit: the window ends are probably not on background (widen the window or use a static Shirley/Smart background)."; }
    else hint=" A shape parameter at a physical limit is often fine (e.g. L frac. 0 = pure Gaussian); fix it there to remove a free parameter.";
    if(shapeP&&!act) act={type:"fix",id:f.id,k:f.k,v:p,label:`Fix ${lab(f.k)} at ${fmt(p,3)}`};
    push({key:"bound|"+f.id+"|"+f.k,id:"bound",sev,ids:f.id==="bg"?[]:[f.id],
      label:`${nm(f.id)} ${lab(f.k)} at ${side} bound (${fmt(bv+add,f.k==="be"?2:3)})`,
      tip:`${nm(f.id)} ${lab(f.k)} = ${fmt(p+add,4)} is within ${th.boundRel*100} % of the range of its ${side} bound ${fmt(bv+add,4)}. The value is set by the bound, not by the data, so it has no σ and the other parameters are biased.`+hint,act});
  });
  // R2 correlated pairs
  if(np&&corr.length===np*np){
    const pairs=[]; for(let a=0;a<np;a++) for(let b=a+1;b<np;b++){ if(atB[a]||atB[b]) continue; const r=corr[a*np+b]; if(!isFinite(r)) continue; const fa=free[a], fb=free[b]; if(!fa||!fb) continue;
      const lim=fa.id==="bg"||fb.id==="bg"?th.corrBg:th.corr; if(Math.abs(r)>lim) pairs.push({a,b,r}); }
    pairs.sort((x,y)=>Math.abs(y.r)-Math.abs(x.r));
    pairs.slice(0,th.corrMax).forEach(({a,b,r})=>{ const fa=free[a], fb=free[b]; let act=null, hint="", sevOv=null;
      if(fa.k===fb.k&&["F","wG","wL"].includes(fa.k)&&fa.id!==fb.id){ const ia=byId[fa.id], ib=byId[fb.id]; const big=comps[ia].area>=comps[ib].area?fa:fb, sm=big===fa?fb:fa;
        act={type:"link",id:sm.id,k:fa.k,ref:big.id,label:`Link ${lab(fa.k)}s`}; hint=" Linking the two widths (fixed ratio) usually resolves it."; }
      else if(fa.k==="be"&&fb.k==="be"){ const ia=byId[fa.id], ib=byId[fb.id]; const big=comps[ia].area>=comps[ib].area?fa:fb, sm=big===fa?fb:fa;
        act={type:"offset",id:sm.id,ref:big.id,label:"Fix separation"}; hint=" Fixing the separation at the current (or literature) value leaves one position free."; }
      else if((()=>{ const unres=id=>{ const c=comps[byId[id]]; return !!(c&&c.doublet&&c.doublet.dE<C.DOUBLET_RESOLVE); }; return (fa.k==="dE"&&unres(fa.id))||(fb.k==="dE"&&unres(fb.id))||(fa.id===fb.id&&unres(fa.id)); })()){
        const q=fa.k==="dE"?fa:fb.k==="dE"?fb:fa; const ref=ctx.recipe&&ctx.recipe.dE?dEof(ctx.recipe,(comps[byId[q.id]]||{}).rstate).v:ctx.dEref;
        if(ref!=null&&isFree(q.id,"dE")) act={type:"fix",id:q.id,k:"dE",v:ref,label:`Fix Δ at ${fmt(ref,2)} eV`}; hint=" The doublet is not resolved (Δ < 1.2 eV), so its two halves trade off against position and width; fix Δ at the literature value and keep the area ratio fixed."; sevOv="info"; }
      else if(fa.k==="area"&&fb.k==="area") hint=" The two areas trade off against each other — their sum is well determined, each alone is not. Report the sum, or constrain widths/positions from literature.";
      else hint=" One of the two parameters should be fixed or linked (literature value) — the data cannot separate them.";
      const inRec=fa.id!=="bg"&&fb.id!=="bg"&&comps[byId[fa.id]].recipe&&comps[byId[fa.id]].recipe===comps[byId[fb.id]].recipe;
      if(inRec){ act=null; hint=" Expected in this recipe (overlapping lines): report the derived fraction — its ± includes this correlation — not single-line areas."; }
      push({key:"corr|"+fa.id+fa.k+"|"+fb.id+fb.k,id:"corr",sev:inRec?"info":(sevOv||"warn"),ids:[fa.id,fb.id].filter(x=>x!=="bg"),
        label:`${nm(fa.id)} ${lab(fa.k)} ↔ ${nm(fb.id)} ${lab(fb.k)} correlated (r ${fmt(r,3)})`,
        tip:`Correlation coefficient ${fmt(r,4)} from the covariance matrix (|r| > ${fa.id==="bg"||fb.id==="bg"?th.corrBg:th.corr}). Changing one can be compensated by the other with almost no change in χ², so both σ are large and the split between them is arbitrary.`+hint,act}); });
  }
  // R3 negligible / ill-determined components
  const setTot={}; comps.forEach(c=>{ if(c.renv) setTot[c.rset]=(setTot[c.rset]||0)+Math.max(0,c.area||0); });
  comps.forEach((c,i)=>{ if(c.renv){ if(!c.ranc||!isFree(c.id,"area")) return;   // envelope recipe: one finding per set (its anchor carries the set's free area)
      const share=(setTot[c.rset]||0)/tot, s=summ[i]&&summ[i].sig?summ[i].sig.area:null, rel=s!=null&&c.area>0?s/c.area:null;
      const ids=comps.filter(q=>q.renv&&q.rset===c.rset).map(q=>q.id), act={type:"dropSet",ids,label:"Remove "+c.rlab};
      if(areaDead.has(c.id)) push({key:"dead|"+c.id,id:"dead",sev:"warn",ids,label:`${c.rlab}: area went to zero`,tip:`The fit drove the ${c.rlab} envelope area to 0 — the data do not need this species. Remove it (simpler model, smaller σ for the others) and fit again.`,act});
      else if(rel!=null&&rel>th.relSig) push({key:"ill|"+c.id,id:"ill",sev:"warn",ids,label:`${c.rlab}: area ill-determined (±${fmt(100*rel,0)} %)`,tip:`σ(area)/area = ${fmt(100*rel,0)} % for the ${c.rlab} envelope: the data cannot tell how much of it there is — usually it overlaps another envelope that already explains the intensity, or it is absent. Remove it unless it is expected; the sum of overlapping species (shown with the fractions) stays well determined.`,act});
      else if(share<th.negShare) push({key:"neg|"+c.id,id:"neg",sev:"info",ids,label:`${c.rlab}: negligible (${fmt(100*share,1)} % of area)`,tip:`The ${c.rlab} envelope carries ${fmt(100*share,2)} % of the fitted area (< ${th.negShare*100} %) — the species is absent within the statistics. Remove it if it is not expected.`,act});
      return; }
    const share=Math.max(0,c.area||0)/tot, s=summ[i]&&summ[i].sig?summ[i].sig.area:null, rel=s!=null&&c.area>0?s/c.area:null;
    if(areaDead.has(c.id)) push({key:"dead|"+c.id,id:"dead",sev:"warn",ids:[c.id],label:`${nm(c.id)}: area went to zero`,tip:`The fit drove the area of ${nm(c.id)} to its lower bound 0 — the data do not need this component. Remove it and fit again (a zero-area line makes the covariance singular).`,act:dropAct(c.id)});
    else if(rel!=null&&rel>th.relSig) push({key:"ill|"+c.id,id:"ill",sev:"warn",ids:[c.id],label:`${nm(c.id)}: area ill-determined (±${fmt(100*rel,0)} %)`,tip:`σ(area)/area = ${fmt(100*rel,0)} % (> ${th.relSig*100} %): the data cannot say how large this component is. Typical causes: it overlaps a larger line with free width/position, or it is not there. Constrain it (literature position/width, link the width) or remove it.`,act:dropAct(c.id)});
    else if(share<th.negShare&&comps.length>1) push({key:"neg|"+c.id,id:"neg",sev:"info",ids:[c.id],label:`${nm(c.id)}: negligible (${fmt(100*share,1)} % of area)`,tip:`This component carries ${fmt(100*share,2)} % of the component area (< ${th.negShare*100} %). If it is not required physically, remove it — it adds free parameters without information.`+(c.recipe?" In a recipe model a negligible group means that state is absent; the simpler recipe option may be enough.":""),act:dropAct(c.id)});
  });
  // R4 merging main lines (both positions free and not tied to each other)
  const freeBe=comps.map(c=>isFree(c.id,"be"));
  for(let i=0;i<comps.length;i++) for(let j=i+1;j<comps.length;j++){ if(!freeBe[i]||!freeBe[j]) continue; if(comps[i].renv&&comps[j].renv) continue; if(comps[i].rsp&&comps[i].recipe&&comps[i].recipe===comps[j].recipe&&comps[i].rsp!==comps[j].rsp) continue;   // envelopes / literature species positions overlap by construction (v1.9.7) if(areaDead.has(comps[i].id)||areaDead.has(comps[j].id)) continue;
    const d=Math.abs(comps[i].be-comps[j].be), w=(fw(i)+fw(j))/2; if(!(d<th.merge*w)) continue;
    const sm=comps[i].area<=comps[j].area?comps[i]:comps[j], bg_=sm===comps[i]?comps[j]:comps[i];
    push({key:"merge|"+comps[i].id+"|"+comps[j].id,id:"merge",sev:"warn",ids:[comps[i].id,comps[j].id],label:`${nm(comps[i].id)} and ${nm(comps[j].id)} merge (ΔBE ${fmt(d,2)} eV)`,
      tip:`The two lines are ${fmt(d,3)} eV apart, less than ${th.merge} × their mean FWHM (${fmt(w,2)} eV): they describe one peak. Remove the smaller one (${nm(sm.id)}), or fix its position from literature relative to ${nm(bg_.id)}.`,act:dropAct(sm.id)}); }
  // R5 stop reason
  if(fin.stop==="stalled"||fin.stop==="maxIter"||fin.conv===false){ const fl=((fin.stall&&fin.stall.flags)||[]).map(q=>`${nm(q.id)} ${lab(q.k)} (${q.why.join(", ")})`);
    push({key:"stop",id:"stop",sev:"warn",ids:[],label:fin.stop==="stalled"?"Fit stalled":"Fit did not converge",tip:(fin.stop==="stalled"?"LM stopped early: only heavily damped steps were possible (SPEC §3.3b).":"The iteration limit was reached before χ² settled.")+(fl.length?" Poorly determined: "+fl.join("; ")+".":"")+" With many overlapping lines, 'Staged fit' (heights → positions → all) reaches a minimum more reliably; also fix or link the flagged parameters.",act:comps.length>=3?{type:"staged",label:"Staged fit"}:null}); }
  // R6 window edges
  const P=fin.mu&&fin.B?fin.mu.map((m,i)=>m-fin.B[i]):null; const ys=C.smoothQ(E,y,C.smoothHW(E));
  let net=0; if(fin.B) for(let i=0;i<n;i++) net=Math.max(net,ys[i]-fin.B[i]); let pmax=0; if(P) for(const v of P) pmax=Math.max(pmax,v);
  const ne=Math.max(3,Math.round(0.3/Math.max(step,1e-6))), kL=Math.min(n-1,Math.max(1,Math.round(th.edgeLen/Math.max(step,1e-6))));
  const sideTest=(lowSide)=>{ const idx=lowSide?[...Array(ne).keys()]:[...Array(ne).keys()].map(k=>n-1-k); const pe=P?idx.reduce((s,i)=>s+P[i],0)/ne:0;
    const e0=lowSide?0:n-1, inner=lowSide?[...Array(kL+1).keys()]:[...Array(kL+1).keys()].map(k=>n-1-k); let mn=Infinity; inner.forEach(i=>{ if(ys[i]<mn) mn=ys[i]; });
    const rise=ys[e0]-mn; const byComp=pmax>0&&pe>th.edgeComp*pmax, byData=net>0&&rise>th.edgeRise*net; return {byComp,byData,pe:pmax>0?pe/pmax:0,rise:net>0?rise/net:0}; };
  [["low",true],["high",false]].forEach(([nmS,lowS])=>{ const t=sideTest(lowS); if(!t.byComp&&!t.byData) return;
    const wmax=comps.reduce((m,c,i)=>Math.max(m,fw(i)),1), by=Math.max(2,+(2*wmax).toFixed(1));
    const sp=ctx.spec||{}; const room=lowS?(sp.be0!=null?E[0]-sp.be0:0):(sp.be1!=null?sp.be1-E[n-1]:0);
    const act=room>step*2?{type:"window",side:nmS,by:Math.min(by,room),label:`Widen window (${nmS} BE) by ${fmt(Math.min(by,room),1)} eV`}:null;
    push({key:"edge|"+nmS,id:"edge",sev:"warn",ids:[],label:`Peak cut at the ${nmS}-BE window end`,
      tip:(t.byComp?`The fitted components still have ${fmt(100*t.pe,0)} % of their maximum at the ${nmS}-BE end (> ${th.edgeComp*100} %). `:"")+(t.byData?`The smoothed data rise by ${fmt(100*t.rise,0)} % of the net peak height towards the ${nmS}-BE end within ${th.edgeLen} eV (> ${th.edgeRise*100} %) — signal continues outside. `:"")+"The background endpoint then sits on a peak, which biases every area. Widen the window so both ends are on flat background."+(act?"":" The spectrum has no more data on this side — measure a wider range next time."),act}); });
  // R7 background above the data
  if(fin.B&&ctx.bgType!=="none"){ const nz=CORE().noiseEstimate(E,y); const kap=nz&&nz.kappa>0?nz.kappa:1; let above=0;
    for(let i=0;i<n;i++){ const s=Math.max(Math.sqrt(kap*Math.max(ys[i],1e-9)),th.bgAboveAmp*pmax); if(fin.B[i]>ys[i]+s) above++; }
    const fr=above/n; if(fr>th.bgAbove){ const sw=["shirley","linear","tougaard","spline"].includes(ctx.bgType);
      push({key:"bgabove",id:"bgabove",sev:"warn",ids:[],label:`Background above data on ${fmt(100*fr,0)} % of the window`,tip:`The background is more than max(1 σ noise, ${th.bgAboveAmp*100} % of the peak height) above the smoothed data on ${fmt(100*fr,0)} % of the window (> ${th.bgAbove*100} %): areas there are negative. Usual causes: a window end on a peak or on noise spikes, or a Shirley without a fixed point on a sloping baseline. Use Smart (clipped to the data) or the in-fit Shirley, or move the window ends onto flat background.`,act:sw?{type:"bg",bg:"smart",label:"Use Smart background"}:null}); } }
  // R8 structured residual (blob of the normalised residual; noise from second differences of the data — model-independent)
  if(fin.mu&&fin.v){ const v=fin.v, r=new Float64Array(n); for(let i=0;i<n;i++) r[i]=(y[i]-fin.mu[i])/Math.sqrt(v[i]);
    const nz=noiseMA1(y,v); const w=Math.max(2,Math.round(th.blobHalf/Math.max(step,1e-6)));
    let best=null; for(let i=w;i<n-w;i++){ let s=0, a=0; for(let j=i-w;j<=i+w;j++){ s+=r[j]; a+=y[j]-fin.mu[j]; } const z=s/Math.sqrt(nz.lrv*(2*w+1)), amp=a/(2*w+1);
      if(pmax>0&&Math.abs(amp)<th.blobAmp*pmax) continue; if(!best||Math.abs(z)>Math.abs(best.z)) best={i,z,amp}; }
    if(ctx.debug) ctx.debug.blob={z:best&&best.z,be:best&&E[best.i],rho:nz.rho,s2:nz.s2};
    if(best&&Math.abs(best.z)>=th.blobZ){ const be=E[best.i], pos=best.z>0;
      let act=null; if(pos){ const anc=comps.reduce((a,c)=>!a||c.area>a.area?c:a,null); act={type:"add",be,F:Math.min(3,Math.max(0.8,anc?fw(byId[anc.id]):1.2)),h:Math.max(0,y[best.i]-fin.mu[best.i]),label:`Add a line at ${fmt(be+sh,1)} eV`}; }
      push({key:"resid",id:"resid",sev:"warn",ids:[],label:pos?`Missing intensity near ${fmt(be+sh,1)} eV`:`Model above data near ${fmt(be+sh,1)} eV`,
        tip:`The residual summed over ±${th.blobHalf} eV around ${fmt(be+sh,2)} eV is ${fmt(Math.abs(best.z),1)} σ ${pos?"above":"below"} zero and ${fmt(100*Math.abs(best.amp)/pmax,1)} % of the peak height (thresholds ${th.blobZ} σ and ${th.blobAmp*100} %; noise level and its neighbour correlation ρ = ${fmt(nz.rho,2)} (MA(1)) are estimated from the data's second differences, so the test does not depend on the model). `+(pos?"A line or a satellite is missing there (or a width/position is wrong). Add one and compare χ²/BIC — or use the region recipe if this is a satellite-rich region.":"The model is too high there: check the line shapes (asymmetry/tails) and the background.") ,act}); } }
  // R9 complex region → recipe or reference spectra
  const ci_=complexInfo(ctx.el,ctx.orb);
  if(ci_&&!(ctx.recipe&&ctx.recipe.id)){ const rc=ci_.recipes[0];
    push({key:"complex",id:"complex",sev:"info",ids:[],label:rc?`${ctx.el} ${ctx.orb}: use the region recipe`:`${ctx.el} ${ctx.orb}: multiplet/satellite region`,
      tip:ci_.why+". "+(rc?(recipeById(rc)&&(recipeById(rc).kind==="par"||recipeById(rc).kind==="emp")?recipeById(rc).how+" "+(recipeById(rc).kind==="emp"?"":"Derived values are shown after the fit."):recipeById(rc)&&recipeById(rc).kind==="env"?"The recipe fits each chemical species as one literature envelope (Biesinger: fixed spacings, area % and FWHMs; one position ± tolerance and one area per species) on the 2p3/2 window with a Shirley background; species fractions are shown after the fit. Pick the species in the recipe row.":"The recipe sets literature lines with doublet links (ΔE, area ratio, equal widths), the window and the background; derived output (e.g. Ce³⁺ fraction) is shown after the fit."):"No literature recipe is built in for this region: fit only the envelope and use the summed area, or compare with measured reference spectra (Store as reference → LCF), which carry the multiplet structure.")+(ci_.closed?` Closed-shell exception: ${ci_.closed}.`:"")+(ci_.note?" "+ci_.note:""),
      act:rc?{type:"recipe",recipe:rc,label:"Apply recipe"}:null}); }
  // R10 peak overlaps with other elements present (v1.9.7, SPEC §17.8)
  if(ctx.present) overlapFindings(ctx.el,ctx.orb,ctx.present,{hv:ctx.hv,lo:E[0]+sh,hi:E[n-1]+sh}).forEach(push);
  out.sort((a,b)=>SEV_ORDER[a.sev]-SEV_ORDER[b.sev]);
  return out;
}
/* ------------------------------------------------------------
   4b) Peak overlaps between elements (v1.9.7, SPEC §17.8) — a short generic table, energies from the library (X-ray Data Booklet BE;
       Auger KE from the library's Avantage table) or the transcribed sources. A rule fires on region (el, orb) when the partner element is
       present (survey / other regions of the sample) and the partner energy lies in the fit window; Auger rules only for Al Kα.
       side = {el, orb (region), sub, be} or {el, auger, ke:[...]} (BE = hν − KE)
   ------------------------------------------------------------ */
const BKL="X-ray Data Booklet (library BE table)", AVK="library Auger KE (Avantage table)";
const OVERLAPS=[
 {id:"BaCo", sev:"warn", a:{el:"Ba",orb:"3d",sub:"3d5/2",be:780.5}, b:{el:"Co",orb:"2p",sub:"2p3/2",be:778.1}, src:BKL,
  tip:"Ba 3d5/2 (780.5 eV) and Co 2p3/2 (778.1 eV) overlap, and so do Ba 3d3/2 (795.7) and Co 2p1/2 (793.2 eV) — in Ba–Co perovskites (e.g. BCFZY, BaCo₀.₄Fe₀.₄Zr₀.₁Y₀.₁O₃) each fit takes up the other element's intensity, which biases both areas, chemical-state splits and at%. Use Ba 4d for Ba (recipe 'Ba 4d', Karslıoğlu 2018), or fit Co 2p and Ba 3d together in one model (Ba doublet ΔE 15.2 eV, Co 2p envelopes)."},
 {id:"LaNi", sev:"warn", a:{el:"La",orb:"3d",sub:"3d3/2",be:853.0}, b:{el:"Ni",orb:"2p",sub:"2p3/2",be:852.7}, src:BKL,
  tip:"La 3d3/2 (853.0 eV) lies on Ni 2p3/2 (852.7 eV) — in La–Ni oxides (LaNiO₃, La₂NiO₄) Ni 2p3/2 cannot be fitted alone: fit La 3d (recipe, which fixes the 3d3/2 to the 3d5/2 by ΔE and area ratio) together with Ni 2p, or use Ni 3p for Ni."},
 {id:"LaMNN", sev:"info", al:true, a:{el:"La",orb:"3d",sub:"3d3/2 satellites"}, b:{el:"La",auger:"M45N45N45",ke:[622.4,623.1]}, src:"Sunding et al. 2011, Table 3 (La(OH)₃ / La₂O₃)",
  tip:"La M45N45N45 Auger (KE 622.4 / 623.1 eV, Sunding 2011 Table 3) sits at BE ≈ 1486.6 − 622.4 ≈ 864 eV under Al Kα, inside the La 3d3/2 satellite region. Sunding et al. treat this overlap in their La 3d analysis; the La 3d recipe adds one free line for it ('MNN Auger'). With Mg Kα the Auger moves out of the La 3d region (its KE does not depend on hν)."},
 {id:"NdOKLL", sev:"warn", al:true, a:{el:"Nd",orb:"3d",sub:"3d",be:980.4}, b:{el:"O",auger:"KLL",ke:[511,492,477]}, src:AVK+"; Haasch et al. 2014 (NdScO₃ note)",
  tip:"O KLL Auger (library KE 511 / 492 / 477 eV → BE ≈ 976 / 995 / 1010 eV under Al Kα) overlaps Nd 3d5/2 (980.4 eV) and 3d3/2 (1003.3 eV) — Haasch et al. (2014) note this for NdScO₃. With oxygen present part of the Nd 3d envelope is O Auger intensity: use Nd 4d, or Mg Kα."},
 {id:"PrCu", sev:"warn", a:{el:"Pr",orb:"3d",sub:"3d5/2",be:928.8}, b:{el:"Cu",orb:"2p",sub:"2p3/2",be:932.7}, src:BKL+"; Haasch et al. 2014 (Pr 3d5/2 component 932.8 eV)",
  tip:"Cu 2p3/2 (932.7 eV) coincides with the Pr 3d5/2 satellite (Haasch 2014 component at 932.8 eV; Pr 3d5/2 main 928.8 eV, Booklet) — with both present neither region can be fitted alone."},
 {id:"FeBa", sev:"info", al:true, a:{el:"Ba",orb:"3d",sub:"3d5/2",be:780.5}, b:{el:"Fe",auger:"L3VV",ke:[703]}, src:AVK,
  tip:"Fe L3VV Auger (library KE 703 eV → BE ≈ 784 eV under Al Kα) lies on the high-BE side of Ba 3d5/2 (780.5 eV) — relevant for Ba–Fe oxides (BCFZY): it adds intensity a Ba doublet cannot take up. Mg Kα moves it."},
 {id:"CoFe", sev:"info", al:true, a:{el:"Fe",orb:"2p",sub:"2p3/2",be:706.8}, b:{el:"Co",auger:"L3VV",ke:[774]}, src:AVK,
  tip:"Co L3VV Auger (library KE 774 eV → BE ≈ 713 eV under Al Kα) falls inside the Fe 2p3/2 region (706.8–720 eV) — with Co present the Fe 2p3/2 envelope fit takes up Auger intensity. Mg Kα moves it."},
 {id:"NiMn", sev:"info", al:true, a:{el:"Mn",orb:"2p",sub:"2p3/2",be:638.7}, b:{el:"Ni",auger:"L3VV",ke:[846]}, src:AVK,
  tip:"Ni L3VV Auger (library KE 846 eV → BE ≈ 641 eV under Al Kα) lies on Mn 2p3/2 (638.7 eV) — with Ni present the Mn 2p3/2 fit takes up Auger intensity. Mg Kα moves it."},
 {id:"RuC", sev:"warn", a:{el:"Ru",orb:"3d",sub:"3d3/2",be:284.2}, b:{el:"C",orb:"1s",sub:"1s",be:284.2}, src:BKL,
  tip:"Ru 3d3/2 (284.2 eV) coincides with C 1s (284.2 eV, Booklet) — C 1s charge referencing and the Ru 3d3/2 area are both affected; fit Ru 3d and C 1s in one model (Ru 3d5/2 at 280.0 eV is free of C) or use Ru 3p."},
 {id:"BiS", sev:"warn", a:{el:"Bi",orb:"4f",sub:"4f5/2",be:162.3}, b:{el:"S",orb:"2p",sub:"2p3/2",be:162.5}, src:BKL,
  tip:"Bi 4f5/2 (162.3 eV) lies on S 2p (2p3/2 162.5 eV, Booklet) — in Bi sulfides fit Bi 4f and S 2p together, or use S 2s for S."},
 {id:"SrP", sev:"info", a:{el:"Sr",orb:"3d",sub:"3d5/2",be:134.2}, b:{el:"P",orb:"2p",sub:"2p",be:135.0}, src:BKL,
  tip:"Sr 3d5/2 (134.2 eV) and P 2p (135 eV, Booklet) overlap — in Sr phosphates fit them together or use P 2s."}];
function overlapFindings(el,orb,present,o){
  o=o||{}; const out=[]; if(!el||!orb) return out; const P=present instanceof Set?present:new Set(present||[]); const hv=o.hv, al=isAl(hv), h=hv||1486.6;
  const inWin=be=>o.lo==null||o.hi==null||(be>=o.lo-1&&be<=o.hi+1);
  const fmt=v=>(+v).toFixed(1);
  OVERLAPS.forEach(R=>{ if(R.al&&!al) return;
    [[R.a,R.b],[R.b,R.a]].forEach(([me,ot])=>{ if(me.auger||me.el!==el||me.orb!==orb) return; if(ot.el!==el&&!P.has(ot.el)) return;
      let be, what; if(ot.auger){ const bes=ot.ke.map(k=>h-k).filter(inWin); if(!bes.length) return; be=bes; what=`${ot.el} ${ot.auger} Auger`; }
      else { if(!inWin(ot.be)) return; be=[ot.be]; what=`${ot.el} ${ot.sub}`; }
      out.push({key:"ovl|"+R.id,id:"overlap",sev:R.sev,ids:[],rule:R.id,label:`${what} (${be.map(fmt).join("/")} eV) overlaps ${el} ${me.sub}`,tip:R.tip+(ot.auger&&hv==null?" (Al Kα assumed.)":"")+" Source: "+R.src+".",act:null}); }); });
  return out; }
/* Diagnosis straight from a problem + fitOnce result (tests and Worker) */
function diagnoseFit(P,r,o){ const C=CORE(); o=o||{}; const fin=o.fin||C.finalize(P,r); return diagnose(Object.assign({E:P.E,y:P.y,fin,bgType:P.bg&&P.bg.type},o),o.th); }

/* ------------------------------------------------------------
   5) One-click actions on a model {lo, hi, bg, comps, cons} (raw BE) — the UI calls this after pushHist (SPEC §17.2)
   ------------------------------------------------------------ */
function applyAction(model,act,o){
  const C=CORE(); o=o||{}; const cs=model.comps, byId={}; cs.forEach(c=>byId[c.id]=c); model.cons=model.cons||{};
  const q=id=>(model.cons[id]=model.cons[id]||{});
  if(act.type==="dropSet"){ const gone=new Set(act.ids||[]); if(!cs.some(c=>gone.has(c.id))) return null; relinkWithout(model,gone); return "Removed "+(act.label||"").replace(/^Remove /,""); }
  if(act.type==="drop"){ const P={comps:cs,cons:model.cons}; C.removeComp(P,act.id); model.comps=P.comps; model.cons=P.cons; return "Removed "+(act.label||"").replace(/^Remove /,""); }
  if(act.type==="link"){ const c=byId[act.id], a=byId[act.ref]; if(!c||!a) return null; const ka=C.getP(a,act.k), kc=C.getP(c,act.k); q(act.id)[act.k]={mode:"link",ref:act.ref,k:ka>0?+(kc/ka).toFixed(6):1,b:0}; return "Linked"; }
  if(act.type==="fix"){ const c=byId[act.id]; if(!c) return null; C.setP(c,act.k,act.v); q(act.id)[act.k]={mode:"fixed"}; return "Fixed"; }
  if(act.type==="widen"){ const c=byId[act.id]; if(!c) return null; const cs_=q(act.id)[act.k]||{mode:"free"}; const [dlo,dhi]=C.defBounds(c,act.k); const lo=cs_.lo!=null?cs_.lo:dlo, hi=cs_.hi!=null?cs_.hi:dhi; q(act.id)[act.k]={mode:"free",lo:+(lo-act.by).toFixed(4),hi:+(hi+act.by).toFixed(4)}; return "Widened"; }
  if(act.type==="offset"){ const c=byId[act.id], a=byId[act.ref]; if(!c||!a) return null; q(act.id).be={mode:"link",ref:act.ref,k:1,b:+(c.be-a.be).toFixed(4)}; return "Separation fixed"; }
  if(act.type==="window"){ const be0=o.be0!=null?o.be0:-Infinity, be1=o.be1!=null?o.be1:Infinity; if(act.side==="low") model.lo=Math.max(be0,model.lo-act.by); else model.hi=Math.min(be1,model.hi+act.by); return "Window widened"; }
  if(act.type==="bg"){ model.bg=Object.assign({},model.bg,{type:act.bg}); if(act.bg!=="active"){ delete model.bg.conv; } return "Background changed"; }
  if(act.type==="add"){ const anc=cs.reduce((a,c)=>!a||c.area>a.area?c:a,null); const c=C.newComp(anc?anc.shape:"GLS",act.be,{F:act.F,el:anc?anc.el:"",line:anc?anc.line:""});
    if(anc&&anc.m!=null) c.m=anc.m; const hh=C.trueFWHM(c.shape,C.shapeParams(c,1)).ymax; c.area=Math.max(1e-6,(act.h||1)/(hh||1));
    q(c.id).be={mode:"free",lo:+(act.be-1).toFixed(4),hi:+(act.be+1).toFixed(4)}; q(c.id).F={mode:"free",lo:0.3,hi:Math.max(4,act.F*2)}; if(anc&&anc.m!=null) q(c.id).m={mode:"link",ref:anc.id,k:1,b:0};
    cs.push(c); return "Added"; }
  return null;
}

/* Remove a group of components; links that pointed into the group are re-anchored on the first remaining user with the same linear
   relation (values are current, so k = v_user/v_anchor for area/widths, b = Δ for positions) — exact, unlike a k = 1 re-anchor. */
function relinkWithout(model,gone){
  const C=CORE(); const cons=model.cons=model.cons||{};
  model.comps=model.comps.filter(c=>!gone.has(c.id)); gone.forEach(id=>delete cons[id]);
  const byRef={}; model.comps.forEach(c=>{ const q=cons[c.id]||{}; for(const k in q){ const l=q[k]; if(l&&l.mode==="link"&&gone.has(l.ref)) (byRef[l.ref+"|"+k]=byRef[l.ref+"|"+k]||[]).push(c); } });
  for(const key in byRef){ const k=key.split("|")[1], users=byRef[key], Y=users[0];
    const v=c=>C.getP(c,k), vY=v(Y);
    if(k==="be"){ const R0=Y.recipe&&recipeById(Y.recipe), tol=R0&&R0.beTol?R0.beTol:0.3; cons[Y.id][k]={mode:"free",lo:+(vY-tol).toFixed(4),hi:+(vY+tol).toFixed(4)}; }
    else cons[Y.id][k]={mode:"free"};
    users.forEach(u=>{ if(u===Y) return; cons[u.id][k]=k==="be"?{mode:"link",ref:Y.id,k:1,b:+(v(u)-vY).toFixed(4)}:{mode:"link",ref:Y.id,k:vY>0?+(v(u)/vY).toFixed(8):1,b:0}; });
    if(k==="F"&&Y.renv){ const R0=recipeById(Y.recipe); if(R0&&R0.wScale) cons[Y.id].F={mode:"free",lo:+(Y.F*R0.wScale[0]).toFixed(5),hi:+(Y.F*R0.wScale[1]).toFixed(5)}; }
    if(k==="area"&&Y.renv){ Y.ranc=1; Y.afloor=0.02; } }
  model.comps.forEach(c=>{ if(c.aref&&gone.has(c.aref)){ const q=cons[c.id]&&cons[c.id].area; c.aref=q&&q.mode==="link"?q.ref:c.id; c.amul=q&&q.mode==="link"?q.k:1; } });
}

/* ------------------------------------------------------------
   6) Staged fit (SPEC §17.4): (1) heights only → (2) positions within their ranges, widths linked → (3) the user's constraints
      (normal fitOnce, stall rule included). The normal fit from the same start is run too; the lower χ²_d (data-variance χ², the
      common objective of both) is kept, so the result is never worse than the normal fit.
   ------------------------------------------------------------ */
const WKEYS=["F","wG","wL"];
function cpParams(dst,src){ dst.comps=CORE().clone(src.comps); ["c0","k","s","_khi"].forEach(k=>{ if(src.bg&&src.bg[k]!=null) dst.bg[k]=src.bg[k]; }); }
function stageCons(P,user,stage){
  const C=CORE(); const cons=C.clone(user||{}); const get=(id,k)=>(cons[id]||{})[k];
  if(stage===1){ P.comps.forEach(c=>{ const q=cons[c.id]=cons[c.id]||{}; C.paramNames(c).forEach(k=>{ if(k==="area") return; q[k]={mode:"fixed"}; }); }); return cons; }
  // stage 2: positions as the user set them; widths linked to the anchor's (current ratios); all other shape/doublet parameters fixed
  P.comps.forEach(c=>{ const q=cons[c.id]=cons[c.id]||{}; C.paramNames(c).forEach(k=>{ if(k==="area"||k==="be"||WKEYS.includes(k)) return; q[k]={mode:"fixed"}; }); });
  WKEYS.forEach(k=>{ const has=P.comps.filter(c=>C.paramNames(c).includes(k)&&!C.derivedParam(c,k)&&(get(c.id,k)||{}).mode!=="fixed");
    if(has.length<2) return; const anc=has.reduce((a,c)=>!a||c.area>a.area?c:a,null); const ua=get(anc.id,k);
    cons[anc.id][k]=ua&&ua.mode==="free"?ua:{mode:"free"};
    has.forEach(c=>{ if(c===anc) return; const va=C.getP(anc,k), vc=C.getP(c,k); cons[c.id][k]={mode:"link",ref:anc.id,k:va>0?vc/va:1,b:0}; }); });
  return cons;
}
function stagedFit(prob,opt){
  const C=CORE(); opt=opt||{}; const t0=Date.now();
  const base=()=>{ const P=C.clone(C.stripProb(prob)); P.E=prob.E; P.y=prob.y; return P; };
  const P0=base(); { const m=C.buildMap(P0); C.applyLinks(P0,m); }
  const user=C.clone(prob.cons||{}); const stages=[];
  // normal path from the same start
  const PN=base(); const rN=C.fitOnce(PN,opt); const sN=C.fitStats(PN,rN);
  // stage 1: heights (areas) only
  const P1=base(); cpParams(P1,P0); P1.cons=stageCons(P1,user,1); const r1=C.fitOnce(P1,opt); stages.push({stage:1,chi2:r1.chi2,iter:r1.iter+(r1.iter1||0),stop:r1.stop});
  // stage 2: positions within their ranges, widths linked
  const P2=base(); cpParams(P2,P1); P2.cons=stageCons(P2,user,2); const r2=C.fitOnce(P2,Object.assign({},opt,{initAreas:false})); stages.push({stage:2,chi2:r2.chi2,iter:r2.iter+(r2.iter1||0),stop:r2.stop});
  // stage 3: the user's constraints, normal fit path (stall rule applies)
  const P3=base(); cpParams(P3,P2); P3.cons=user; const r3=C.fitOnce(P3,Object.assign({},opt,{initAreas:false})); const s3=C.fitStats(P3,r3); stages.push({stage:3,chi2:r3.chi2,iter:r3.iter+(r3.iter1||0),stop:r3.stop});
  const useStaged=s3.chi2d<sN.chi2d*(1-1e-9);
  const P=useStaged?P3:PN, r=useStaged?r3:rN;
  return {P,r,staged:{used:useStaged?"staged":"normal",chi2dStaged:s3.chi2d,chi2dNormal:sN.chi2d,chi2Staged:r3.chi2,chi2Normal:rN.chi2,stopStaged:r3.stop,stopNormal:rN.stop,stages,ms:Date.now()-t0}};
}
/* Worker job: staged fit → finalize (same summary as finalJob) + staged report */
function stagedJob(prob,opt){ const C=CORE(); const s=stagedFit(prob,opt); const out=C.finalize(s.P,s.r); out.staged=s.staged; return out; }

const XPS_GUIDE={GUIDE_TH,RECIPES,RECIPES_ENV,RECIPES_PAR,HAASCH,OVERLAPS,overlapFindings,parModel,parLines,parDefaults,parFlags,envSets,envDefaults,envModel,relinkWithout,complexInfo,recipesFor,recipeById,recipeModel,recipeDerived,noiseMA1,d2Corr,diagnose,diagnoseFit,applyAction,stageCons,stagedFit,stagedJob};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_GUIDE; else root.XPS_GUIDE=XPS_GUIDE;
})(typeof self!=="undefined"?self:this);
