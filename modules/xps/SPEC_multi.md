# XPS series · multivariate · ARXPS · imaging — specification (v1.0, 2026-10-28)

Companion to `SPEC.md` (single spectra). This file defines `src/xps_multi.js`, `src/xps_arxps.js`, `src/xps_image.js`
and the VAMAS / multi-column additions to `src/xps_io.js`. The independent Python references (`test/ref/ref_multi.py`, `ref_arxps.py`,
`ref_image.py`) were written from this file with numpy/scipy. If code and this file disagree, fix this file first, then the code.

Only published equations are used (no Avantage / CasaXPS code or data). Every equation names its source. Where a source could not be
checked verbatim, the section says so and states what was implemented instead.

**Notation.** A data matrix **D** is r energies × c spectra, stored as an array of c columns (one per etch level, angle or pixel).
"Counts per unit" `cfac` = dwell × scans (1 when the data are already counts). Poisson variance of a CPS value y is `max(y, 1/cfac)/cfac`.

---

## 1. Series processing (depth profiles, angle series, line scans) — `xps_multi.js`

### 1.1 Assembly
`seriesFromSpectra(spectra)` groups spectra by (sample, region). The coordinate is VAMAS experimental variable 0 (`meta.vamas.expNum[0]`,
e.g. *Etch Time* s, *Etch Level*, *Angle* deg); with no experimental variable but varying take-off angles, the analyser take-off polar
angle; otherwise the block index. `commonGrid` puts all spectra of a series on one energy axis (identical axes are kept; otherwise linear
interpolation onto the overlap at the finest step).

### 1.2 Region integration per level
Background B from `XPS_CORE.bgStatic` (SPEC §2: linear / Shirley / smart / Tougaard, endpoint means of n points), area
`A = ∫(y − B) dE` (trapezoid). Uncertainty: A is written as a linear function of the data, `A ≈ Σ aᵢ yᵢ`, with
`aᵢ = wᵢ − (1/n)[i ∈ low end]·Σₖ wₖ(1 − tₖ) − (1/n)[i ∈ high end]·Σₖ wₖ tₖ`, `tₖ = (Eₖ − E_lo)/(E_hi − E_lo)`, w = trapezoid weights.
This is exact for a linear background and is used as the approximation for the iterative backgrounds. `σ_A² = Σ aᵢ² max(yᵢ, 1/cfac)/cfac`.
*Validated* against Monte Carlo (400 Poisson realizations, peak on a Shirley step): mean σ / MC sd = 1.04 (linear), 0.99 (Shirley).
The Shirley background's own dependence on the peak is not in the linearization; for this case the effect is within the MC precision (≈ 4 %),
but it is an approximation, not an identity.

### 1.3 Peak models per level
`modelSeries(E, Y, fitFn)` calls a user/core fitter `fitFn(E, y, j) → {components:[{name, area, sArea}]}` for every level, so the
`XPS_CORE` fitter (or a fixed template) can be applied to every level without this module knowing the fitter.

### 1.4 Composition vs level
`compositionSeries(rows, areas, opts)` applies `XPS_CORE.quantify` (SPEC §4: `n = A/(S·T·λ)`, at% with delta-method σ, n.d. rule
A ≤ 0 or A < 2σ) at every level, with the same transmission and IMFP options as single spectra.

### 1.5 Information-depth ("escape-depth") correction of a sputter profile
**What the term means here.** Two different corrections carry this name. (a) The λ(KE) factor of the sensitivity
(`n = A/(S·T·λ)`, SPEC §4; the energy dependence of λ is the escape-depth correction of Briggs & Seah, *Practical Surface Analysis* Vol. 1,
via Seah & Dench 1979) — already part of quantification and applied per level by §1.4. (b) The deconvolution of the exponential
information-depth response from a sputter depth profile, implemented here. I could not consult the Briggs & Seah text to quote an equation
number for (b), so the equation is derived:

Hofmann's MRI depth-resolution model (S. Hofmann, Surf. Interface Anal. 21 (1994) 673) contains the information-depth term
`g_λ(z) = (1/λₑ) exp(−(z − z₀)/λₑ)`, z ≥ z₀, `λₑ = L cosθ`. With mixing and roughness neglected the measured signal is
`I(z₀) = ∫_{z₀}^{∞} c(z) exp(−(z − z₀)/λₑ)/λₑ dz`. Differentiating, `dI/dz₀ = (I − c)/λₑ`, hence the exact inverse

  **c(z₀) = I(z₀) − λₑ · dI/dz₀**

Implementation `infoDepthCorrect(z, X, λₑ, {hw})`: X = normalized intensities n_k (not at%) of each element vs depth z (nm, from
`depthScale(t, rate)`), derivative by a local quadratic least-squares fit over 2hw + 1 points (fewer at the ends), corrected n clipped at 0,
then renormalized to at%. Limits: assumes constant λₑ per element and no sputter mixing/roughness; amplifies noise (choose hw);
depth scale must be known.

---

## 2. Multivariate analysis — `xps_multi.js`

### 2.1 SVD and preprocessing
Thin SVD by one-sided Jacobi (Hestenes) on the side with fewer vectors (accurate small singular values); for very large matrices
(image stacks) the r × r Gram matrix with cyclic Jacobi eigen-decomposition (`svdGram`). Sign convention: the element of largest |value| of
each left vector is positive. Preprocessing (`preprocess`): scale `none`; `area` (each spectrum ÷ its sum); `poisson` —
**D̃ = G^{-½} D H^{-½}**, G = diag(mean spectrum), H = diag(mean of each spectrum) (Keenan & Kotula, Surf. Interface Anal. 36 (2004) 203);
centre `none` or `mean` (row means). Reconstruction undoes centring and scaling.

### 2.2 Number of factors: eigenvalue statistics (Malinowski)
λⱼ = sⱼ², R = max(r, c′), C = min(r, c′), c′ = c (c − 1 when mean-centred).
- Real error `RE(n) = √( Σ_{j>n} λⱼ / (R (C − n)) )`, indicator `IND(n) = RE(n)/(C − n)²`, n_IND = argmin
  (Malinowski, Anal. Chem. 49 (1977) 612; *Factor Analysis in Chemistry*, 3rd ed., Wiley 2002).
- F-test: `F(1, C − n) = [Σ_{j=n+1}^{C} (R − j + 1)(C − j + 1)] / [(R − n + 1)(C − n + 1)] · λₙ / Σ_{j=n+1}^{C} λⱼ`, significance level
  = P(F_{1,C−n} > F) (regularized incomplete beta); n_F = count of leading factors significant at α = 5 % before the first non-significant one
  (Malinowski, J. Chemometrics 3 (1988) 49).
These assume homoscedastic noise. For Poisson data with a large intensity range they under-count (validation §2.7).

### 2.3 Number of factors: noise limit (recommended)
The Avantage-style "fit within the noise" idea, made statistically exact for Poisson data with a weighted low-rank fit
(maximum-likelihood PCA, Wentzell, Andrews, Hamilton, Faber & Kowalski, J. Chemometrics 11 (1997) 339):
`min Σᵢⱼ wᵢⱼ (Dᵢⱼ − (A Bᵀ)ᵢⱼ)²`, `wᵢⱼ = cfac / max(D̂ᵢⱼ, 1/cfac)` (Pearson weights from the current model, refreshed each sweep),
alternating weighted least squares (B column-wise, then A row-wise), start = rank-n reconstruction of the Poisson-scaled SVD, stop when
|Δχ²| ≤ 1e-10 χ² (max 300 sweeps). `χ²ₙ` is the converged weighted residual, `dofₙ = (r − n)(c − n)`;
**n_noise = smallest n with χ²ₙ ≤ dofₙ + 2√(2 dofₙ)**. n = 0 uses the data as variance.
Method `pca` (unweighted reconstruction, Pearson χ²) is kept for speed but is biased high when the dynamic range is large: unweighted
factors leak high-variance noise into low-count channels (observed χ²/dof = 1.13–1.20 at the true rank for 30 000-count peaks on a
100-count background).

### 2.4 Target testing (TFA)
Non-centred PCA only. Target x (length r) row-scaled like the data; `t = S⁻¹Uᵀx`, `x̂ = U Uᵀ x` (n factors);
`AET = √(Σ(x − x̂)²/r)`, `REP = RE(n)·‖t‖`, `RET = √max(AET² − REP², 0)`, **SPOIL = RET/REP**; < 3 acceptable, 3–6 moderately
acceptable, > 6 rejected (Malinowski, Anal. Chim. Acta 103 (1978) 339). Background note: a background that is identical at every level is
nearly a combination of the concentration profiles and is not resolved as a factor; then true component spectra fail the test. Subtract the
background first (`subtractBackground`) or use a background that differs between levels.

### 2.5 Pure-component recovery
- **OPA** (Cuesta Sánchez, Toft, van den Bogaert & Massart, Anal. Chem. 68 (1996) 79): picks the n most dissimilar spectra
  (max det of the Gram matrix of unit-normalized spectra; first reference = mean spectrum).
- **SIMPLISMA** purest energy channels (Windig & Guilment, Anal. Chem. 63 (1991) 1425): μᵢ, σᵢ over spectra, α = f·max μ (f = 3 %),
  `λᵢ = √(σᵢ² + (μᵢ + α)²)`, `yᵢⱼ = Dᵢⱼ/(λᵢ√c)`, weight `wᵢ = det(YᵀY)` of the selected rows plus row i (1 for the first),
  purity `pᵢ = wᵢ σᵢ/(μᵢ + α)`.
- **ITTFA** (Gemperline, J. Chem. Inf. Comput. Sci. 24 (1984) 206; Vandeginste et al., Anal. Chim. Acta 173 (1985) 253): needle in the level
  domain at the OPA levels, `x ← V_n V_nᵀ x`, negatives → 0, ÷ max, until max change < 1e-10; spectra by row-wise NNLS.
- **MCR-ALS** with non-negativity in both modes (Tauler, Chemom. Intell. Lab. Syst. 30 (1995) 133): alternate `C = NNLS(S, D)`,
  `S = NNLS(C, Dᵀ)`, each spectrum normalized to unit area, optional closure; stop when the lack of fit `LOF = 100·√(Σe²/Σd²)` changes
  < 1e-10 relative. Default start: SIMPLISMA profiles. **Limit: rotational ambiguity** — many (S, C) pairs give the same LOF; OPA-started
  solutions of the 3-state test had 11.6 at.% fraction error at the same LOF as the 1.7 at.% SIMPLISMA solution. Subtract the background first.
- NNLS everywhere: Lawson & Hanson (*Solving Least Squares Problems*, 1974, ch. 23) active-set in Gram form (`nnlsGram`), with optional
  sign-free variables kept in the passive set.

### 2.6 Fitting levels to reference spectra (Avantage "fit to reference levels" idea)
`y ≈ Σₖ aₖ refₖ(E − δ) + Σ_{m<bg} b_m P_m(E)`, aₖ ≥ 0, b sign-free, P₀ = 1, P₁ = (E − Ē)/(half range).
Weights `w = cfac/max(y, 1/cfac)` (pass 1), then `cfac/max(μ, 1/cfac)` (pass 2, IRLS). `σ(a) = √diag((AᵀWA)_PP⁻¹)` over the passive set
(non-zero a plus all b); `φ = χ²/(N − p)` reported (multiply σ by √φ when φ > 1). Optional common shift δ: grid of step `shiftStep`
over ±maxShift, then parabolic refinement of χ²(δ) (linear interpolation of references, constant extrapolation). Amount fractions use the
reference areas; their σ ignores coefficient covariance (stated).

### 2.7 Validation (test/test_multi.js; numbers from the full run)
Reference agreement (numpy/scipy): SVD σ and vectors ≤ 6e-15 · σ₁; PCA in 6 preprocessing combinations, RE/IND/F/SL ≤ 7e-14 relative,
n_IND and n_F identical; MLPCA χ² table ≤ 1e-14; target-test statistics ≤ 1e-8; OPA, SIMPLISMA, ITTFA identical; NNLS vs scipy ≤ 1e-15;
reference fit (12 cases, bvls) ≤ 3e-14; integration ≤ 1e-15; composition ≤ 1e-16; information-depth correction ≤ 2e-14.
Synthetic truth: see the run log for the current numbers — at release:
- factor count (K = 2, 3, 4 states + level-dependent background, 10 Poisson realizations each): noise criterion correct 10/10, 10/10, 10/10
  at 3 000-count peaks and 9/10, 8/10, 10/10 at 30 000; at 300 counts the weakest factor is below the noise edge (σ/edge 0.7–1.8, not
  required). F-test (5 %) fails at 30 000 counts (heteroscedastic noise); IND 3/10 for K = 4 at 30 000.
- target test: true state spectra SPOIL 0.0–0.4, flat background 0.4; 1 eV-shifted or half-width spectra 53 and 26 (rejected).
- MCR-ALS on background-subtracted data: max fraction error 2.4 (K = 2), 1.7 (K = 3), 8.6 at.% (K = 4, four overlapping O 1s states).
- reference fit with true references (900 coefficients): ±1.96σ coverage 96.0 %, mean z 0.07, sd z 0.98; energy shifts 0.23 / −0.31 / 0.05 eV
  recovered within 0.003 eV.
- information-depth correction: noise-free sigmoid interface, max error 0.51 → 0.065; with 1 % noise the mean of the max at% error
  63.5 → 11.5 at.% (hw = 3).
- synthetic SDP VAMAS file (SiO₂/Si, 12 levels): composition from Shirley areas within 0.4 at.% of the generator.

---

## 3. Readers — `xps_io.js` additions
- **VAMAS (ISO 14976)** — every spectrum now carries `meta.vamas` = {mode (NORM, SDP, MAP …), scan, nReg, nPos, nX, nY, expVars[{label, unit}],
  expVals (strings), expNum (numbers), blockIndex, x, y (MAP/MAPDP), fov, linescan (MAPSV/MAPSVDP/SEM), takeOffPolar, takeOffAz, srcPolar, srcAz,
  tilt, tiltAz, rotation, sputter {energy, current, widthX, widthY, polar, azimuth, mode} (SDP-type modes), sputIon, add (additional numerical
  parameters), tag (display text)}. `readFile` returns `vamas` (header summary) for VAMAS files. Spectrum names are unchanged
  (CasaXPS writes a *Data Set* experimental variable into every file).
- **Fix:** the six first/last-linescan values are read only for MAPSV, MAPSVDP and SEM (ISO 14976; the earlier code read them for MAP,
  which would misalign every MAP block).
- **Multi-column text** `textMultiParse(text, name, {axis, hv})`: energy column + many spectra; the last text line before the numbers gives
  column names; the first number in each name is the coordinate (e.g. "Etch Time 30 s" → 30). Energy BE ascending; KE axes converted.
  `textParse` (single spectra) is unchanged.
- Synthetic test files (not measured data): `test/ref/make_vamas.py` → `test/data/synthetic/` (SDP depth profile with sputter fields,
  NORM angle series with *Angle* variable and per-block take-off angles, MAP 8 × 6 with field of view, multi-column text) + `*_truth.json`.

---

## 4. Angle-resolved XPS — `xps_arxps.js`
Angles θ are emission angles from the surface normal, lengths in nm. Straight-line attenuation (elastic scattering neglected); results
above ≈ 60–65° are increasingly biased (warning issued for θ > 65°).

### 4.1 Attenuation lengths
- Seah & Dench, Surf. Interface Anal. 1 (1979) 2: elements `λ_m = 538/E² + 0.41(aE)^½` monolayers, inorganic compounds
  `2170/E² + 0.72(aE)^½`, λ = a·λ_m nm; organics `λ_d = 49/E² + 0.11E^½` mg m⁻², λ = 1000 λ_d/ρ nm (ρ in kg m⁻³). `a³ = M/(ρ N_A n)`.
- From an IMFP (e.g. `XPS_CORE.tpp2m`): `L = λ_in (1 − 0.735 ω)`, ω = single-scattering albedo (Jablonski & Powell; NIST SRD 82
  user guide, eq. E.3). ω must be supplied.
- The Cumpson–Seah CS2 formula (Surf. Interface Anal. 25 (1997) 430) is **not implemented**: I could not verify its exact form from an
  accessible primary source (the NIST guide text available to me was garbled at that equation). Users can enter L directly.

### 4.2 Single overlayer (Hill et al., Chem. Phys. Lett. 44 (1976) 225)
`I_o = I_o∞(1 − e^{−d/(L_o cosθ)})`, `I_s = I_s∞ e^{−d/(L_s cosθ)}`, `R = I_o/I_s`, `R∞ = I_o∞/I_s∞` (R∞ includes the
λ, density and sensitivity ratios). One angle: L_o = L_s ⇒ `d = L cosθ ln(1 + R/R∞)`; otherwise R(d) is monotonic and solved by bisection
to 1e-15 relative (the equation the Thickogram — Cumpson, Surf. Interface Anal. 29 (2000) 403 — solves graphically).
Several angles: weighted least squares `Σ((Rᵢ − R(d, θᵢ))/σᵢ)²` in d (optionally R∞) by Gauss–Newton with step halving, analytic
derivatives; `σ = √diag((JᵀWJ)⁻¹)` and `σ·√max(φ,1)`. Hill plot: `ln(1 + R/R∞)` vs `1/cosθ`, slope d/L (only for L_o = L_s).

### 4.3 Multilayer (≤ 3 layers on a substrate)
Species k in layer m (substrate m = N): `I_k(θ) = I_k∞ [1 − e^{−d_m/(L_{k,m} cosθ)}]^{(m<N)} Π_{j<m} e^{−d_j/(L_{k,j} cosθ)}`.
Fit of the free thicknesses with any angle-dependent instrument factor removed: per angle `G_a = Σ_k wI M / Σ_k w M²` (profiled),
residuals `(I − G_a M)/σ`; box-constrained Levenberg–Marquardt (d ∈ [0, 1000] nm, forward-difference Jacobian); covariance from the
profiled Jacobian; dof = K·A − p − A.

### 4.4 Relative depth plot
`RDP_k = ln( Σ_{grazing} I_k / Σ_{normal} I_k )`, default sets = the largest and the smallest angle; larger = nearer the surface; a common
angle factor shifts every RDP equally (only order/differences matter). Public description: Thermo Fisher application note "Angle resolved XPS"
(document D16069); see also Cumpson, J. Electron Spectrosc. Relat. Phenom. 73 (1995) 25 for what ARXPS can and cannot resolve.

### 4.5 Maximum-entropy depth profile (Livesey & Smith, J. Electron Spectrosc. Relat. Phenom. 67 (1994) 439; Smith & Livesey, Surf. Interface Anal. 19 (1992) 175)
Layers j = 0 … J−1 with tops `z_j = j·z_max/(J − 1)` (last semi-infinite), default `z_max = 3·max L·max cosθ`, J = 20.
Composition `c_kj = softmax_k(u_kj)` (closure Σ_k c = 1, positivity). Kernel `W_kaj = L_k cosθ_a (e^{−z_j/(L_k cosθ_a)} − e^{−z_{j+1}/(L_k cosθ_a)})`,
model `M_ka = Σ_j c_kj W_kaj`, `χ² = Σ ((I_ka − G_a M_ka)/σ_ka)²` with G_a profiled as in 4.3,
entropy `S = Σ (c − m − c ln(c/m))` (Skilling form; c floored at 1e-300 inside the logarithm), default model m = composition at the most
normal angle, constant with depth. Minimize `Q = χ²/2 − αS` over u with L-BFGS (memory 12, Armijo backtracking, analytic gradient
`∂Q/∂u_lj = c_lj (g_lj − Σ_k c_kj g_kj)`, `g = ∂Q/∂c = −Σ_a r_ka G_a W_kaj/σ_ka + α ln(c/m)`; envelope theorem for G).
**α**: discrepancy ("historic MaxEnt") criterion χ²(α) = N = K·A − A, found by continuation from α = 1e6 down by half-decades (each solve
warm-started from the previous), then 40 bisection steps in ln α warm-started from the larger-α end. Below the discrepancy α the problem is
ill-posed (solutions depend on the start); those solutions are never used as starts. If χ² = N is not reached the result is flagged
`notReached` (model error or σ underestimated).

### 4.6 Validation (test/test_arxps.js)
Reference agreement: Seah–Dench exact; overlayer solution vs brentq ≤ 4e-15; angle fits vs scipy least_squares (d, σ, R∞, χ²) ≤ 4e-9;
multilayer forward ≤ 6e-16, fit d ≤ 2e-11, σ ≤ 3e-8; RDP exact; MaxEnt objective ≤ 1e-12, gradient ≤ 1e-11, fixed-α solutions
(α = 10, 30, 300) |Δc| ≤ 5e-7, discrepancy α 4.2999 vs 4.2999.
Synthetic truth (Poisson counts, 300 realizations): overlayer 0.5 / 1.8 / 5 nm — mean 0.5000 / 1.7999 / 4.9995 nm, ±1.96σ coverage
94.7 / 95.7 / 94.3 %; end-to-end from the synthetic
VAMAS angle series (Shirley subtraction → reference fit → R(θ) → fit) d = 1.800 ± 0.001 nm (true 1.8); Hill plot gives 1.90 nm because
L_o ≠ L_s (use the exact fit). Multilayer (contamination 0.6 nm / oxide 1.5 nm, unknown angle factor, 150 realizations) coverage 95.3 / 94.0 %, mean z −0.05 / −0.15.
MaxEnt (1 % noise, 10 angles 0–70°, 10 realizations): uniform composition stays within 0.03; A overlayer on B, t = 0.8 / 1.5 / 3 nm →
equivalent thickness 0.91 ± 0.09 / 1.63 ± 0.11 / 3.12 ± 0.11 nm, 50 % crossing 0.60 / 1.33 / 2.84 nm (interfaces are smeared and the crossing is pulled toward the surface by ≈ 0.2 nm);
buried Gaussian layer (centre 2.2 nm) recovered amount +8 %, centroid 2.65 nm (+0.45 nm). Failure modes: a 1 nm layer at 7–8 nm
(> 2.3 L) is not located (the amount appears spread toward the default); σ underestimated 20× is flagged `notReached`.

---

## 5. Spectral imaging — `xps_image.js`
Stack `{nx, ny, nE, E, data: Float32Array pixel-major data[p·nE + e], p = y·nx + x, cfac}`; constructors from images at energies
(`fromImages`), energy-major arrays (`makeStack({layout:"energy"})`) and per-position spectra (`fromSpectra`, VAMAS MAP blocks).
- `windowImage`: Σ (or mean) of channels with lo ≤ E ≤ hi.
- `areaMap`: per-pixel area with background none / linear (endpoint means of n points; computed as a dot product with the §1.2
  coefficients) / Shirley (SPEC §2 iteration); σ map as §1.2.
- `atomicMaps`: n = A/(S·T·L) per pixel, at% and delta-method σ, n.d. rule of SPEC §4, NaN where Σn ≤ minTotal. `ratioMap`.
- `roiSpectrum`: mean spectrum of a mask/rectangle, σ of the mean `√(Σy/cfac)/n_pix`. `binStack`: b × b sums (Poisson preserved).
- `subtractBackgroundStack`: per-pixel linear/Shirley subtraction; afterwards use unweighted NNLS (the Poisson weights of the subtracted
  data are not the counting variance).
- `pcaStack`: §2.1 preprocessing on the fly, r × r Gram, Jacobi eigen-decomposition, scores `Uᵀx/s`; Malinowski statistics from all r
  eigenvalues. `pcaReconstruct` (noise reduction). `noiseFactorsStack`: §2.3 MLPCA on an evenly strided subset of 2048 pixels (default) or
  method `pca` on the full stack.
- `nnlsMap`: per-pixel §2.6 without shift (bg 0/1/2 sign-free; weight none or Poisson IRLS, 2 passes); coefficient, σ, reduced χ² and
  amount-fraction maps.

Validation (test/test_image.js): every map vs numpy ≤ 1e-11 (Shirley 1e-9 tolerance, observed 5e-16); PCA eigenvalues ≤ 1.4e-14 of λ₁,
loadings/scores identical, IND ≤ 1.3e-14, reconstruction ≤ 4e-8 (Float32 storage); per-pixel NNLS vs scipy bvls ≤ 1e-14.
Synthetic truth: Poisson-weighted state maps ±1.96σ coverage 95.6 / 94.8 % (1 920 pixels), mean reduced χ² 1.01; at% map z-scores mean
−0.01, sd 1.01; factor count 3 (2 states + varying background) from the MLPCA subset; PCA(3) reconstruction lowers the RMS error vs truth
4.2×; synthetic VAMAS MAP file → state maps within 3.5 % / 6.2 % of the generator areas.
Performance, 256 × 256 × 100, one thread (node 22, this container; a Worker runs the same code): window 11 ms, linear area 58 ms,
Shirley area 0.38 s, PCA (Poisson, 8 factors) 1.2 s, MLPCA noise test (2048-pixel subset) 1.4 s, PCA reconstruction 0.17 s,
NNLS 0.56 s (unweighted) / 1.2 s (Poisson IRLS), Shirley subtraction 0.37 s. Peak memory ≈ 3 × the stack (Float32 26 MB + copies).

---

## 6. API for the UI integrator
All functions are synchronous, pure and DOM-free. Arrays are `Float64Array` (stack data `Float32Array`). Load order in a page or Worker:
`xps_core.js` → `xps_multi.js` → `xps_arxps.js` → `xps_image.js` (globals `XPS_CORE`, `XPS_MULTI`, `XPS_ARXPS`, `XPS_IMAGE`; in node
`require`). Each of `XPS_MULTI` and `XPS_IMAGE` has `run({fn, args})` for a generic Worker dispatcher.

```js
// Worker (same pattern as the core worker: concatenate the four source texts into one Blob)
self.onmessage = e => { const {mod, fn, args, id} = e.data;
  try { const r = self[mod].run ? self[mod].run({fn, args}) : self[mod][fn](...args); postMessage({id, r}); }
  catch (err) { postMessage({id, error: String(err.message || err)}); } };
// main thread: worker.postMessage({mod:"XPS_IMAGE", fn:"nnlsMap", args:[stack, refs, {bg:1, weight:"poisson"}], id}, [stack.data.buffer])
```

**Series / depth profile (`XPS_MULTI`)**
| function | input | output |
|---|---|---|
| `seriesFromSpectra(spectra, {coord?, iVar?})` | XPS_IO spectra | `[{key, region, sample, coordLabel, coordUnit, coord, spectra}]` |
| `commonGrid(spectra)` | spectra of one series | `{E, Y:[col], resampled}` |
| `integrateSeries(E, Y, {lo, hi, bg:{type,n,…}, cfac})` | | `[{area, sArea, converged}]` per level |
| `modelSeries(E, Y, fitFn)` | `fitFn(E,y,j) → {components:[{name,area,sArea}]}` | per level `{components, raw}` |
| `compositionSeries(rows, areas, {imfp})` | rows `[{label, el, S, BE, hv, PE, tf}]`, areas[row][level] `{area, sArea}` | `out[level][row] = {at, sAt, nrm, sNrm, nd}` |
| `depthScale(t, rate)`; `infoDepthCorrect(z, X, lamE, {hw})` | X[element][level] = normalized n | `{n, at, warn}` |
| `subtractBackground(E, cols, {lo, hi, bg})` | | `{E, cols, i0, i1}` |

**Multivariate (`XPS_MULTI`)**
| function | output |
|---|---|
| `pca(cols, {scale:"none"\|"poisson"\|"area", center:"none"\|"mean", nKeep, method, alpha})` | `{r, c, U, s, V, eig, varPct, malinowski:{table:[{n,RE,IND,F,SL,REV}], nIND, nF}, pp, method, warn}` |
| `reconstruct(P, n)` | n-factor data in original units (noise reduction) |
| `noiseFactors(cols, P, cfac, {method:"mlpca"\|"pca", nMax, z, full})` | `{table:[{n, chi2, dof, ratio, limit, within}], nNoise, method}` |
| `wlra(cols, n, {cfac})` | `{rec, chi2, iter}` |
| `targetTest(P, n, x)` | `{t, xhat, AET, REP, RET, SPOIL, verdict, cos}` |
| `opa(cols, n)`, `simplisma(cols, n, f)` | indices (levels / energy channels) |
| `ittfa(cols, P, n, {needles})` | `{C, S, needles, iter}` |
| `mcrAls(cols, n, {init:"simplisma"\|"opa"\|"ittfa"\|[spectra], closure, maxIter, tol})` | `{S (unit area), C, lof, iter, warn}` |
| `fitReferences(cols, refs, {E, weight:"poisson"\|"none", cfac, bg:0\|1\|2, maxShift, shiftStep})` | per spectrum `{coef, sig, bgCoef, sigBg, chi2, dof, phi, mu, shift, fraction, sFraction}` |
| `nnlsGram(G, h, n, free)`, `nnlsCols(cols, b, wts, free)`, `svdCols`, `symEig` | numerical building blocks |

**ARXPS (`XPS_ARXPS`)**
| function | output |
|---|---|
| `seahDench(E, {kind, a, rhoKgM3})`, `atomSize(M, ρ, n)`, `ealFromImfp(λ, ω)` | nm |
| `overlayerThickness(R, θ, {Lo, Ls, Rinf})` | d |
| `fitOverlayer(θ[], R[], σR[], {Lo, Ls, Rinf, fitRinf})` | `{d, sd, sdScaled, Rinf, sRinf, chi2, dof, phi, perAngle, warn}` |
| `hillPlot(θ[], R[], σR[], {Lo, Rinf})` | `{x, y, slope, d, sd}` |
| `multilayerI(model, d, θ)`, `fitMultilayer(model, θ[], I[k][a], σ[k][a], {fix, d0, scale})` | `{d, sd, sdScaled, chi2, dof, phi, G, model}` |
| `relativeDepth(θ[], I[k][a], σ, {normal, grazing, nSet})` | `{rows:[{k, rdp, s, rel}], order}` |
| `maxentProfile(θ[], L[k], I[k][a], σ[k][a], {J, zmax, m, alpha})` | `{z, dz, c[k][j], alpha, chi2, N, S, G, curve, notReached, model, warn}` |

**Imaging (`XPS_IMAGE`)** — see §5 for definitions.
`makeStack({nx, ny, E, data, layout, cfac})`, `fromImages(images, nx, ny)`, `fromSpectra(spectra)`, `windowImage(st, lo, hi, mode)`,
`areaMap(st, {lo, hi, bg, n})` → `{area, sArea, converged}`, `atomicMaps(items, {minTotal})` → `{at, sAt, labels}`, `ratioMap(a, b)`,
`roiSpectrum(st, mask|rect)` → `{E, y, sum, n, sigma}`, `binStack(st, b, mode)`, `subtractBackgroundStack(st, {lo, hi, bg, n})`,
`pcaStack(st, {scale, center, nKeep})` → `{U, s, eig, varPct, scores, malinowski, …}`, `pcaReconstruct(st, P, n)`,
`noiseFactorsStack(st, P, {method, nSub, nMax})`, `nnlsMap(st, refs, {bg, weight})` → `{coef, sig, bgCoef, redchi, fraction, refAreas, note}`.

**Readers (`XPS_IO`)**: `readFile` (VAMAS → `r.vamas`, each spectrum `meta.vamas` as §3), `textMultiParse(text, name, {axis, hv})`.

### Running the checks
```bash
python3 -I test/ref/make_vamas.py         # synthetic VAMAS / text files + truth
python3 -I test/ref/ref_multi.py          # ≈ 3 s
python3 -I test/ref/ref_arxps.py          # ≈ 5 s
python3 -I test/ref/ref_image.py          # ≈ 3 s
node test/test_multi.js                   # full ≈ 1–2 min (QUICK=1 for a short run)
node test/test_arxps.js                   # full ≈ 1–2 min
node test/test_image.js                   # ≈ 20 s (includes the 256×256×100 timing)
```
