# XPS spectral pre-processing — specification (SPEC_proc v1.0)

This document is the definition. `src/xps_proc.js` (JS engine, global `XPS_PROC`) and `test/ref/proc_ref.py`
(independent Python reference) are both written from it. If code and document disagree, fix the document first.
Only published methods are used (references in §9); no vendor code or data.

Conventions follow `SPEC.md` §0: energy axis ascending inside the engine (BE ascending; KE is allowed where stated),
intensities in CPS or counts, areas in intensity·eV. "Uniform grid" means `|ΔE_i − Δ| ≤ 1e-6·Δ`.
Nothing in this module is used for quantification by default: deconvolution and smoothing are for display,
peak finding, start values and resolution checks. Quantitative positions/areas come from fitting the raw data
(SPEC.md §3), optionally with the PSF of §4 in the model.

---------------------------------------------------------------------------------------------------------
## 1. Noise model (`noiseModel`, `resolveVar`)

Second difference `D2_i = y_{i−1} − 2y_i + y_{i+1}` (i = 1..N−2); for independent noise `Var D2_i ≈ 6σ_i²`.
Local level `ℓ_i` = mean of `y` over the 5-point window `[i−2, i+2]` clipped to the data.
- Poisson-like (`kind:"poisson-like"`, default when all `y > 0`): `var_i = κ·max(ℓ_i, f)`,
  `κ = median_i[ D2_i² / (6·max(ℓ_i, f)) ] / 0.454936423119573` (median of χ²₁), `f = 1e-3·max|y|`.
  κ = 1 for counts; for CPS data κ = 1/dwell-equivalent.
- constant (`kind:"constant"`): `var_i = s²`, `s = 1.482602218505602 · MAD(D2)`, `s² ← s²/6`.
`resolveVar(y, noise)`: array → used as is; number → constant variance; `"poisson"` → `max(y_i,1)`; otherwise auto.
Limitation: curvature of sharp lines inflates D2; the median makes this robust when lines occupy < ~30 % of the region.

## 2. Resampling (`resample(E, y, opt)`)

Every method is a linear map written as sparse weight rows `y'_j = Σ_i W_ji y_i`, so the variance of independent
inputs propagates exactly as `var'_j = Σ_i W_ji² var_i` (PCHIP: first order, §2.3). Output samples are
correlated (neighbouring rows share inputs): do not treat resampled data as independent in a χ² fit without
accounting for this.
Default output grid for `{step}`: `E'_j = E_0 + jΔ'`, `j = 0..⌊(E_{N−1} − E_0)/Δ' + 1e-9⌋`. `method:"auto"` →
`rebin` if Δ' ≥ Δ, else `pchip`. Points outside the input range → NaN.

### 2.1 Flux-conserving rebin (down-sampling) — Carnall (2017)
Channel edges `e_0 = x_0 − (x_1 − x_0)/2`, `e_i = (x_{i−1}+x_i)/2`, `e_N = x_{N−1} + (x_{N−1} − x_{N−2})/2`
(same for the output grid). `y'_j = Σ_i y_i |B_i ∩ B'_j| / Σ_i |B_i ∩ B'_j|` (mean over the covered part of the
output bin; bins with zero cover → NaN; `cover_j` returned). Conserves `Σ y_i |B_i|` exactly when the output bins
cover the input. y is treated as an intensity density (CPS), not counts per bin.
Example: ×4 on aligned grids, each output bin = 3 whole + 2 half input channels, `Σw² = 3/16 + 2/64 = 0.21875`.

### 2.2 Linear (`linear`) — standard 2-point weights.

### 2.3 PCHIP (`pchip`, up-sampling, non-uniform OK) — Fritsch & Carlson (1980), Fritsch & Butland (1984)
Exactly the scipy `PchipInterpolator` rule: slopes `m_k = (y_{k+1} − y_k)/h_k`; interior
`d_k = 0` if `m_{k−1}·m_k ≤ 0` (sign change or zero), else `d_k = (w1+w2)/(w1/m_{k−1} + w2/m_k)` with
`w1 = 2h_k + h_{k−1}`, `w2 = h_k + 2h_{k−1}`; ends `d = ((2h_0+h_1)m_0 − h_0 m_1)/(h_0+h_1)`, set to 0 if
`sign d ≠ sign m_0`, else to `3m_0` if `sign m_0 ≠ sign m_1` and `|d| > 3|m_0|`. Cubic Hermite on each interval.
Variance: `d_k` is 1-homogeneous in the slopes, so its gradient rows reproduce the value exactly (Euler);
`var'` uses these rows (exact first-order propagation; at zero/sign-changing slopes PCHIP is not
differentiable and the one-sided branch is used).

### 2.4 Lanczos (`lanczos`, uniform input) — Duchon (1979)
`y(x) = Σ_i y_i L(u − i) / Σ_i L(u − i)`, `u = (x − x_0)/Δ`, `L(t) = sinc(t)·sinc(t/a)` for `|t| < a`
(`sinc t = sin πt / πt`), taps `i = ⌊u⌋−a+1 .. ⌊u⌋+a` inside the data; default a = 3.

### 2.5 Sampling adequacy
For a Gaussian line of FWHM w sampled at Δ, the relative spectral amplitude at the Nyquist frequency is
`aliasLevel(w,Δ) = exp(−π² σ² / (2Δ²))`, σ = w/(2√(2 ln 2)); `suggestStep(w, tol) = πσ / √(2 ln(1/tol))`
(tol = 1e-3 → Δ ≈ w/3.0). Point interpolation (`pchip`, `lanczos`) of data with Δ > w/2 aliases → use `rebin`.

## 3. Smoothing and derivatives (display / peak finding only)

### 3.1 Savitzky–Golay (`savgol(y, window, order, {deriv, delta, mode, cval, var})`) — Savitzky & Golay (1964), Steinier et al. (1972)
Equals `scipy.signal.savgol_filter(y, window, order, deriv, delta, mode, cval)` for odd `window`, `order < window`.
Coefficients: least-squares polynomial of degree `order` on `u_j = (j − c)/s`, `c = (window−1)/2`,
`s = max(1, c)`, computed through a modified-Gram–Schmidt QR (twice) of the Vandermonde matrix; the
`deriv`-th derivative at position `p` has weights `Σ_{m≥deriv} m!/(m−deriv)! ((p−c)/s)^{m−deriv} B_m / (s^deriv δ^deriv)`,
`B = R⁻¹Qᵀ`. `deriv > order` → zeros.
Interior: `y'_i = Σ_j c_j y_{i−c+j}`. Boundaries (scipy.ndimage conventions):
`mirror` (d c b | a b c d | c b a), `nearest`, `constant` (cval), `wrap`, and `interp` (default) = polynomial fitted
to the first/last `window` points evaluated at the first/last `c` positions.
### 3.2 Gaussian (`gaussSmooth(y, σ_pts, {mode, truncate})`)
Equals `scipy.ndimage.gaussian_filter1d(y, σ, mode, truncate=4)`: kernel `exp(−j²/2σ²)`, `|j| ≤ ⌊4σ + 0.5⌋`,
normalised to sum 1; default mode `reflect` (d c b a | a b c d | d c b a).
### 3.3 Wrappers
`smooth(E, y, {method:"sg"|"gauss", window|widthEV, order, fwhmEV, mode, var})`,
`derivative(E, y, {window=7, order=2, deriv=1})` = SG derivative with `delta = ΔE` (d^n y / dE^n, BE ascending).
Noise of a linear filter output: `Σ_j w_ij² var_j` (returned with `var`), white-noise reduction `noiseFactor(rows,i)`.

## 4. Point-spread function (`psfKernel(psf, Δ)`)

Bin-averaged kernel `h_j = ∫_{(j−½)Δ}^{(j+½)Δ} psf(x) dx`, `j = −R..R`, renormalised to `Σh = 1`; `h_j` moves
intensity from E to E + jΔ.
- `{type:"gauss", fwhm}`: `h_j = Φ((j+½)Δ/σ) − Φ((j−½)Δ/σ)` (erfc forms, no cancellation), `R = max(1, ⌈6σ/Δ + ½⌉)`.
- `{type:"voigt", wG, wL}`: Voigt = Re w(z)/(σ√2π), Faddeeva w by Weideman (1994) N = 40; bin integral by composite
  8-point Gauss–Legendre with `⌈4Δ/min(wG,wL)⌉` panels; radius `X = max(6σ, wL/(π·10⁻³))` (Lorentz tail mass ≤ 1e-3),
  capped at `max(10·FWHM_V, 6σ)`, `R = ⌈X/Δ + ½⌉`; truncated mass renormalised (`mass` returned).
- `{type:"measured", E, y, nb=3}`: subtract the straight line through the means of the first/last `nb` points, clip
  < 0, centre at the centroid, flux-conserving rebin (§2.1) onto offsets jΔ.
- `{type:"array", h}`: odd-length kernel already on the grid.
- From a Fermi-edge fit (§6) or reference line (§7): `result.psf` = `{type:"gauss", fwhm: wG}`.

## 5. Deconvolution (`deconvolve(E, y, opt)`)

### 5.1 Operator
Uniform ascending grid required. Half-sample-symmetric boundaries (scipy `reflect`, ρ):
`(H o)_i = Σ_j h_j o_{ρ(i−j)}`, adjoint `(Hᵀ r)_k = Σ_i Σ_j h_j r_i [ρ(i−j) = k]`. For symmetric h, H is symmetric,
diagonalised by the DCT-II with eigenvalues `λ_m = Σ_j h_j cos(π m j / N)`; `Σ(Ho) = Σo` exactly. (A clamped/edge
extension was tried first and made the iteration unstable at the edges; it is not used.)
Points within R of an edge (`edgePoints`) depend on the reflection assumption.

### 5.2 Iterations (o⁰ = y clipped to the constraint set)
- `jansson` (default): `o ← clip_[a,b]( o + r(o) · Hᵀ(y − Ho) )`, `r(o) = r₀ max(0, 1 − 2|o − (a+b)/2|/(b−a))`
  — Jansson (1970; Jansson, Hunt & Plyler 1970) relaxation with the reblurring of Kawata & Ichioka (1980).
  Defaults `a = 0` (positivity), `b = 4·max y` (`bFactor`), `r₀ = 1.5`. `reblur:false` gives the original
  `o + r(o)(y − Ho)` form (warned: in oversampled data it amplifies noise at frequencies where λ ≈ 0, because the
  Van Cittert filter `(1 − (1−λ)^k)/λ → k` there; with reblurring the filter → 0).
- `landweber`: `o ← P₊(o + r₀ Hᵀ(y − Ho))` (projected Landweber; r₀ = 1).
- `vancittert`: `o ← P₊(o + r₀(y − Ho))` (Van Cittert 1931; for comparison only).
- `rl`: `o ← o · Hᵀ(y / Ho) / Hᵀ1` (Richardson 1972; Lucy 1974; Poisson ML-EM, Shepp & Vardi 1982),
  `o⁰ = max(y, 10⁻⁶ max y)`, ratio := 1 where `Ho ≤ 1e-300`. Scale-invariant, flux-conserving, positive.

### 5.3 Stopping (noise-based)
`χ²_k = Σ_i (y_i − (Ho_k)_i)² / var_i`. With `iter` given, exactly `iter` iterations. Otherwise stop at the first k with
`χ²_k ≤ target` (`stop:"discrepancy"`, Morozov 1966), or when `χ²_{⌊k/2⌋} − χ²_k < 1` for k ≥ 50 (`"plateau"`:
doubling k gains less than one χ² unit), or after 20 consecutive increases (`"stalled"`), or at `maxIter` = 20000.
k = 0 is returned (with a note) when `χ²_0 ≤ target` already.
- `stopRule:"classic"`: target = τN (τ = 1).
- `stopRule:"projected"` (default, §5.4): target = `χ²_⊥ + τ(N − N_⊥)`.

### 5.4 Noise-only modes (`passband`)
`P_m = Y_m² / mean(var)`, Y = orthonormal DCT-II of y (`Y_m = s_m Σ_i y_i cos(πm(i+½)/N)`, `s_0 = √(1/N)`,
`s_m = √(2/N)`), so `E P_m = 1` for pure noise. `P̄_m` = mean of P over `[m−w, m+w]`, `w = max(2, ⌊N/100⌋)`;
`m_c` = first m with `P̄_m ≤ 1.5`. Noise-only set `⊥ = {m : |λ_m| < ε} ∪ {m ≥ m_c : |λ_m| < ½}` (ε = `passTol` =
0.01; λ from the symmetrised kernel). `χ²_⊥ = Σ_{m∈⊥} P_m`.
Reason: in oversampled spectra most DCT modes have λ ≈ 0; their noise can never be fitted, so `χ²_k` has a
realisation-dependent floor ≈ χ²_⊥, and the classic target N is unreachable for a sizeable fraction of noise
realisations (observed 3/12 with the doublet of §8 B1). Subtracting the measured unfittable part applies the
discrepancy principle to the fittable subspace only. With heteroscedastic noise `mean(var)` is the white level of the
DCT coefficients (approximation).

### 5.5 Diagnostics returned
- `resolution`: the same algorithm, same k, same bounds, applied to a noise-free blurred point source on a pedestal
  `median(y)` with blurred height `max y − median y`; `fwhmEff` = its FWHM above the pedestal; `gain = FWHM_PSF/fwhmEff`
  (`limitedBySampling` when fwhmEff < 1.5Δ). For non-linear methods this is the response at the data's level.
- `noise`: Monte-Carlo (M = `mc` = 6) re-deconvolution of `Ho_k + N(0, var)` with the same k; `sigma_i` = sample sd,
  `gain = √(Σσ_o²/Σvar)` over points ≥ R from the edges, `roughness = ‖D2 o‖/‖D2 y‖`. Warning when
  `gain > maxNoiseGain` (default 3).
- `areaIn/areaOut/areaRatio`, `chi2Hist`, `k`, `stop`, `target`, `passband`, `warnings`.

### 5.6 Known behaviour (measured in §8)
- Discrepancy stopping is conservative: the truth error keeps decreasing for 4–16× more iterations, at the cost of
  2–3× more noise amplification. It is the safe default for display; `iter` lets the user go further.
- Partial restoration of a blended pair pushes the two maxima apart (the effective resolution kernel has negative
  side lobes): +0.08–0.09 eV on a 0.70 eV doublet under a 0.80 eV PSF. Use deconvolved spectra to *see* and to *start*
  fits; take positions/areas from a fit of the raw data.
- Component areas from splitting at the valley are biased (ratio 2.15–2.41 for a true 2:1); total area is conserved.
- Jansson/Landweber may reach `maxIter` for wide PSFs at high count rates (slow modes with small λ); a warning says so.
- The Gaussian part of the PSF is assumed shift-invariant over the region.

## 6. Fermi edge (`fermiFit(E, y, opt)`, `fermiModel(E, p, opt)`)

`y(E) = b₀ + b₁(E − E_F) + ∫ (A + S(u − E_F)) f(u) G_σ(E − u) du`,
`f(u) = 1/(1 + exp(−s(u − E_F)/kT))`, s = +1 on a BE axis (occupied side at E > E_F), −1 on KE;
σ = wG/(2√(2 ln 2)); k = 8.617333262e-5 eV/K (CODATA 2018), T default 300 K (fixed).
Evaluation: f = H_s + e with H_s the step. Step part in closed form (z = (E−E_F)/σ, d = E − E_F):
s=+1: `A Φ(z) + S(dΦ(z) + σφ(z))`; s=−1: `A Φ(−z) + S(dΦ(−z) − σφ(z))`. The remainder
`e(v) = −s·sign(v)/(1 + e^{|v|/kT})` is integrated on `|v| ≤ 36kT` with 8-point Gauss–Legendre panels of width
`min(kT, σ)` (skipped when `|d| > 36kT + 9σ`).
Fit: Levenberg–Marquardt (Marquardt 1963) on `(E_F, wG, A, S, b₀, b₁)` (A, S, b₀, b₁ analytic Jacobian; E_F, wG
central differences, h = 1e-6(|p| + scale)). Start: E_F at the maximum of the SG-smoothed slope, A from the end
plateaus, wG from step height / maximum slope. Noise: `var` array (unscaled covariance); `"poisson"` (two passes:
weights 1/max(y,1) then 1/max(model,1), unscaled); default `"auto"` (same two passes, covariance × χ²_red).
`cov = (JᵀWJ)⁻¹`, 95 % interval = ±1.959964 σ. Returns `p, sigma, ci95, cov, chi2, dof, redchi2, model, EF, resolution
(= wG ± σ), psf, warnings` (edge too close to a boundary, unresolved wG < 2Δ, χ²_red > 3 in auto mode).
wG is the total Gaussian broadening (source + analyser + any other Gaussian contribution) at the edge's KE.

## 7. PSF from a reference line (`peakPSF(E, y, {wL, fitWL, noise})`)

`y = A·V(E − E₀; wG, wL) + b₀ + b₁(E − E₀)`; wL (core-hole lifetime width) fixed to the user value unless `fitWL`.
LM as in §6 (central-difference Jacobian), same noise options. `psf = {type:"gauss", fwhm: wG}`.
Systematics (warned): a too small wL inflates wG (0.633 vs 0.55 eV for wL 0.20 vs 0.30); the metallic DS asymmetry
of e.g. Ag 3d is ignored, so wG is an upper estimate for asymmetric lines. No literature wL is built in.

## 8. Validation (`node test/test_proc.js`, 76 checks, ~170 s; `QUICK=1` ~70 s)

A — engine vs `test/ref/proc_ref.json` (`python3 -I test/ref/proc_ref.py`, ~10 s, numpy 2.x / scipy 1.18):

| quantity | reference | max error |
|---|---|---|
| erf / erfc (55 points, |x| ≤ 6) | scipy.special | 3.9e-15 abs / 9.3e-16 rel |
| Voigt | scipy.special.voigt_profile | 3.1e-15 of peak |
| Savitzky–Golay, 472 cases (5 modes, window 5–31, order 0–5, deriv 0–3) | scipy.signal.savgol_filter | 8.2e-11 of max |
| SG variance | scipy impulse responses | 4.4e-15 |
| Gaussian smoothing, 40 cases (5 modes, kernel wider than signal) | scipy.ndimage.gaussian_filter1d | 6.9e-16 |
| PCHIP (non-uniform, flat run, sign changes) | scipy.interpolate.PchipInterpolator | 2.6e-16 |
| PCHIP variance (away from slope kinks) | finite-difference Jacobian of scipy | 1.6e-10 |
| rebin (+ variance, partial edges) | cumulative-integral formulation | 9.2e-16; area conserved 2.1e-16 |
| Lanczos a = 3, 5 | direct formula | 2.7e-13 |
| Gaussian / Voigt kernels | ndtr differences / quad(voigt_profile) | 4.9e-15 / 2.3e-15 |
| H, Hᵀ (reflect) | dense numpy matrices (np.pad symmetric) | < 1e-13 |
| noise κ, s² | numpy | < 1e-12 |
| Jansson, Landweber, Van Cittert, RL (10 cases, fixed k and both stop rules) | dense numpy iteration | identical k; o to 3.4e-14 |
| N⊥, m_c, χ²⊥ | scipy.fft.dct | identical / < 1e-9 |
| Fermi model (BE/KE, wG ≪ kT … ≫ kT, T = 30/300 K) | scipy.integrate.quad of the full integral | 2.6e-13 of max |
| Fermi fit / σ | scipy.optimize.least_squares (lm) | |Δp| ≤ 1.0e-7 σ; σ rel 4.6e-8 |
| Voigt reference-line fit / σ | least_squares | 6.2e-7 σ; 2.3e-7 |

B — synthetic truth (Poisson noise, fixed seeds):
- B1 doublet: two Voigt lines (wG 0.20, wL 0.0625 eV, FWHM ≈ 0.25) at ±0.35 eV, areas 3000:1500 on 300 counts,
  blurred by a 0.80 eV Gaussian (one maximum, FWHM ≈ 0.86 eV), Δ = 0.04 eV, 12 realisations, default stop:

  | method | resolved | k | separation (true 0.700) | main FWHM | point-source gain | noise gain | doublet area/true | Σo/Σy |
  |---|---|---|---|---|---|---|---|---|
  | jansson | 12/12 | 429–20000 | 0.784 ± 0.035 | 0.341 | 2.78 | 4.1 | 1.0025 | ≤ 3.2e-4 |
  | rl | 12/12 | 95–11454 | 0.782 ± 0.040 | 0.372 | 2.84 | 2.7 | 1.0018 | ≤ 1.1e-15 |
  | landweber | 12/12 | 180–20000 | 0.794 ± 0.047 | 0.485 | 2.01 | 3.5 | 0.9975 | ≤ 2.4e-4 |

- B2 stopping (one realisation): χ² falls monotonically to the target (χ²/N = 1.004 Jansson k = 839; 1.003 RL k = 146).
  RMSE vs truth: raw 1205 → 730 (Jansson) / 820 (RL) at the stop → 553 / 609 at 16·k / 4·k, while the noise gain
  rises 2.07 → 6.05 (Jansson) and 1.65 → 4.37 (RL); the over-iteration warning fires at 16·k and not at k/2.
  Featureless Poisson data stop at k = 0. Auto κ on Poisson data = 1.028.
- B3 Fermi edge (E_F 0.05 eV, wG 0.35 eV, 300 K, A 2000, S −300, b₀ 50, Δ 0.02 eV, 250 realisations):
  95 % coverage E_F 97.6 %, wG 98.0 % (Poisson weights; sd(z) 0.93/0.92), 96.4 %/96.4 % with auto scaling;
  bias −0.10 meV (E_F), −0.26 meV (wG). KE-axis edge at 1486.0 eV recovered within 1σ.
- B4: ×4 rebin keeps the area to 5e-7 and gives var/μ = 0.2183 (theory 0.21875); up-sampling a 0.9 eV line from
  0.25 eV: max error PCHIP 3.3 % (flattening at the maximum), Lanczos-4 0.29 %; SG derivative (21 pt, order 4) of a
  Gaussian 6.8e-5 of max.
- B5: measured PSF (off-grid, sloped baseline) keeps FWHM 0.6000 and centroid 7e-14 eV; deconvolution with it equals the
  analytic-PSF result to 6.6e-5; reference-line wG coverage 91 % (n = 100, auto noise).

## 9. API summary (`XPS_PROC`, browser global / `module.exports`, no DOM, Worker-safe)

```js
noiseModel(y, {kind, floor}) → {kind, kappa, s2, var}
resample(E, y, {step | Enew, method:"auto"|"rebin"|"pchip"|"lanczos"|"linear", a, var, fwhm, keepRows})
  → {E, y, var, method, cover, areaIn, areaOut, stepIn, stepOut, warnings}
aliasLevel(fwhm, step), suggestStep(fwhm, tol=1e-3)
savgol(y, window, order, {deriv, delta, mode:"interp"|"mirror"|"nearest"|"constant"|"wrap", cval, var}) → y | {y, var}
gaussSmooth(y, sigmaPts, {mode:"reflect"|…, truncate:4, cval, var}) → y | {y, var}
smooth(E, y, {method:"sg"|"gauss", window|widthEV, order, fwhmEV, mode, var}); derivative(E, y, {window, order, deriv})
psfKernel({type:"gauss",fwhm} | {type:"voigt",wG,wL} | {type:"measured",E,y,nb} | {type:"array",h}, step)
  → {h, R, step, fwhm, fwhmSampled, mass, type}
deconvolve(E, y, {psf | kernel, method:"jansson"|"rl"|"landweber"|"vancittert", noise|var, iter, maxIter:20000,
  tau:1, stopRule:"projected"|"classic", passTol:0.01, r0, bounds:{a,b}, bFactor:4, reblur:true, positive:true,
  mc:6, seed, maxNoiseGain:3, gainTest:true})
  → {E, o, model, resid, k, stop, chi2, chi2Hist, target, passband, resolution:{fwhmPSF, fwhmEff, gain},
     noise:{sigma, gain, roughness}, areaIn, areaOut, areaRatio, edgePoints, var, kernel, warnings}
fermiFit(E, y, {T:300, axis:"BE"|"KE", noise:"auto"|"poisson"|var[], p0, fix:{EF,wG,A,S,b0,b1}})
  → {p, sigma, ci95, cov, chi2, dof, redchi2, model, EF:{value,sigma,ci95}, resolution:{fwhm,sigma,ci95}, psf, warnings}
fermiModel(E, {EF, wG, A, S, b0, b1, T}, {T, axis})
peakPSF(E, y, {wL, fitWL:false, noise, p0}) → {p:{E0,A,wG,wL,b0,b1}, sigma, redchi2, model, psf, psfVoigt, warnings}
helpers: erf, erfc, faddeeva, voigt, gauss, lorentz, fwhm, gridInfo, uniformGrid, rng(seed){next,normal,poisson}, lm
```
UI notes: run `deconvolve` in a Worker for N > 500 or large k (cost ≈ 2·N·(2R+1) per iteration, plus `mc+1`
re-runs and an O(N²) DCT once); show `warnings`, `k`/`stop`, `resolution.gain` and `noise.gain` next to the plot;
draw `noise.sigma` as a band; mark `edgePoints` at both ends.

## 10. References
- P. A. Jansson, J. Opt. Soc. Am. 60, 184 (1970); P. A. Jansson, R. H. Hunt, E. K. Plyler, J. Opt. Soc. Am. 60, 596 (1970);
  P. A. Jansson (ed.), *Deconvolution of Images and Spectra*, 2nd ed., Academic Press (1997).
- S. Kawata, Y. Ichioka, J. Opt. Soc. Am. 70, 768 (1980) — reblurring in iterative restoration.
- P. H. van Cittert, Z. Phys. 69, 298 (1931). L. Landweber, Am. J. Math. 73, 615 (1951).
- W. H. Richardson, J. Opt. Soc. Am. 62, 55 (1972); L. B. Lucy, Astron. J. 79, 745 (1974);
  L. A. Shepp, Y. Vardi, IEEE Trans. Med. Imaging 1, 113 (1982).
- V. A. Morozov, Soviet Math. Dokl. 7, 414 (1966); M. Bertero, P. Boccacci, *Introduction to Inverse Problems in Imaging*,
  IOP (1998) — discrepancy principle, semi-convergence.
- A. Savitzky, M. J. E. Golay, Anal. Chem. 36, 1627 (1964); J. Steinier, Y. Termonia, J. Deltour, Anal. Chem. 44, 1906 (1972).
- F. N. Fritsch, R. E. Carlson, SIAM J. Numer. Anal. 17, 238 (1980); F. N. Fritsch, J. Butland, SIAM J. Sci. Stat. Comput. 5, 300 (1984).
- A. C. Carnall, "SpectRes: a fast spectral resampling tool in Python", arXiv:1705.05165 (2017).
- C. E. Duchon, J. Appl. Meteorol. 18, 1016 (1979).
- J. A. C. Weideman, SIAM J. Numer. Anal. 31, 1497 (1994).
- D. W. Marquardt, J. Soc. Ind. Appl. Math. 11, 431 (1963).
- W. Hörmann, Insurance Math. Econ. 12, 39 (1993) — PTRS Poisson sampler (tests only).
- S. Hüfner, *Photoelectron Spectroscopy*, 3rd ed., Springer (2003) — Fermi-edge broadening model.
