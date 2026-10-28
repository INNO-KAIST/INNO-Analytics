# XPS reference implementation (Python) — `test/ref/`

An independent Python 3 (numpy 2.x / scipy 1.18) implementation of `../../SPEC.md` v1. I wrote it from the
spec alone and did not read the JS engine (`src/`). Its purpose is adversarial validation: the JS engine and
this code should agree within the tolerances listed below. Where they disagree, either the spec is ambiguous or
one of the two implementations is wrong.

| file | role |
|---|---|
| `xps_ref.py` | library: lineshapes, backgrounds, model/fit, covariance and statistics, TPP-2M, transmission, at%, noise κ |
| `make_fixtures.py` | regenerates `fixtures/*.json` (`python3 make_fixtures.py all`, ~1 min) |
| `validate.py` | self-checks and third-party comparison; writes `fixtures/validation.json` (~3 min) |
| `fixtures/lineshapes.json` | 41 lineshape cases |
| `fixtures/backgrounds.json` | 40 real cases (8 regions × 5 types) + 7 synthetic cases |
| `fixtures/fits.json` | 7 synthetic fits with known truth + 3 real fits (C1s, O1s, Ta4f) |
| `fixtures/quant.json` | TPP-2M values, a TPP-2M literature check, toy at% for 4 T × 3 λ options |
| `fixtures/validation.json` | results of `validate.py` |
| `ref_casa_v2.py` → `ref_casa_v2.json` | CasaXPS 2.3.27 V2 readings (`casa_readings_v2.md`) re-simulated with scipy: CasaXPS Linear endpoints, unweighted LSQ on channels i0+1…i1−1, RMS = √Σr², SGL as height mix, LA(1,1,n) with wG = n·F/144.87 (`../../SPEC_casa_v2.md`; checked by `test_xps_casa.js` §5) |
| `ref_shapes2.py` → `ref_shapes2.json` | CasaXPS Cookbook 2026 line shapes (TLA, GL(m,n)/SGL(m,n) + T/PHI tails, ST/SB modifiers, DS(α,n[,l]), TP) written from the Cookbook formulas in x_c = −x/F; quad / FFT + cumulative Simpson / mpmath (`../../SPEC_shapes2.md`; checked by `test_shapes2.js`) |

All arrays are stored at full double precision. Energies are binding energies (BE) in ascending order, as the
spec requires internally. Real data are taken from `../data/bct600/*.txt` (sorted ascending, treated as CPS with
c = 1). Every fixture contains the sliced `E` and `y` arrays, so the loader and region slicing cannot be a source
of disagreement.

---------------------------------------------------------------------------------------------------------
## 1. Fixture schemas

### lineshapes.json
`{spec, generator, note, cases:[case]}`. Each case has these fields:

| key | meaning |
|---|---|
| `id` | e.g. `ls21_DS` |
| `shape` | `GLS`, `GLP`, `V`, `DS`, `LA` or `GLT` |
| `params` | GLS/GLP `{F,m}`, V `{wG,wL}`, DS `{alpha,wL,wG}`, LA `{F,alpha,beta,wG}`, GLT `{F,m,T,tau}` |
| `E0` | always 0. `x` is E − E0 |
| `x`, `f` | 201 points. `f` is unit-area (normalised over ℝ) **except DS**, see below |
| `area_R` | ∫ℝ f dx computed by adaptive quad; `null` for DS with α > 0 (diverges) |
| `fwhm_true` | numerical FWHM: maximum found by Brent's method, then half-maximum crossings by brentq at xtol 1e-14 |
| `x_max`, `f_max` | position and value of the maximum (DS, LA and GLT maxima are shifted away from 0) |
| `region` | DS only: normalisation window `[xmin,xmax]` relative to E0. `x` spans exactly this window |
| `norm` (DS) | `"region"` for α > 0. For α = 0 it is `"R …"`, the literal reading (see A1) |
| `f_region_norm` (DS, α = 0) | the same curve normalised over the region (the other reading of A1) |
| `region_integral_unnormalised` (DS) | ∫region (ds⊗G) — the normalising constant before dividing |
| `fwhm_OL` (V) | Olivero–Longbothum value |
| `J` (GLP) | J(m) |
| `N` (LA) | N |

### backgrounds.json
`{spec, cases:[case]}`. Each case has these fields:

| key | meaning |
|---|---|
| `id`, `file` | `file` is the source txt; `null` for synthetic cases |
| `type` | `linear`, `shirley`, `smart` or `tougaard` |
| `n_end` | 5 for real cases. 5 or 1 for synthetic cases |
| `params` | tougaard: `{kind:"U2",C:1643}` or `{kind:"U3",C:551,D:436}` |
| `E`, `y`, `B` | ascending BE inside the region, data, background |
| `endpoints` | `{y_lo,y_hi,E_lo,E_hi}` |
| `info` | shirley/smart: `{iterations, converged}`. tougaard: `{s, I_hi}` |
| `ys`, `smooth_h` (smart) | smoothed data used in the clamp, and the half-window h |
| `B_expected`, `max_abs_dev_from_expected` (synthetic) | the analytic expectation and the deviation from it |

The synthetic cases are as follows:
* `syn_flat_*`: y_lo = y_hi, so every type gives the constant y_lo.
* `syn_shirley_gauss`: y = 300 + 0.08·cumtrapz(P) + P with n_end = 1. This is an exact fixed point of the
  discrete Shirley iteration.
* `syn_linear_ramp`: a ramp plus a Gaussian. The linear background recovers the ramp exactly.

### fits.json
`{spec, note, cases:[case]}`. Each case has these fields:

* `id`, `kind` (`synthetic` or `real`), `file`, `region:[E[0],E[-1]]`, `c` (count factor), `E`, `y`
* `model`:
  ```
  { "background": {"type":"linear|shirley|smart|tougaard|none","n_end":5,"params":{...}}
                | {"type":"active","c0":P,"k":P,"s":P},
    "components": [ { "name", "shape", "E0":P, "A":P, <shape params>:P,
                      "doublet": null | {"delta":P, "r":P, "q":P} } ] }
  P = {"mode":"free","value":start,"lo":lo,"hi":hi} | {"mode":"fixed","value":v}
    | {"mode":"link","ref":<component name>,"k":1,"b":0}   # = ref's SAME parameter * k + b
  ```
  * A is the total area of a doublet. The partner sits at E0 + delta, has area A·r/(1+r), and its widths
    (F, wG, wL) are multiplied by q.
  * A ≥ 0 always holds, and active k ≥ 0.
  * The DS region is [E[0], E[-1]] of the fit data.
* `truth` (synthetic only):
  * `params` — the true values, in the same layout as `result.params`.
  * `mu` — the noiseless model.
  * `B` — the true background.
  * `seed` — the numpy `default_rng` seed. Data are `Poisson(c·mu)/c`.
* `result` — the final values after IRLS pass 2:
  * `params[comp][pname]` — every parameter, including fixed and linked ones.
  * `sigma[comp][pname]` — σ for free parameters only.
  * `free` and `theta` — the ordered free vector. `sigma_vec` and `cov` (p×p, `null` for excluded rows) use the
    same order.
  * `chi2`, `phi = chi2/(N−p)`, `N`, `p` (the number of all free parameters), `rho1`.
  * `ar1_factor = sqrt((1+ρ₁)/(1−ρ₁))` when ρ₁ > 0, otherwise 1.
  * `AICc`, `BIC`, `DW`, `runs_z`, `n_runs`.
  * `at_bound` — empty in every fixture.
  * `chi2_pass1` and `theta_pass1`.
  * `chi2_with_v_mu_final` — a diagnostic: χ² evaluated with v = max(μ_final, v_min)/c.
  * `uniqueness` — the fraction of starts within Δχ² < 1 of the best. It depends on the RNG, so do not compare it.
  * `n_starts`, `kappa_noise` (section 6), `hessian_cond`, `hessian_singular`.
  * `mu`, `B`, `v` — the weights actually used in pass 2.

### quant.json
* `hnu` — 1486.6 eV.
* `tpp2m_materials` — λ in Å at several kinetic energies for Au, Si, SiO₂, Ta₂O₅, Co₃O₄, BaTiO₃ and PMMA.
  Each entry includes the derived Ep, β, γ, U, C and D.
* `tpp2m_literature_check` — see §3.
* `toy_peaks` — four peaks with A, σ_A, S and BE.
* `toy_cases` — 12 combinations of T ∈ {none, power x=0.5, thermo a=[0,−0.45,0.05,−0.004], table (log–log)} and
  λ ∈ {none, power n=0.6, tpp2m (Ta₂O₅-like matrix)}, all at PE = 50 eV. Each result row gives KE, T, λ, n, at%
  and σ_at%. σ_at% comes from the delta method with independent areas:
  ∂at_i/∂n_j = 100(δ_ij Σn − n_i)/(Σn)².

---------------------------------------------------------------------------------------------------------
## 2. Numerical methods (reference side)

* **V** — `scipy.special.wofz`. The limits are handled exactly: wG = 0 gives L and wL = 0 gives G.
* **GLP** — J(m) = (π/√b)·erfcx(√(a/b)), or √(π/a) when m = 0.
* **GLT** — closed form.
  * G⊗e is the ex-Gaussian `(1/2τ)exp(−x²/2s²)erfcx((s/τ−x/s)/√2)`. For z < 0 it is rewritten with
    `erfcx(z)=2e^{z²}−erfcx(−z)` to avoid 0·∞.
  * L⊗e = (1/πτ)·Im[e^{z}E₁(z)] with z = −(x+iγ)/τ.
  * Both are checked against quad (agreement < 1e-12).
* **DS**
  * ds = Re[e^{−iπα/2}(γ−ix)^{α−1}]. The tail lies on the x > 0 side, as the spec requires.
  * The antiderivative is Re[i·e^{−iπα/2}(γ−ix)^α/α], or atan(x/γ) when α = 0. It is used for the region
    integral, so the inner integral is exact.
  * Gaussian convolution: fixtures use adaptive quad (epsrel 1e-12). Fits use a trapezoid on ±10σ with
    h = min(σ,γ)/40, which is exponentially accurate for an analytic integrand (agreement with quad < 1e-12).
* **LA** — N is analytic via `gammaln`. The convolution is done by quad (fixtures) or by a trapezoid with
  h = min(σ,F/2)/60 (fits). These agree to < 2e-9 relative to the maximum. The remaining error comes from the
  2nd-derivative kink at x = 0.
* **Backgrounds** are exactly the discrete operations in the spec.
  * Cumulative trapezoid on the data grid.
  * Shirley starts from B₀ = y_lo. It stops when max|ΔB| < 1e-7(|Δy|+1e-12), or after 200 iterations.
  * smart applies `min(B, ys)` after every update.
  * Tougaard I(E_k) is the trapezoid over E′ ∈ [E_0, E_k]. I(E_hi) is found by `np.interp`.
* **Smoothing (ys)**
  * A local quadratic LSQ fit on the window [i−h, i+h], truncated at the ends of the region.
  * h = max(2, round_half_up(0.25/ΔE)), where ΔE is the median |ΔE| of the region.
* **Fit**
  * `scipy.optimize.least_squares(method='trf')` with bounds, `x_scale='jac'` and ftol = xtol = gtol = 1e-15.
  * Starts: the spec start, plus 24 Latin-hypercube starts (12 for the slower DS and LA cases) over the
    nonlinear free parameters. Areas for every start come from NNLS weighted by the pass-1 variance.
  * Pass 1 uses v = max(ys, v_min)/c and keeps the best χ².
  * Pass 2 uses v = max(μ₁, v_min)/c. It is re-solved from the 4 best pass-1 solutions and the best is kept.
* **Jacobian** for covariance: central differences with one Richardson step, so the error is O(h⁴).
  * C = (JᵀV⁻¹J)⁻¹·φ.
  * Parameters within 1e-6·(hi−lo) of a bound would be dropped from the covariance. No fixture has one.
* **TPP-2M** — the spec formula, with λ in Å.
* **Transmission**
  * power: (PE/KE)^x.
  * thermo: exp Σ a_k ln(KE/PE)^k.
  * table: log–log interpolation, clamped outside the table.
* **Noise κ** — the section 6 formula, using the same `ys`.

---------------------------------------------------------------------------------------------------------
## 3. Validation results (`validate.py`, details in `fixtures/validation.json`)

* **Unit area** — max |∫ℝf − 1| over every integrable lineshape is about 3e-14.
* **FWHM**
  * GLS, GLP (m = 0, 1) and GLT (T = 0) give FWHM = F to < 1e-12.
  * For the Voigt, the numerical FWHM agrees with Olivero–Longbothum to ≤ 2.37e-4 relative. A scan of
    wL/(wL+wG) puts the worst case at 0.0237 %, near a ratio of 0.225, so the spec's "≤0.02 %" is slightly
    optimistic. This affects display only.
* **Limits** (max relative deviation): every limit below agrees to ≤ 2e-8 (the wG, wL = 1e-7 limits) or to
  ≤ 1e-12 (the others).
  * V(wG→0) → L and V(wL→0) → G
  * LA(1,1,0) → L and LA(1,1,wG) → V
  * DS(α=0)/π → V and DS(α=0, wG=0)/π → L
  * GLP m = 0 → G and m = 1 → L
  * GLT T = 0 → GLS
* **Convolutions** (fine grid or closed form vs quad): DS ≤ 3e-16, LA ≤ 2e-9, GLT tail ≤ 1e-13. The DS region
  integral (analytic inner integral vs nested quad) agrees to ≤ 1e-10.
* **Shirley**
  * Converges in 5–13 iterations on every real region.
  * The fixed-point residual is < 1e-7·|Δy|.
  * The `syn_shirley_gauss` background is reproduced to 8e-7 absolute, which is the stopping tolerance.
  * The flat and ramp cases are reproduced exactly.
* **Third-party comparison**, n_end = 1 because the packages only use single end points:
  * **lmfitxps 1.x `shirley_calculate`** agrees with this Shirley on all 8 real regions to ≤ 1e-13 (relative to
    the step height).
  * **pynxtools-xps `Shirley`** differs by 0.14–0.33 of the step height. It is a different algorithm, not a
    tolerance issue: it integrates y − B with B starting at 0, i.e. without removing the end-point offset, so its
    k and B are not the Proctor–Sherwood ones. Reported as INFO only (corrected 2026-10-15).
  * **lmfitxps `tougaard_calculate`** (U2: tc=1643, tcd=1, td=0; U3: tc=551, tcd=−1, td=436) agrees with this
    Tougaard as well, except where s < 0 (C 1s windows), which SPEC v1.1 clamps to 0 and lmfitxps does not. On a uniform grid, its rectangle sum is identical to the trapezoid, because K(0) = 0 and
    y(E_lo) − y_lo = 0. Its iterative scaling stops at 1e-6, so agreement is ~1e-6.
* **Fits**
  * A bounded Powell polish started from every reported optimum never lowers χ² by more than 1e-6.
  * |∇χ²|·σ ≈ 0 for every interior parameter.
  * A re-run with a different seed and 100 starts (30 for DS) finds no lower χ².
  * The DS model evaluated with quad and with the fine grid agrees to ~1e-12.
* **Truth recovery** (synthetic fits, within about 2σ).
  * Recovered:
    * syn3 (Voigt doublet + active background, q = 1.142 ± 0.013 vs 1.15)
    * syn4 (DS α = 0.121 ± 0.011 vs 0.12; region area 8085 ± 220 vs 8000)
    * syn5 (LA)
    * syn6 (active Shirley)
    * syn7 (GLP + GLT; τ 0.98 ± 0.05 vs 1.2 is 4.7σ off, a GLT T/τ correlation)
  * syn1 and syn2 use **static** backgrounds defined from end points, and they underestimate areas by 1.8–5 %.
    This is a property of the spec's background model, not a fitting bug: the Lorentzian tails outside the
    window and inside the end-point averages are lost (see P4).
* **TPP-2M**
  * The equations are identical to Tanuma, Powell & Penn, *Surf. Interface Anal.* **35** (2003) 268, eqs. (4)
    and (6).
  * Ep reproduces the printed Table 1 values for all 15 elements to ≤ 0.005 eV, which is the print rounding.
  * λ reproduces the published TPP-2M values for Si, Cu and Au at 1, 3 and 5 keV, as quoted in Jaber et al.,
    *Jordan J. Phys.* **16** (2023) 381, Table 1, column "Ref. [22]": Si 2.386 / 5.787 / 8.868 nm against
    2.39 / 5.79 / 8.87 nm, and Cu and Au within the 2-significant-figure rounding of the print.
  * Against the optical-data IMFPs of TPP VII (Table 4), the RMS deviation is 1–25 % for most elements, with
    diamond at 73 %, graphite at 50 % and Cs at 38 %. This matches the known TPP-2M failures for these solids
    (TPP IX, *SIA* 43 (2011) 689, which reports 12.3 % average RMS, or 9.2 % without diamond, graphite and Cs).
  * Caveat: I read the table values through a text extractor, so verify them against the paper before relying
    on them for anything more than a sanity check.

---------------------------------------------------------------------------------------------------------
## 4. Suggested JS-vs-reference tolerances

| fixture | quantity | tolerance |
|---|---|---|
| lineshapes GLS/GLP/V | f | 1e-10 × max f (analytic; Faddeeva implementations differ ~1e-13–1e-10) |
| lineshapes DS/LA(wG>0)/GLT | f | **1e-3 × max f** (spec 1.8 discretisation) |
| lineshapes | fwhm_true | 1e-6 abs (analytic shapes), 2e-3 rel (convolved) |
| lineshapes | area_R | 1e-9 (analytic), 1e-3 (convolved) |
| backgrounds linear | B | 1e-12 rel |
| backgrounds shirley/smart | B | 1e-6·\|y_hi−y_lo\| (stopping rule 1e-7; iterate count may differ by 1) |
| backgrounds smart | ys | 1e-9 rel (the same polynomial LSQ) |
| backgrounds tougaard | B | 1e-9 rel |
| fits (GLS/GLP/V only) | χ² | rel 1e-6; params \|Δθ\| ≤ 0.02σ; σ rel 1e-3; ρ₁, AICc, BIC abs 1e-3 |
| fits with DS/LA/GLT | χ² | rel 5e-3; params \|Δθ\| ≤ 0.2σ; σ rel 2 % (the JS shape error is up to 1e-3 of the peak) |
| fits | uniqueness, n_starts, seconds | not comparable |
| quant | λ, T, at%, σ_at% | 1e-10 rel |

If the JS χ² is **lower** than the reference χ² by more than the tolerance, read it as a reference failure (a
missed global minimum), not as a JS failure. Report it.

---------------------------------------------------------------------------------------------------------
## 5. Spec ambiguities (where the reference makes a literal choice)

* **A1 — DS normalisation at α = 0.** "α > 0 → region area" implies ℝ normalisation (= Voigt) at α = 0. That
  makes f discontinuous at α → 0⁺, because the region integral < the ℝ integral. A fit with α free and lo = 0
  sees a jump there. The reference is literal, and both variants are in the fixture.
* **A2 — DS region.** I took [Emin, Emax] as the first and last data points of the fit region, not E_lo/E_hi.
  The doublet partner uses the same region.
* **A3 — Smoothing.**
  * "간격" (spacing) is taken as the median |ΔE|.
  * `round` is round-half-up (JS `Math.round`; Python `round` is banker's rounding).
  * Smoothing is applied to the **region-sliced** data, not to the whole spectrum. The edge values depend on
    this choice.
* **A4 — smart with y_hi = y_lo.** The reference returns y_lo **without** the min(·, ys) clamp, because the
  zero-step rule is stated for shirley.
* **A5 — Shirley and Tougaard end points.**
  * Q(Emax) integrates to the last point, while y_hi is an n-point average. So B(Emax) = y_hi but B(E_hi) ≠ y_hi.
  * Tougaard is the reverse: B(E_hi) = y_hi via interpolation, and B(Emin) = y_lo.
* **A6 — Tougaard trapezoid.** It includes the E′ = E node, where K(0) = 0.
* **A7 — IRLS.** The reported χ² is the pass-2 objective with fixed weights v = μ₁ (from pass 1), not Pearson
  with μ₂. `chi2_with_v_mu_final` is provided as well. The spec does not say where pass 2 starts.
* **A8 — Bounds.** p counts all free parameters, even when some are dropped from the covariance at a bound.
  "At a bound" has no tolerance in the spec (the reference uses 1e-6 of the range). No fixture hits a bound.
* **A9 — Links.**
  * "Same parameter of the reference component" also covers doublet `delta`.
  * σ of linked parameters is not defined (the reference gives null; |k|·σ_ref would be natural).
* **A10 — NNLS start weights** are unspecified. The reference uses the pass-1 v.
* **A11 — hν** is not specified (1486.6 eV was used, from the VAMAS header). Table-T extrapolation is not
  specified either (the reference clamps).
* **A12 — "Uniqueness"** depends on the LHS RNG, so it is not reproducible across engines.

## 6. Physics concerns about the spec

* **P1 — LA side convention.**
  * The spec puts α on the **low**-BE side.
  * Published CasaXPS metal templates use α < β with the tail on the **high**-BE side, for example Ni metal
    LA(1.1, 2.2, 10) in Biesinger et al., *Appl. Surf. Sci.* 257 (2011) 2717. That implies that in CasaXPS α
    controls the high-BE side.
  * If so, CasaXPS parameters typed into this engine would give a mirrored asymmetry.
  * CasaXPS's third LA argument is also an integer Gaussian-width index, not a FWHM in eV.
  * The claim "같은 형태" (same form) should be checked against the CasaXPS manual before release.
* **P2 — Doublet q scales wG as well.** Coster–Kronig broadening is a lifetime (Lorentzian) effect. The
  instrumental Gaussian is common to both partners, so q should multiply wL (and F for GL shapes), not wG.
* **P3 — DS "area" is region-dependent.**
  * at% from DS areas changes with the window width, and is not comparable to ℝ-normalised areas of other shapes.
  * This needs a stated convention for quantification, for example a fixed tail cut-off in multiples of the
    width.
* **P4 — Static end-point backgrounds bias areas** when Lorentzian or GL tails extend outside the window. This
  is −1.8 % in syn1 and up to −5 % for the smaller peaks in syn2. The active Shirley recovers the truth (syn6).
  Worth a UI warning or a recommended window rule.
* **P5 — Tougaard s can be negative.** When y_hi < y_lo, s < 0 (C1s regions here, s ≈ −3900), which is an
  unphysical negative loss cross-section. Tougaard on ~10 eV windows is also outside its intended use, which
  needs ≥ 30–50 eV above the peak. Suggest clamping s ≥ 0, or warning.
* **P6 — AICc/BIC with model-dependent weights.** In pass 2, v = μ₁ differs between competing models. So
  "the same variance model" is not strictly satisfied, and ΔAICc between models carries a weight-change term.
  Comparing models with pass-1 (data-based) weights would be cleaner.
* **P7 — κ (section 6) is biased by curvature.**
  * The 2nd difference includes signal curvature, which is ~20 % at the tops of 1.2 eV peaks with 0.05 eV steps.
  * The median of d² from ~100–300 points has a 12–25 % sampling error.
  * On Poisson synthetics with known c, κ·c ranged 0.64–1.31.
  * Keep it display-only, as the spec says.
* **P8 — TPP-2M Nv convention.**
  * TPP-2M is very sensitive to Nv, the valence electrons per formula unit, which is a convention for d-metals
    and compounds.
  * Its stated validity range is 50–2000 eV (later tested to 200 keV).
  * It fails for diamond, graphite and alkali metals (above).
  * If S comes from Avantage/ALTHERMO-type factors, those already include T and λ, so applying T and λ again
    would count them twice.

---------------------------------------------------------------------------------------------------------
## 7. Resolved in SPEC v1.1 (2026-10-15, after this review)

* **P1** LA sides swapped: α is now the high-BE side (matches the metal templates, e.g. Ni LA(1.1, 2.2)).
* **P2** Doublet q multiplies the lifetime width only (wL for V/DS, F for GL shapes and LA); wG is common.
* **P3** DS areas stay region integrals; the quant panel warns when a DS component is included.
* **P4** The fit panel warns when a static background meets large Lorentzian tails (active Shirley recommended).
* **P5** Tougaard s < 0 is clamped to 0 with a warning (both references updated).
* **P6** AICc/BIC use χ²_d (data-only variance). Multi-start selection uses the pass-1 χ² as well — the JS
  engine used the pass-2 χ² across starts, which compares different weights (found with `fits.json`).
* **P8** ALTHERMO1/Scofield ratios show no KE trend (log-slope 0.01 over 744 lines), so ALTHERMO1 does not
  contain T·λ; T and λ are still applied.
* **A11** Table transmission: clamped outside the table (spec now says so).
* Variance floor: v = max(ŷ, v_min)/c exactly as in the spec (the JS engine had v_min·c).
