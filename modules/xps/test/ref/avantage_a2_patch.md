# Avantage Smart / Shirley backgrounds, Area(P), Area(N) and Atomic % — K-Alpha real data "a2"

**Applied in v1.9.2** (2026-10-29): the core patch is in `src/xps_core.js` (applied once; `--check` reports nothing to do), the §5 proposals are done (xps_io prop 415/416, UI tfFor, 'Avg. (eV)', Shirley pad, §4 warnings) — SPEC.md §2c and §4. `test_real_a2.js` now tests `src/` directly. Original instructions:

```
python3 -I test/ref/avantage_a2_apply.py            # edits src/xps_core.js (exact-string anchors, idempotent)
python3 -I test/ref/avantage_a2_apply.py --check    # dry run
python3 -I test/ref/avantage_a2_apply.py <copy>/xps_core.js
node test/test_real_a2.js                           # works before and after the merge (patches a temp copy of src/)
```

The apply script replaces four anchors in `xps_core.js`: the `bgShirleyAvantage` block, the `bgSmartAvantage` + `avantagePad` block,
one `makeProblem` line and the export list. If an anchor is missing or not unique, it prints the anchor and exits 1 without writing.
A step whose new text is already present is skipped. The INNO default Smart and Shirley are **not touched**. Only `conv:"avantage"` changes.

Data: `test/data/real/a2/<two sample folders>` (private, gitignored). This file contains numbers and region names only; the two samples are called S1 and S2 (S2 is the one with a Co2p region; `test_real_a2.js` prints the folder names locally).
Avantage settings as reported: Smart and Shirley use Background Average at Start&End = 0.50 eV. Shirley uses Iterations 50 and Convergence 0.01.
The instrument is a K-Alpha with Al Kα mono, PE 50 and a 0.05 eV step. The TXF was removed at source.

---------------------------------------------------------------------------------------------------
## 1. Result

| | Smart, 12 regions | Shirley, 12 regions |
|---|---|---|
| before (v1.8 `conv:"avantage"`) | Area(P) median 0.46 %, Co3p −6.3 % | typical 3–7 % |
| **after (this patch)**: max\|ΔB\|/H, median | 4.2e-4 | 1.4e-5 |
| **after**: \|ΔArea(P)\|, median · max | 0.006 % · 0.26 % (S2 O1s) | 0.005 % · 2.78 % (S2 Co3p) |
| regions meeting max\|ΔB\|/H ≤ 1.5e-3 **and** \|ΔA\| ≤ 0.1 % | 10/12 | 11/12 |

Here H is the peak height above the Avantage background, and Area(P) is the trapezoid integral over the window channels.

| region | Smart ΔB/H | Smart ΔA | Shirley ΔB/H | Shirley ΔA |
|---|---|---|---|---|
| S1 Ba3d | 1.0e-4 | −0.004 % | 1.5e-5 | −0.005 % |
| S1 Ba4d | 1.2e-6 | −0.002 % | 1.1e-6 | −0.002 % |
| S1 C1s  | 5.5e-4 | −0.006 % | 2.8e-6 (linear) | −0.010 % |
| S1 Co3p | 1.8e-4 | 0.061 % | 6.3e-4 | −0.054 % |
| S1 Ta4f | 7.6e-4 | 0.021 % | 9.7e-6 | 0.007 % |
| S1 O1s  | **2.5e-3** | −0.033 % | 1.8e-5 | −0.003 % |
| S2 Ba4d  | 3.0e-4 | 0.004 % | 1.1e-6 | 0.003 % |
| S2 C1s   | 3.5e-6 | 0.004 % | 3.5e-6 (linear) | 0.004 % |
| S2 Co2p  | 4.2e-4 | 0.008 % | 6.5e-4 | −0.055 % |
| S2 Co3p  | 1.3e-3 | −0.005 % | **6.2e-3** | **−2.78 %** |
| S2 O1s   | **5.7e-3** | **−0.26 %** | 1.4e-5 | 0.001 % |
| S2 Ta4f  | 1.0e-4 | 0.002 % | 5.9e-6 | −0.002 % |

At% check: INNO `quantify` with TF none and TPP-2M matches Avantage Atomic % to ≤ 0.0054 at% (all four tables). This uses Avantage's Area(P).
With INNO's own Area(P) from this patch, the difference is ≤ 0.056 at%.
**The double transmission correction is confirmed.** Applying the VGD TF again moves at% by up to **9.76 at%** with the current "thermo" (ln) form.
With the Avantage log10 form the shift is up to 5.85 at%. See §5.

---------------------------------------------------------------------------------------------------
## 2. Findings (Avantage code was not used; definitions only from the Avantage manual)

### 2.1 Facts about the export
- VGD y equals the xlsx counts/s to within 4.8e-6.
- The non-zero Backgnd channels are exactly the channels nearest End BE and Start BE (24/24).
- Area(P) is the trapezoid of (y − Backgnd) over those channels. The difference divided by ∫|y − B| is < 5e-4 for 24/24 regions.
  S2 Co3p Shirley shows a net relative difference of 1.6e-3, only because its net area is small after positive and negative parts cancel (§4).

### 2.2 Endpoint rule (48/48 ends, within print precision; max rel. 3.3e-6)
`h = round(W / (2·step))` channels, where W is the "Background Average at Start & End" setting. For 0.50 eV at 0.05 eV this gives h = 5 (11 points).
The endpoint value is the **plain mean of y over channels c−h … c+h**. This mean includes channels outside the window and is **truncated at the spectrum edge**.
Here the window is channels 1…N−2, so only one outside channel exists. The low-BE end is therefore mean(y[0..6]), 7 points.
Examples:
- S1 Ba3d: 56038.31 (Avantage 56038.3); the old 9-point value was 55446.0.
- S2 Ba4d: 13842.77 (Avantage 13842.8).
The same rule holds for Smart and Shirley, which have identical end values.

Kratos data (`test/data/real/kara`): the step there is also 0.05 eV, but its ends are 9-point (h = 4). The step size therefore does **not** explain the old 9-point finding.
Under this rule those files were processed with W ≈ 0.40 eV (any 0.35–0.45 eV), not 0.50 eV. Check this in Avantage if the files are reopened.
Without `avgEV`, `ah` or `rule:"v19"`, the v1.8 code runs unchanged (§6), so existing projects, the kara test and the SPEC §2c reference tests do not move.

### 2.3 Smoothed data S used inside the backgrounds
- **Smart:** S is a ±h moving average over the **window channels only, with the low-BE end channel excluded** (channels i0+1 … i1).
  Every interior Smart touch point (background = data; 68 points) equals this S to print precision.
  Example: S1 Ba4d has kinks at channels 4 and 6. These equal mean(y[2..9]) and mean(y[2..11]), which is S at the truncated lower bound.
  With the full-spectrum S the kinks are not touch points.
- **Shirley:** the integrand is a ±2h (21-point) moving average with the same truncation.
  A free 29-tap kernel regression of ΔB_k on y_{k±14} (Ba4d) returns a flat kernel on exactly −10…+10.
  The fit rms is 0.025 for 21 points, against 0.33 for 19 or 23 points.

### 2.4 Shape: not a converged Shirley
Inside every Avantage segment, the background satisfies, to print precision:

  u = p·u₁ + (1−p)·G(q·u₁ + (1−q)·u)

- G(v) = Δ·Q[max(0, s − v)] / Q_N, where Q is the inclusive cumulative sum from the low-BE end.
- s = S − (low end value) and Δ = high end − low end.
- u₁ = G(0) is the first ("one-step") Shirley.
- **p = 0.0700, q = 0.25** (fit; `bg.avP`, `bg.avQ`).

Evidence:
- **Unconstrained fit of u on {u₁, u₂ = G(u₁), G(u)}.** The u₁ weight is 0.0694–0.0700 in all 10 non-linear Shirley regions.
  Max residual 0.05–0.08 counts in Ba4d and Ta4f, which is print precision.
- **The v1.8 "λ" puzzle.** The old Kratos fit needed an integrand S − λB − (1−λ)y_lo with λ = 0.926–0.932, and λ was not constant.
  This is the same law: 1 − λ = p, and the q term makes the apparent λ vary.
- **Shirley clamp.** In the O1s and Ta4f tails, Avantage's Shirley lies above the data, yet it still rises in proportion to (S − a), with no dependence on B.
  This is the p·u₁ term, while the max(0, ·) clamp zeroes the G term. In Ba4d the low-end channels with S < a have exactly zero increments.
- **Plain Shirley iterations do not match.** These were tested on raw y, S5 and S10, with constant or linear start, rectangle, trapezoid or exclusive sums, and with or without clamp, for N = 1…60.
  Best results were 2e-2 to 1e-1 of Δ, so plain iterations (stopped at any N) are excluded.
- **Mechanism (not confirmed).** The closest mechanistic candidate found is a running mean of about 14 clamped iterates: input = mean of previous outputs, output = mean.
  Its weights give p ≈ 1/14 and q ≈ H₁₃/13 ≈ 0.245, but it reproduces the data only to 1e-4 – 3e-3 of Δ.
  The fixed-point form above is more accurate and is what the patch uses. Avantage's actual iteration and stop rule were **not identified**.
  No convergence threshold on area or on background change selects a consistent N across regions.

### 2.5 Smart segmentation
- The Smart segmentation is a recursion on [i, j]:
  - If the high-BE end value is not above the low-BE end value, the segment is **linear**. This is the manual rule, and it matches C1s, the Co3p low-BE parts and others exactly.
  - Otherwise the segment uses the §2.4 shape with that segment's end values.
  - The segment splits at the channel of maximum (B − S) > 0, with value S there, and both halves recurse.
- Touch sets are identical to Avantage in S1 Ba4d, S1 Co3p and S2 C1s.
  S1 Ba3d (459), S2 Co3p (398) and S2 Ta4f (222) each miss one Avantage touch next to a window end or another touch. The other differences are:
  - Extra INNO splits within 1–14 channels of a window end: S2 Ba4d 2,3; S1 Ta4f 386,390; S2 Co2p 818; S1 C1s 3. Avantage leaves the background up to 20–100 counts above S there.
    The same happens in the Kratos data, at 1–34 channels from an end. No length, depth, threshold or area-convergence rule reproduced which violations Avantage tolerates.
    The tested candidates were relative versus absolute size, first, last or middle violation, and a greedy global split with an area stop.
    The area effect is ≤ 0.06 %.
  - **O1s tails (both samples):** some Avantage segments (S1 [235–277] and [281–367], S2 [225–268] and [268–365]) are not reproduced by any tested segment law. Errors are 13–38 counts against a step of 1.3–4.4k.
    This gives S1 O1s ΔB/H 2.5e-3 (ΔA −0.03 %) and S2 O1s ΔB/H 5.7e-3 (ΔA **−0.26 %**). These are the only Smart regions outside the target.
    In these segments the background runs 40–120 counts under a noisy plateau, so a hidden internal touch is likely, but none was found.

### 2.6 Shirley
- If b ≤ a, the background is linear (both C1s, exact).
- Otherwise the background is the §2.4 shape over the whole window, with **no splitting**. It therefore passes above the data freely (§4).
- The "Iterations 50 / Convergence 0.01" settings are **not** reproduced as such. The fixed-point form ignores `bg.it` and `bg.cl`.
  Under v19 the UI's "Iter./Conv." inputs have no effect (§6).

---------------------------------------------------------------------------------------------------
## 3. Recheck of the old Kratos-derived Smart rule (`test/test_real_kara.js`, 26 regions)

| rule | \|ΔArea(P)\| median | max | max\|ΔB\|/H median · max |
|---|---|---|---|
| v1.8 (9-point S, converged Shirley on S) | 0.41 % | 4.12 % (3 Fe2p) | 3.9e-3 · 5.3e-2 |
| **v1.9.2 with `ah:4` (h = 4)** | **0.044 %** | 4.98 % (1 C1s), 4.49 % (1 Fe2p) | 1.9e-3 · 2.3e-2 |
| v1.9.2 with h = 5 | 0.37 % | 6.61 % | — |

- With the new rule, 20 of 26 Kratos regions are within ±0.3 %. Fe2p in files 2–4 is now 0.29–0.39 % (was −3.6 to −4.1 %).
- The two outliers are file-1 C1s and Fe2p. Avantage leaves 13–17 channels above S near both window ends there, and INNO splits at them (§2.5).
- `test_real_kara.js` is unchanged by the patch (44/44) because it passes no `avgEV`/`ah` and so runs v1.8. `test_real_a2.js` §A7 checks the v19 numbers above on the same files.
- **h = 4 confirms the 9-point endpoint finding.** h = 5 is clearly worse, so the old 9 points come from the setting, not from the step.
- The old rule is the default when no `avgEV`/`ah` is given, and `bg.rule:"v18"` forces it. It runs the v1.8 code unchanged and uses only the nearest 4 pad channels.
  Recommendation: once the UI exposes "Avg. (eV)" (§5 proposal 3), use v19 for Kratos files too, with `ah:4` or the actual Avantage setting.

---------------------------------------------------------------------------------------------------
## 4. Regions where Avantage's curve is physically questionable

"Above S" counts the interior channels where the exported Avantage background lies above the ±5-channel smoothed data. "Neg." is the negative part of ∫(y − B) divided by its positive part.
The last column is a plain trapezoid Shirley on y: 50 iterations, stopping when the relative area change is ≤ 0.01.

| Avantage Shirley | above S (longest run) | max excess (cps) | neg./pos. area | plain Shirley 50 / 0.01 | Area(P) Shirley / Smart |
|---|---|---|---|---|---|
| S1 Ba3d | 161 / 458 (158) | 9252 | 3.6 % | converges, 3 it | 0.88 |
| S1 Co3p | 138 / 397 (99) | 527 | 9.7 % | 5 it | 0.73 |
| S1 Ta4f | 159 / 397 (108) | 1595 | 9.4 % | 4 it | 0.78 |
| S1 O1s  | 235 / 397 (235) | 4697 | 12.3 % | 4 it | 0.80 |
| S2 Co2p  | 204 / 817 (118) | 9510 | 4.3 % | 5 it | 0.85 |
| **S2 Co3p** | **283 / 397 (174)** | 748 | **55.5 %** | **does not converge** (total Q_N ≤ 0 in 14 of 50 iterations) | **0.28** |
| S2 O1s   | 242 / 397 (236) | 4836 | 11.6 % | 4 it | 0.80 |
| S2 Ta4f  | 228 / 397 (117) | 1966 | 13.8 % | 4 it | 0.70 |
| (Ba4d ×2, C1s ×2 linear) | ≤ 116 | ≤ 357 | ≤ 2.7 % | Ba4d 2 it | 0.90–1.00 |

All Avantage Shirley curves are monotonic. Only the two C1s curves decrease, because the manual's linear fallback applies.
The Smart curves stay ≤ S except for 1–4 channels near the window ends (max 102 cps, §2.5). None is above the data over an extended range.

- **S2 Co3p Shirley is pathological.** The high-BE plateau (≈ 9.3–9.6k) sits below the 0.50 eV end average (10033), and the low-BE side (≈ 7.0–7.3k) sits below its end (7586).
  - The background is above the smoothed data on 72 % of the window. The net area 5212 is the remainder of +11.7k and −6.5k (cps·eV), only 28 % of Smart.
  - A standard Shirley iteration oscillates and never meets 0.01. This is the manual's "fail to converge if intensity falls within the background range" case.
  - The INNO model reproduces it to ΔB/H 6.2e-3 and ΔA −2.8 % (−145 cps·eV of 5212). The weights p and q are fitted on well-behaved regions, so the exact non-converged state is not captured.
- **Recommendation.**
  - Keep `conv:"avantage"` reproducing this behaviour for cross-checking. Do not fix it there, so that numbers can be compared with Avantage.
  - Have the UI warn when either of these holds:
    - `B.aboveFrac > 0.2` (the patch returns `B.above`, `B.aboveFrac` and `B.maxAbove` for Shirley);
    - the negative part of the net area exceeds 5 %.
  - Show a stronger warning ("Avantage's Shirley fails to converge here; use Smart") when aboveFrac > 0.5 or the net area is < 50 % of the Smart area.
  - By that rule, S2 Co3p, S1 O1s, S2 O1s and S2 Ta4f (aboveFrac 0.58–0.72 in the patch's own S(±2h)) get the strong warning.
  - The INNO default backgrounds are unchanged.

---------------------------------------------------------------------------------------------------
## 5. Area(N), Atomic % and the double transmission correction

- **ECF.** Area(N) = Area(P)/(SF·TXFN·ECF) with TXFN = 1 on every row.
  ECF/λ_TPP-2M(KE = 1486.68 − Peak BE) is constant to a max log deviation of 1.8e-4 (all four tables), within Area(N) rounding.
  The material is ρ 6.22, Nv 10, M 183.3, Eg 0, the same one fitted to the Kratos files. INNO's default TPP-2M material gives up to 0.029 at%, and KE^0.6 gives up to 1.58 at%.
- **Atomic %.** Atomic % = Area(N)/ΣArea(N) to ≤ 0.005 at%. INNO `quantify` reproduces it to ≤ 0.0054 at% with:
  - Avantage Area(P);
  - SF from the table;
  - `tf:{model:"none"}`;
  - `imfp:{model:"tpp2m", mat as above}`.
- **Double TF.** The VGD stores prop 415 = true (Acq_TXFN_Correction) and 416 = 3, alongside TXFN coefficients 4.288, 0.8466, −0.7506 and 0.1002.
  - The UI's `tfFor` in "auto" mode returns `sp.meta.tf`, so `quantify` divides by T(KE) a second time.
  - T(KE) runs from 11.55 (Ta4f) to 22.75 (Ba3d, Co2p), a factor of 2 across the elements.
  - The at% error is up to **9.76 at%** with the "thermo" ln form (for example S2 Co2p 38.14 → 28.49), and up to 5.85 at% with the Avantage log10 form.
  - **Confirmed:** with a TXF-removed VGD, the TF must not be applied again.

### Proposal (not part of the apply script; xps_io.js and xps_ui.js not edited)
1. `xps_io.js`, in `vgdParse` after `const tfc=…`:
   ```js
   const txfnCorr=gv(415), txfnMode=gv(416);   // Acq_TXFN_Correction (bool), mode (3 seen with "TXF Shape Correction, retaining Intensity")
   ```
   Then add the following to `meta`:
   ```js
   tfRemoved:txfnCorr===true, tfMode:txfnMode
   ```
   Keep `meta.tf` (the coefficients stay informative) and add `src:"… (removed at acquisition, prop 415)"` when `tfRemoved`.
   Document in the header comment: `415 Acq_TXFN_Correction · 416 TXFN mode`.
2. `xps_ui.js` `tfFor(sp)`, first line:
   ```js
   if(q.model==="auto"&&sp.meta&&sp.meta.tfRemoved) return {model:"none",note:"VGD: TXF removed at acquisition (prop 415)"};
   ```
   Show that note in the quant panel. Optionally warn if the user forces "profile" or "power" on such a spectrum.
3. `xps_ui.js` `bgOf` (line ~142): pass the pad for Avantage Shirley too, with the right width:
   ```js
   const av=p.bg.conv==="avantage"&&(p.bg.type==="smart"||p.bg.type==="shirley");
   … _pad:C.avantagePad(sp.y,R.i0,R.i1,C.avantageHalf(p.bg,R.E))
   ```
   Replace the Avantage Shirley "Iter./Conv." inputs, which v19 ignores, with one "Avg. (eV)" input bound to `bg.avgEV` for both Smart and Shirley. Show the resulting `h` and the point count.
   Leave it empty to keep the v1.8 rule. Its tooltip should say to enter the Avantage "Background Average at Start & End" value.
   Add the warnings from §4 next to the existing "Avantage 규약 … Linear" note (line ~974).

---------------------------------------------------------------------------------------------------
## 6. Exact replacement (what `avantage_a2_apply.py` inserts)

### API
- v19 is used **only if** `bg.avgEV > 0` (eV), `bg.ah` (channels) or `bg.rule:"v19"` is set (`avantageV19(bg)`). Otherwise, or with `bg.rule:"v18"`, the v1.8 code runs unchanged.
  This keeps the UI's Iter./Conv. Shirley and the SPEC §2c reference tests in `test_xps_casa.js` bit-identical.
- h = round(avgEV/(2·step)), or `ah`. Without either, h is 4.
- `bg.avP` and `bg.avQ` override p and q. `bg.avHs` overrides the smoothing half-width.
- The returned array carries:
  - Smart: `splits`, `smoothed`, `half`, `rule`, `ends`;
  - Shirley: `above`, `aboveFrac`, `maxAbove`, `fallback`, `iter`, `converged`.
- Both work for ascending and descending BE arrays. The result is invariant under reversal to < 1e-6.

### Anchor 1
Replace the whole `/* Avantage 규약 Shirley (SPEC §2c; …` comment and `function bgShirleyAvantage(E,y,bg){…}` (10 lines) with:

```js
/* Avantage 규약 Shirley v1.8 (SPEC §2c; Thermo Avantage Data System manual, Shirley Background): 반복 수 상한 it(기본 10, 1 = 단순 한 번),
   넓이 변화가 cl(상대) 이하이면 멈춤; 높은 BE(낮은 KE) 끝이 낮은 BE 끝보다 낮거나 같으면 Linear 로 바꿈. 끝점·적분은 INNO §2 그대로.
   v1.9.2 부터는 bg.rule === "v18" 일 때만 쓴다 (기본은 아래 bgShirleyAvantage — K-Alpha 실측으로 정한 규약). */
function bgShirleyAvantageV18(E,y,bg){ … unchanged body of the old bgShirleyAvantage … }
function avantageV19(bg){ … }          // rule "v19" | avgEV>0 | ah given (and rule ≠ "v18")
function avantageHalf(bg,E){ … }        // h from bg.ah | bg.avgEV/(2·step) | 4
function avantageOriented(E,y,bg,fn){ … }  // run on ascending-BE arrays, map back
function avantageEnds(y,pl,ph,h){ … }   // ±h means incl. pad, truncated
function avantageSmooth(y,hs){ … }      // ±hs mean over channels 1…N−1 (low-BE end excluded)
function avantagePQ(s,D,p,q){ … }       // fixed point u = p·u₁ + (1−p)·G(q·u₁ + (1−q)·u)
function bgShirleyAvantage(E,y,bg){ … } // v19 if avantageV19(bg), else V18
```

### Anchor 2
Replace the `/* Avantage 규약 Smart (v1.8, …` comment, `function bgSmartAvantage(E,y,bg){…}`, the `avantagePad` comment and `function avantagePad(y,i0,i1){…}` with
`bgSmartAvantageV18` (the old body, using the 4 nearest pad channels), the new `bgSmartAvantage`, and `avantagePad(y,i0,i1,h)`, where h defaults to 4.

### Anchor 3
In `makeProblem`, replace this line:
```js
  if(P.bg&&P.bg.conv==="avantage"&&P.bg.type==="smart") P.bg._pad=avantagePad(spec.y,r.i0,r.i1);
```
with:
```js
  if(P.bg&&P.bg.conv==="avantage"&&(P.bg.type==="smart"||P.bg.type==="shirley")) P.bg._pad=avantagePad(spec.y,r.i0,r.i1,avantageHalf(P.bg,r.E));   // v1.9.2: 끝값 ±h 평균에 창 밖 h 채널
```

### Anchor 4
In the export list, replace `bgShirleyAvantage,bgSmartAvantage,avantagePad,` with:
```
bgShirleyAvantage,bgSmartAvantage,avantagePad,bgShirleyAvantageV18,bgSmartAvantageV18,avantageV19,avantageHalf,avantagePQ,avantageSmooth,
```

The full new text, verbatim, is the `STEPS` table in `avantage_a2_apply.py`. The core code is reproduced here:

```js
function avantageHalf(bg,E){
  if(bg&&bg.ah!=null&&isFinite(bg.ah)) return Math.max(0,Math.round(+bg.ah));
  if(bg&&bg.avgEV>0&&E&&E.length>1){ const st=Math.abs(E[E.length-1]-E[0])/(E.length-1); if(st>0) return Math.max(0,Math.round(bg.avgEV/(2*st)+1e-9)); }
  return 4;
}
function avantageSmooth(y,hs){   // 창 안, 낮은 BE 끝(0번) 채널 제외
  const N=y.length, S=new Float64Array(N), c=new Float64Array(N+1); for(let i=0;i<N;i++) c[i+1]=c[i]+y[i];
  for(let i=0;i<N;i++){ const lo=Math.max(1,i-hs), hi=Math.min(N-1,i+hs); S[i]=hi>=lo?(c[hi+1]-c[lo])/(hi-lo+1):y[i]; }
  return S;
}
function avantagePQ(s,D,p,q){
  const n=s.length-1, lin=()=>Float64Array.from(s,(_,t)=>D*t/(n||1));
  if(n<1) return lin();
  const G=v=>{ const C=new Float64Array(n+1); for(let t=1;t<=n;t++) C[t]=C[t-1]+Math.max(0,s[t]-v[t]); if(!(C[n]>0)) return lin(); const k=D/C[n]; for(let t=0;t<=n;t++) C[t]*=k; return C; };
  const u1=G(new Float64Array(n+1)), u=Float64Array.from(u1), w=new Float64Array(n+1); let it=0, conv=false;
  for(;it<5000;it++){ for(let t=0;t<=n;t++) w[t]=q*u1[t]+(1-q)*u[t]; const g=G(w); let m=0;
    for(let t=0;t<=n;t++){ const nu=0.5*u[t]+0.5*(p*u1[t]+(1-p)*g[t]); m=Math.max(m,Math.abs(nu-u[t])); u[t]=nu; }
    if(m<=1e-10*Math.abs(D)+1e-12){ conv=true; break; } }
  u.iter=it; u.converged=conv; return u;
}
// Smart: ends ±h (pad) · S = avantageSmooth(y,h) · recursion: linear if V[j] ≤ V[i], else avantagePQ(S−V[i], V[j]−V[i]) ; split at argmax(B−S) > 0 with V[k]=S[k]
// Shirley: ends ±h (pad) · S = avantageSmooth(y,2h) · linear if b ≤ a, else a + avantagePQ(S−a, b−a) over the whole window (no split)
```

Run time: ≤ 0.23 s per region (S1 O1s Smart, 27 splits). Most regions take a few milliseconds.

---------------------------------------------------------------------------------------------------
## 7. Verification (patch applied to a copy of `src/`, Node 22)

| suite | unpatched copy | patched copy |
|---|---|---|
| test/test_real_a2.js | 22/22 (it patches its own temp copy) | 22/22 |
| test/test_real_kara.js | 44/44 | 44/44 (identical output: no `avgEV`/`ah` → v1.8) |
| test/test_xps.js | 153/153 | 153/153 |
| test/test_xps_casa.js | 111/111 | 111/111 (SPEC §2c Avantage-Shirley reference cases unchanged) |
| test/test_casa_v3.js | 69/69 | 69/69 |
| test/test_proc.js | 76/76 | 76/76 |
| test/test_arxps.js | 32/32 | 32/32 |
| test/test_image.js | 23/23 | 23/23 |
| test/test_multi.js (QUICK) | 50/50 | 50/50 |
| test/test_shapes2.js | 76/77 | 75/77 |

The `test_shapes2` failures are wall-clock limits ("≤ 10 s"; 10.2–11.2 s) on the LM-stall stress cases. They occurred while two suite runs shared the CPU, they appear in the unpatched copy too, and they do not touch Avantage backgrounds.
An earlier draft that made v19 the default without `avgEV` failed four §2c reference cases in `test_xps_casa.js`. This is why v19 is opt-in.
