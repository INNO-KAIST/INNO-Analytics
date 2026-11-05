# Line Flatten Implementation Plan

> **For agentic workers:** Use superpowers:subagent-driven-development and per-WBS verification. The controller coordinates and integrates; implementer, reviewer and validator are separate.

**Goal:** Provide 0/1/2-order line flattening with particle exclusion and synchronized non-destructive before/after previews, verified against real AFM data.

**Architecture:** Extend the existing AFM numerical core and operation registry rather than changing the offline single-HTML architecture. Use Float64 for computed planes, explicit fit-mask polarity and replay-safe mask identity; preview must use a cancellable non-committing computation and a shared display range.

**Tech Stack:** Existing JavaScript, Node, Python, Playwright, NumPy and existing test-only Gwyddion 2.71. No product runtime or build tool is added.

**Spec:** `docs/PRD.md` §10, approved 2026-11-05.

## Global Constraints
- Preserve source IBW, original sample count, calibration, units, existing ROI/undo/history and other modules.
- Nonzero finite particle-mask pixels are excluded from fitting but corrected in the chosen apply area. Invalid/misaligned masks fail explicitly.
- Scan direction X/Y and degree 0/1/2 must be explicit; per-line deficient or singular fits fail before any commit, and NaN/Inf data are preserved.
- Preview never mutates data, history, save state, ROI or source; stale asynchronous responses cannot be displayed or applied.
- User-provided AFM files are read-only; no external data transfer or deployment.
- Existing folder is not a Git repo; use an isolated temporary copy and retain verifiable backups before product integration.
- A WBS item must pass its tests and independent review before proceeding. Stop after two unsuccessful attempts to fix the same error.

## Review Focus
1. Mask polarity, calibration/plane alignment, missing or changed masks during replay.
2. Deficient/singular lines, finite input with overflow, large constant offsets and nonfinite values.
3. Stale preview following image/plane/mask/parameter changes, cancellation and overlapping operations.
4. 8/16-bit stack output, undo/replay, original IBW and no resampling/precision truncation.
5. Same LUT/range and synced view in preview, memory limits and real-data numerical evidence.

### Task 1: Line fit numerical engine — LF-01, LF-04, LF-05
**Files:** `modules/image/src/core/45_afm_level.js`; new `modules/image/test/test_line_flatten.js` and independent NumPy reference tests.
**Interface:** Extend `IMG_CORE.afmLevel(p,o)` with `method:"lineFlatten"`, `degree:0|1|2`, `axis:"x"|"y"`, `fitMask`/`applyMask` of 1=include. Return `{w,h,bit:32,data:Float64Array}`; never mutate input. Keep existing `rows` behavior.
- [ ] Write and execute failing tests: degree 0 equivalence to `rows`; analytic per-line `a[y]+b[y]*x+c[y]*x*x`, independent asymmetric exclusion masks and apply mask, axes, NaN/±Inf, singular/insufficient fit, overflow, empty regions and large baseline.
- [ ] Implement normalized, numerically stable per-line least squares without silently changing order. Guard memory and invalid inputs.
- [ ] Run new numerical tests, `node modules/image/test/test_afm.js`, `node modules/image/test/test_img.js`; independent reviewer checks all finite/nonfinite cases and compatibility.

### Task 2: Particle masks and safe commit/replay — LF-02, LF-04
**Files:** `modules/image/src/ui/05_ops.js`, supporting UI and new `test_line_flatten_mask_ui.js`.
**Interface:** UI chooses an existing binary/threshold/particle image as exclusion mask; finite nonzero means exclude. Reuse existing threshold and ROI editing tools to create/refine masks. Validate size, calibration and c/z/t alignment before conversion to fitMask. Record cryptographically strong fingerprint plus plane/geometry/polarity, and reject missing/changed/reordered masks on replay rather than silently reselect.
- [ ] Failing browser tests cover default no-mask, exclusion preserving a known particle, independent ROI/apply regions, invalid/nonfinite or mismatched masks, changed/deleted/reloaded mask, batch replay and immutable source.
- [ ] Implement through `setPlanes`, undo and `runWithUnit`; preserve existing 8/16-bit multistack separate-result behavior and recorded output semantics.
- [ ] Run focused browser and numerical tests; independent reviewer verifies mask identity, data preservation and replay.

### Task 3: Before/after preview — LF-03, LF-04
**Files:** Image Analysis UI/CSS under `modules/image/src/ui`, `modules/image/src/img.css`; new `test_line_flatten_preview_ui.js`.
**Interface:** Provide side-by-side before/after visualization with shared LUT/range and synchronized navigation, optional exclusion-mask overlay and explicit Apply/Cancel. Preview worker results carry image id/version, c/z/t, mask identity and parameters to reject stale responses. The apply path must recompute or validate an exact matching token.
- [ ] Write failing browser tests for pixel equivalence to final operation, no mutation or history/save writes, synchronized navigation, LUT/range equality, mask overlay, cancel, 8/16-bit multistack and undo.
- [ ] Test rapid parameter/image/plane/mask changes, stale worker completion and memory-limit failure; no stale result may display or commit.
- [ ] Implement with existing worker and canvas conventions; verify using desktop and narrow viewport with no horizontal overflow. Independent reviewer inspects interaction and concurrency.

### Task 4: Independent real-data validation and integration — LF-05, LF-06
**Files:** test-only validation scripts and `modules/image/SPEC.md`, `docs/사용 가이드.md`, `docs/HISTORY.md`, `docs/PROGRESS.md`; existing generated HTML.
- [ ] Independently compare 0/1/2-order and masked/unmasked real Asylum IBW channels (AC topography, SKPM, CAFM where usable) to per-line NumPy SVD least squares, pixel-by-pixel, with specified finite/nonfinite handling and tolerance `1e-8 nm + 1e-10*max(abs(reference))` where units are nm. Report modes and missing cases. If Gwyddion 2.71 built-in row correction can be driven safely with equal semantics, compare separately; otherwise record it as unverified.
- [ ] Run approved existing Node, AFM, and browser regression commands, including Fiji where available; record skipped comparisons, output and actual exit codes. Run all added focused tests.
- [ ] Independently review source changes and safety; back up actual project files, verify no intervening changes, inject single HTML twice and compare hashes. Run the relevant tests against actual integrated HTML in an isolated profile.
- [ ] Update documentation and requirement statuses. Do not publish or push without explicit user request.
