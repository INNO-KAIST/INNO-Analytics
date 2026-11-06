# AFM CSV Export Implementation Plan

> For agentic workers: use superpowers:subagent-driven-development and test-driven-development; verify each WBS before advancing.

**Goal:** Export AFM numeric planes in matrix or coordinate CSV with reproducible metadata.

**Architecture:** Extend the existing offline Image Analysis module, reuse its current pixel planes, original IBW source metadata and ZIP conventions. Export a ZIP containing CSV(s) and JSON metadata; explicit validation prevents unknown original or calibration from being silently inferred.

**Tech Stack:** Existing JavaScript, Node.js, Python and Playwright; no product runtime or build-tool additions.

**Spec:** E:/Develop/INNO Analytics/docs/PRD.md §12, approved 2026-11-06.

## Global Constraints
- Read source IBW read-only; no private data transmission, deployment, or raw file commits.
- Export current values or verified original source only; consistent snapshot, no source/history changes; -0/NaN/Inf policy and IEEE-754 roundtrip.
- Pixel-center coordinates (0,0) top-left, X right/Y down; anisotropic calibration if validated, else px; no invented physical origin.
- ZIP safe paths, unique names, complete selection validation, enforce memory/size limits; no partial download on errors.
- Per-WBS failing tests → implementation → focused regression → independent review → progress record; stop after two consecutive failed fixes of same issue.

## Review Focus
- Original source SI units, geometry history and stale/missing calibration must not be presented as verified original values.
- C/Z/T plane order, rows and pixel coordinate alignment, signed zero, floating boundaries and nonfinite.
- Invalid or unsupported mixed selections, path traversal, duplicate filenames, worker race and in-flight changes.
- Large images/ZIP32 limits, cancellation, no mutation or loss of current data/history.
- Independent CSV roundtrip and real IBW evidence; Excel rounding; original file and existing analysis unaffected.

### Task 1: CSV serialization and safe ZIP — CSV-01,02,03,04
**Files:** existing modules/image/src code per injector architecture; new focused Node tests under modules/image/test.
- [ ] Write failing tests for matrix and coordinate values, coordinates/scaling, Float64 roundtrip/-0/NaN/Inf, safe CSV cells, unknown metadata, pathological filenames, deterministic safe ZIP and explicit budget failures.
- [ ] Implement pure serialization and package preparation with existing offline conventions; do not mutate input. Run focused tests and core regression.
- [ ] Independent reviewer validates physical units, Excel limitations, source protection and resource limits; record result.

### Task 2: Source selection and export UI — CSV-01,02,03,04
**Files:** modules/image/src/ui and focused browser tests.
- [ ] Fail-first tests cover current vs verified original SI data, unavailable original/geometry/calibration failure, non-square and stack ordering, mixed RGB, duplicate names, asynchronous mutation, memory caps and no pixel/history changes.
- [ ] Implement selection UI, file packaging/download, correct metadata and abort behavior; independent review and browser regression.

### Task 3: Independent real IBW, full regression and integration — CSV-05
**Files:** test-only independent parser/reference, SPEC.md, docs/HISTORY.md, docs/사용 가이드.md, generated HTML and PROGRESS.md.
- [ ] Independently reparse both CSV formats and JSON for actual IBW source/current and test every pixel/coordinate, compare measured values and units; no network/raw writes.
- [ ] Run new + approved numerical/browser tests and relevant Fiji where available; record skips; two-pass HTML injection/hash.
- [ ] Independent final review, compare project baseline, backup exact overwritten files, integrate, repeat tests on actual HTML; update documentation and progress. No push/deploy without explicit request.

### Task 2A: Color provenance state across image lifetime (user-approved 2026-11-06)
- Write failing tests from split-color, split-snapshot, numeric-snapshot and mixed-input reproductions; add OIR/POIR, NPY, TIFF, PNG, ND2, CZI, LIF, AFM parsers, clone, geometric operations, channel split, RGB-to-gray, undo, legacy sessions and snapshots.
- Define validated numeric / color-components / display-gray / unknown semantics. Preserve single color components and use conservative multi-input propagation. Explicit grayscale is image intensity, never restored AFM quantity. Persist in .innoimg and Undo; old/malformed origin metadata is unknown.
- Handle ordinary TIFF/PNG and external/legacy reimports honestly: absent provenance cannot establish physical measurement origin; do not assert perfect recognition. Never relabel unverified content as measured AFM values.
- Independently review source format, old-file compatibility, serialization, altered files and real IBW roundtrip before resuming integration. User requests #6/#7 remain separate. After user-approved restart, same-error consecutive-fix loop count restarts; stop/report if threshold reached again.
