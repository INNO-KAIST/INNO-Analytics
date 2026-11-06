# AFM Export Naming Implementation Plan

> For agentic workers: use subagent-driven development and test-driven development, with independent review and verification at each WBS.

**Goal:** Use original filename and actual channel for PNG and CSV exports, distinguishing duplicates and planes.

**Architecture:** Reuse current Image Analysis metadata and offline JS; extract real source and channel without guessing, apply shared safe naming and ZIP collision resolution, retain metadata filename references and all numeric values.

**Tech Stack:** Existing JS, Node, Python and Playwright; no new runtime.

**Spec:** E:/Develop/INNO Analytics/docs/PRD.md §13, approved 2026-11-07.

## Global Constraints
- SKPM_AAA.ibw + Height Retrace → SKPM_AAA_Height Retrace.png; preserve actual channel spelling (do not invent spaces).
- Distinguish known duplicate images, same channel names and C/Z/T; avoid unsafe/reserved paths, case-fold/NFC collisions, extension confusion, separators, bounded names; metadata references must match ZIP entries exactly.
- Browser download folder collisions cannot be inspected or guaranteed; rely on browser policy. Preserve raw files and data, do not upload private raw or tests' sensitive data.
- Do not break CSV numeric/source precision, provenance, prior Image Analysis functionality. User has explicitly authorized push/deployment after verification.
- WBS item tests and independent review before next; same error two successive unsuccessful fixes => stop/report.

## Review Focus
- IBW filename with periods and extension, label with slashes, single-channel vs multichannel, SOURCE parent removed/reloaded.
- Duplicate labels across images, same image/planes, case-insensitive and Unicode normalization collisions; accurate metadata mapping.
- Flat ROI/zoom/corrections must not affect values, image identity or channel selection; legacy data may need fallback.
- Windows reserved names and file/path/ZIP security; offline browser file-origin and duplicate download behavior.
- Distinct stdout-free privacy checks for deployment; isolated fresh clone because existing deploy tool can reset and clean.

### Task 1: Safe naming and PNG/CSV connection — NM-01,02,03,04
**Files:** modules/image/src/img_io.js, ui/90_main.js, optionally shared existing UI/core and new tests.
- [ ] Failing tests for actual IBW, `Height Retrace`, no source, no channel, periods, C/Z/T and same name collisions; original filename vs channel, duplicate/snapshot/reload, forbidden/reserved/Unicode filenames; ZIP CRC/metadata and PNG download filename.
- [ ] Implement shared naming with actual source provenance and safe fallback. Ensure ZIP two-pass collision resolution and metadata correspondence; keep C/Z/T suffix only where needed and preserve values/units.
- [ ] Run focused JS/browser + existing CSV provenance/IBW tests, independent reviewer verifies scope, all collision paths and no silent metadata loss.

### Task 2: Independent real IBW, full regression, documentation & integration — NM-04
**Files:** relevant tests, modules/image/SPEC.md, docs/HISTORY.md, docs/사용 가이드.md, generated HTML, docs/PROGRESS.md.
- [ ] Real read-only IBW test actual multi-channel filenames/CSV values, PNG exact and other-source cases; raw SHA unchanged.
- [ ] Run AGENTS.md Node image/AFM, added CSV/IBW/browser tests, adapted Windows Fiji if available (honest skips); verify double injection hash.
- [ ] Update docs, independent complete review, compare project baseline, retain backup manifest, integrate and verify actual product HTML.

### Task 3: GitHub deployment — NM-04
- [ ] Create a brand-new isolated dry-run clone path; confirm INNO-KAIST/INNO-Analytics, main index and source diff, no secrets/private paths or raw data; do not modify any prior clone.
- [ ] On successful preflight use separate brand-new isolated publish path, push authorized main+source, verify both remote HEADs and public Pages content SHA-256 against actual integrated local HTML. Record limitations and report commit links.
