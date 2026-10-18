# INNO Analytics — UI language rule (2026-10-17)

The user decided: **the UI skeleton is English; explanatory prose stays Korean.**

## English (translate)
- Tab names (main nav): `XPS · EM Image · Image Analysis · EDS Maps · EDS Spectrum · Diffraction · EDS Interference · Energy Diagram · DRE Kinetics` (2026-10-19: EM Metadata became the **Metadata** view inside EM Image — sub-switch `Image · Metadata`; nav stays one row with short names / `More ▾` — modules/nav/nav_fit.js)
  - EM Image = the image view of the opened electron-microscopy file (Velox EMD, DM3/DM4, Bruker BCF, MRC, TIFF): S/TEM, SEM.
  - Image Analysis = the ImageJ-style general analysis (optical microscopy, AFM, SEM, TEM — any image).
- View / panel / card titles, section headings, sub-tab names.
- Menu items, button labels, toggle/checkbox labels, select options, form field labels, units labels.
- Table column headers, legend and axis titles, empty-state titles, tooltip (`title=`) names of controls.
- Dialog titles and their buttons (Open, Cancel, Apply, Close, Export CSV …).

## Korean (keep)
- Explanatory hints and help sentences, notes/toasts/warnings/errors, status sentences, info-panel prose, the user guide (docs/사용 가이드.md), SPEC/HISTORY.
- A short Korean hint *next to* an English control is fine when it explains something non-obvious.

## Style
- Title Case for tab, panel and section titles (`Peak Fitting`, `Analyze Particles`). Sentence case for buttons and field labels (`Export CSV`, `Fit range`, `Add component`).
- Use the established term of the field / reference software: ImageJ menu names in Image Analysis (`Analyze Particles…`, `Convert to Mask`), Avantage/CasaXPS terms in XPS (`Peak Fit`, `Background`, `Charge Correction`), Velox terms for EM.
- Do not translate: element symbols, units, file-format names, physical-quantity symbols.
- Never change: ids, data-* keys, state/storage keys, saved-file formats and their field names, CSV column names that the app itself re-imports, `.ijm`/macro strings, code identifiers.
- Keep labels short enough for the mobile nav (≤ ~14 characters per tab label; use a short form via CSS/attribute on phones if needed).
