# Energy-level dataset for an interactive energy-level diagram tool

Compiled on 2026-09-29. Energies are in cm⁻¹. Wavelengths are vacuum values (1e7/E) unless a source quotes something else.

## Files

| file | content |
|---|---|
| `ln3_laf3.json` | 4f^N multiplet barycenters of Ln³⁺ (Ce … Yb) in LaF3 up to 52,000 cm⁻¹: 652 multiplets with S, L, J, E, provenance, eigenvector composition, Stark ranges and notes |
| `ln2_and_5d.json` | Divalent Ln (Sm²⁺, Eu²⁺, Tm²⁺, Yb²⁺) and Ce³⁺/Eu²⁺ 4f–5d data: NIST free-ion levels, host data with sources, Dorenbos relations |
| `tm_d_ions.json` + `tm_d_ions_formulas.md` | d^n ions: NIST-derived free-ion B and C, closed-form formulas, host examples, Tanabe–Sugano grids for d²–d⁸ |
| `dn_ligand_field.py` | Exact d^n octahedral ligand-field solver with S and O-irrep classification, built from first principles and validated |
| `fit_free_ion_BC.py`, `free_ion_BC_fit.json` | Fit of B, C (and α) to NIST term centroids |
| `curated_transitions.json` | 63 UCNP/PA-relevant transitions. Each has a typical λ, a LaF3 barycenter λ and a LaF3 Stark window |
| `host_phonons.json` | Maximum phonon energies, with the secondary source and its primary reference |
| `pipeline/` | Scripts that produced the Ln³⁺ data. They contain absolute scratch paths to a qlanth clone; to re-run, `git clone https://github.com/zia-lab/qlanth` and edit the paths in `qlh.py`, `common.py` and `stage*.py` |

## 1. Ln³⁺ in LaF3 (`ln3_laf3.json`)

**Primary source.** W. T. Carnall, G. L. Goodman, K. Rajnak, R. S. Rana, J. Chem. Phys. 90, 3443 (1989), together with the Stark-level tables of ANL-88-8. The tables were used as digitized in the qlanth repository (`data/Carnall.xls`, one sheet per ion, with Exp and Calc columns for every crystal-field level). OSTI (for ANL-88-8), crossref and PMC were not reachable from this sandbox, so the qlanth digitization was the machine-readable route.

**Why a recomputation was needed.** The qlanth sheet fills multiplet labels down each block. Where multiplets overlap, whole blocks carry a single label. Examples: Pr 1I6/3P1, Nd 4F5/2+2H9/2, Nd 4G5/2+2G7/2, Nd 2K13/2/4G9/2, Sm and Dy clusters. Barycenters cannot be taken reliably from the label column alone.

**Method.**
1. Rebuilt Carnall's full effective Hamiltonian from qlanth's Python sparse matrices (`pyql/sparsef1..7.py`, Lizarazo-Ferro et al., PRB 113, 075127 (2026)) using Carnall's Table-1 parameters (`data/carnall.py`). Conventions: M2 = 0.56 M0, M4 = 0.38 M0, P4 = 0.75 P2, P6 = 0.50 P2, spin–spin off, and hole–particle conjugation for N > 7 as in qlanth (flip ζ, T^i, B^k_q; t2Switch = 0).
2. Checked the result against Carnall's own Calc column. After matching every row one-to-one (Hungarian assignment plus a constant offset), the rms differences are:

   | ion | rms (cm⁻¹) | ion | rms (cm⁻¹) |
   |---|---|---|---|
   | Ce | 0.3 | Tb | 14 |
   | Pr | 0.7 | Dy | 12 |
   | Nd | 4 | Ho | 9 |
   | Pm | 4 | Er | 16 |
   | Sm | 9 | Tm | 6 |
   | Eu | 10 | Yb | 70 |
   | Gd | 14 | | |

   Each ion's `meta.validation` entry holds its full statistics. The Yb figure is large because Table 1 of the paper does not include Yb crystal-field parameters, so qlanth reuses the Tm ones.
3. Assigned every crystal-field state to the free-ion (intermediate-coupling) multiplet onto which it projects most strongly. This gives, for each multiplet, its Stark levels and the matching Carnall rows (Exp and Calc).
4. Chose the recommended value `E_cm1` in this order:
   - `exp`: the mean of all observed Stark levels, used only when every Stark level of the multiplet was observed;
   - otherwise `calc`: the mean of Carnall's calculated Stark levels;
   - otherwise (above Carnall's table range, or not tabulated) the qlanth recomputation;
   - for multiplets where crystal-field J-mixing makes the assignment ambiguous (`n_stark_assigned ≠ n_stark_expected`), the intermediate-coupling free-ion energy. This equals the projection-weighted barycenter.

   All intermediate values are kept, so the tool can pick its own convention: `E_exp_mean_observed`, `E_carnall_calc_bary`, `E_calc_bary` and `E_freeion_calc`.
5. Energy zero is Carnall's scale: the lowest Stark level of the ground multiplet sits at about 0. Barycenters are therefore higher than the "lowest Stark level" numbers often quoted in the literature.

**Labels.**
- `term` is Carnall's conventional label when it is consistent with the eigenvector. Otherwise it is the largest seniority-summed SL component.
- `term_unique` disambiguates repeated SLJ labels, for example Er 2H11/2(1) at 19,339 and 2H11/2(2) at 51,363.
- `composition` gives the eigenvector in the seniority-indexed basis (e.g. `2H2` is the second 2H term).
- Notes flag the mixing cases you asked about: Pr 1I6/3P1; Er 2H11/2 and 4S3/2 thermal coupling; Er 2G9/2 (Carnall's label) is the level that UCNP papers call 2H9/2 (≈410 nm); Tm 3F2/3F3 and the 3H4/3F4 naming swap; Nd 4G5/2+2G7/2 and 4F5/2+2H9/2.
- 27 levels have row labels that disagree with the eigenvector. They are flagged in `note`, and most are fill-down or typo artefacts in the sheet (e.g. Eu "7F5" for 7F6, Er "4F7/2" for 4D7/2).

**Coverage.** Levels are included up to 52,000 cm⁻¹. The count of recommended values taken from experiment is in parentheses.

| ion | multiplets (exp) | ion | multiplets (exp) |
|---|---|---|---|
| Ce | 2 (2) | Tb | 83 (4) |
| Pr | 13 (6) | Dy | 89 (18) |
| Nd | 39 (25) | Ho | 60 (5) |
| Pm | 80 (0; calc only) | Er | 34 (23) |
| Sm | 117 (24) | Tm | 12 (3) |
| Eu | 100 (5) | Yb | 2 (1) |
| Gd | 21 (12) | | |

Carnall's tables stop at about 25k (Pm, Eu), 39k (Tb, Dy, Ho), 47–51k (Pr, Sm, Gd, Er), and 68–75k (Nd, Tm). Above those points, values are recomputed.

**Cross-check against the reference values in the brief.** Values are barycenter minus reference.

| reference | difference (cm⁻¹) | note |
|---|---|---|
| Tm 1G4 | +42 | |
| Tm 1D2 | +61 | |
| Tm 3P0 | +104 | |
| Tm 3F4 | +129 | |
| Tm 3H5 | +130 | |
| Tm 3H4 | +145 | |
| Tm 1I6 | +278 | lowest Stark level 34,769 matches 34,800 |
| Er 4I13/2 | +201 | lowest Stark 6,604 |
| Er 4F9/2 | +253 | lowest 15,391 |
| Er 2H11/2 | +239 | lowest 19,266 |
| Er 4F7/2 | +220 | lowest 20,656 |
| Er 4G11/2 | +206 | lowest 26,526 |
| Er 2H9/2 (≈24,500) | — | is Carnall's 2G9/2 at 24,748 (+248) |
| Yb 2F5/2 | +200 | the 10,250 reference is the lowest Stark level (10,260 observed) |
| Eu 5D0 | +13 | |
| Tb 5D4 | +42 | |
| Gd 6P7/2 | +96 | |

Every difference above 200 cm⁻¹ comes from comparing a barycenter with an absorption onset or lowest Stark level. None of them is a data error.

## 2. Divalent ions and 5d levels (`ln2_and_5d.json`)

- **Free-ion levels from the NIST ASD:** Ce IV 5d 2D3/2 = 49,737 cm⁻¹; Eu III lowest 4f⁶5d = 33,856 cm⁻¹; Sm III 7F_J and 5D0 = 14,533, lowest 4f⁵5d = 26,284; Tm III 2F5/2 = 8,774, lowest 4f¹²5d = 22,897; Yb III lowest 4f¹³5d = 33,386.
- **Dorenbos relations** (two secondary sources that quote them): D(2+) = 0.64·D(3+) − 0.233 eV; E(Eu²⁺) = 0.64·E(Ce³⁺) + 0.53 eV; ε_CFS(Eu) = 0.77·ε_CFS(Ce). The free-ion Eu²⁺ reference is 34,000 cm⁻¹.
- **Host data:**
  - CaF2:Eu²⁺: absorption 3.06 eV, emission 2.92 eV (≈425 nm).
  - LaF3:Eu²⁺: anomalous emission at 600 nm.
  - α-NaYF4:Ce³⁺: absorption 300 nm, emission 332 nm.
  - Ce:LaF3: laser at 286 nm, tunable 275–315 nm.
  - CeF3: 5d features at 5.01–6.45 eV.
  - SrF2:Sm²⁺: 5D0→7F0 at 683.8 nm and 5D0→7F1 at 696.7 nm. Assignments differ between two papers; see the note in the file.
  - CaF2:Sm²⁺: d–f emission at 710 nm.
  - Tm²⁺ in MF2: ground-state splittings and the 1120/1180 nm lines.
  - Yb²⁺ bands: 360 and 310 nm.
- **Gaps:**
  - Sm²⁺ crystal-field 7F_J/5D_J levels in SrF2, CaF2 and BaFCl (Wood & Kaiser 1962; Axe & Sorokin 1963) were not retrieved.
  - The first 5d band maximum of Ce³⁺ in LaF3 (≈250 nm, 40,000 cm⁻¹) is not verified from a primary source.
  - Ce³⁺ in β-NaYF4 and Eu²⁺ in NaYF4 were not found.

## 3. d ions (`tm_d_ions.json`, `tm_d_ions_formulas.md`)

- **Solver.** The exact solver reproduces the textbook free-ion d³ term formulas, the closed-form d³ and d⁸ quartet and triplet formulas, and the Dq-independent 10B+5C level of d⁵.
- **Free-ion B from NIST.** Using (⁴P−⁴F)/15:

  | ion | B (cm⁻¹) | textbook value |
  |---|---|---|
  | Cr³⁺ | 917 | 918 |
  | Co²⁺ | 971 | 971 |
  | Mn⁴⁺ | 1065 | – |
  | Ni²⁺ | 1056 | 1041 |

  Full B, C (and Trees α) fits are also given.
- **Host examples derived here:**
  - Ruby: Dq = 1800, B = 619, C = 3248, from the canonical band positions, which were not re-fetched.
  - K2SiF6:Mn⁴⁺: Dq = 2198, from 455 nm and the 620.5 nm ZPL (Sijbom et al. 2016). B needs the ⁴T1 band, which was not sourced.
- **Gaps:** host Dq, B and C for Cr³⁺ fluorides (K2NaScF6, LiCaAlF6), GdAlO3:Cr³⁺, Mn²⁺ fluorides (ZnF2, KMgF3), and Co²⁺ and Ni²⁺ hosts.

## 4. Curated transitions (`curated_transitions.json`)

Each transition has three wavelengths:
- `typical_nm`: from the brief, or a commonly reported value that is flagged as not individually source-verified;
- `lambda_bary_nm_LaF3`: computed from the LaF3 barycenters;
- `lambda_window_nm_LaF3`: the span of all Stark-to-Stark combinations.

Of the 59 transitions that have a typical value, 58 fall inside their LaF3 window (±1%). The exception is Gd 6I7/2 at 275 nm, where "275 nm" refers to the whole 6I_J group (272–278 nm).

## 5. Phonon energies (`host_phonons.json`)

| host | cm⁻¹ | host | cm⁻¹ |
|---|---|---|---|
| β-NaYF4 | 360 (another source gives 370) | Y2O3 | 600 |
| LaF3 | 350 | YAG | 860 |
| LiYF4 | 460 | YVO4 | 950 |
| CaF2 | 320 (Raman) | YBO3 | 1050 |
| SrF2 | 284 | silicate glass | ~1000 |
| BaF2 | 240 | germanate glass | ~900 |
| PbF2 | 257 | tellurite glass | ~750 |
| | | bismuthate glass | ~440 |
| | | Ga-La-S glass | ~450 |

Each value carries its citation. Phosphate glass, NaGdF4, NaLuF4 and KY3F10 were **not sourced** and are left null.

## Caveats

- The divalent-ion, d-ion and phonon numbers taken from the web were read through an automated page reader (WebFetch). Primary tables should be spot-checked before publication.
- NIST ASD numbers are verbatim table rows but came through the same reader.
- Values recalled from standard knowledge are explicitly labelled as such: ruby band positions, Dorenbos 6.12 eV, and the Layne phosphate value.
