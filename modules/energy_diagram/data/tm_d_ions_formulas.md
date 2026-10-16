# d^n ions in octahedral fields — formulas, matrices and parameters

Companion to `tm_d_ions.json` and the executable reference `dn_ligand_field.py`.

## 1. Model

H = Σ_{i<j} e²/r_ij (Racah A, B, C) + V_oct(Dq), no spin–orbit.

- Racah ↔ Slater: A = F0 − 49F4, B = F2 − 5F4, C = 35F4, where F2 = F²/49 and F4 = F⁴/441.
- Octahedral one-electron matrix in the complex |l=2, m⟩ basis (eg higher):
  ⟨±2|V|±2⟩ = Dq, ⟨±1|V|±1⟩ = −4Dq, ⟨0|V|0⟩ = 6Dq, ⟨2|V|−2⟩ = ⟨−2|V|2⟩ = 5Dq,
  so e_g = +6Dq, t2g = −4Dq, and Δo = 10Dq.
- States are classified by S (from ⟨S²⟩) and by the O-group irrep from the characters of C3[111], C2[001], C4[001] and C2'[110] acting on the orbital part. Accidental degeneracies are decomposed through the character table.
- Tetrahedral: pass Dq_T = −(4/9)Dq_oct (negative Dq). Labels then refer to Td; drop the g subscripts.

References: Y. Tanabe & S. Sugano, J. Phys. Soc. Jpn. 9, 753 (1954) and 9, 766 (1954). A. B. P. Lever, *Inorganic Electronic Spectroscopy*, 2nd ed. (Elsevier, 1984). M. G. Brik & A. M. Srivastava, ECS J. Solid State Sci. Technol. 7, R3079 (2018) (d³ review). The review is cited but was not re-read in this session because the site was blocked.

## 2. Closed forms (each checked against the exact diagonalization to 1e−6 cm⁻¹)

**d³ (Cr³⁺, Mn⁴⁺, V²⁺), octahedral, energies relative to ⁴A2g:**

| state | energy |
|---|---|
| ⁴T2g(F) | 10Dq |
| ⁴T1g(F) | 7.5B + 15Dq − ½·√(225B² − 180B·Dq + 100Dq²) |
| ⁴T1g(P) | 7.5B + 15Dq + ½·√(225B² − 180B·Dq + 100Dq²) |

Inverting these gives a working method: Dq = E(⁴T2)/10, then B comes from E(⁴T1(F)). The ²E, ²T1 and ²T2 levels have no closed form. Their strong-field blocks are 4×4 (²E), 5×5 (²T1) and 5×5 (²T2), plus ²A1 and ²A2 at 1×1 each; these sizes were counted from the exact calculation. C therefore has to come from a numerical fit to E(²E), for example with `derive_d3()` in `pipeline/build_tm.py`.

**d⁸ (Ni²⁺), octahedral, relative to ³A2g:** ³T2g = 10Dq. ³T1g(F) and ³T1g(P) take the same form as the d³ quartets: 7.5B + 15Dq ∓ ½√(225B² − 180B·Dq + 100Dq²). This was verified numerically.

**d⁵ (Mn²⁺, Fe³⁺):** ⁴A1g(G) and ⁴Eg(G) lie at exactly 10B + 5C above ⁶A1g, independent of Dq (verified). Block sizes are ⁶A1 1; ⁴A1 1, ⁴A2 1, ⁴E 2, ⁴T1 3, ⁴T2 3; ²A1 4, ²A2 3, ²E 7, ²T1 8, ²T2 10.

**Free-ion d³ terms relative to ⁴F:** ⁴P = 15B, ²G = 4B + 3C, ²P = ²H = 9B + 3C, ²F = 24B + 3C, ²D pair = 20B + 5C ∓ √(193B² + 8BC + 4C²). The code reproduces all of these.

## 3. Free-ion Racah parameters (from NIST ASD term centroids, computed here)

| ion | B from (⁴P−⁴F)/15 or (³P−³F)/15 | B, C fit (no Trees) | B, C, α fit (Trees) |
|---|---|---|---|
| Cr³⁺ (Cr IV) | 917.2 | 945.5 / 3590.5 (rms 640) | 975.3 / 3476.7 / α 72.3 (rms 208) |
| Mn⁴⁺ (Mn V) | 1065.0 | 1102.1 / 4277.4 (rms 749) | 1137.8 / 4141.1 / α 85.9 (rms 211) |
| Mn²⁺ (Mn III) | from ⁴G = 10B + 5C and ⁴D = 17B + 5C: B = 786.7, C = 3795.7 | 887.2 / 3449.6 (rms 755) | 918.5 / 3255.5 / α 72.7 (rms 161) |
| Co²⁺ (Co III) | 970.8 | 989.7 / 4033.7 (rms 717) | 1023.4 / 3905.2 / α 81.1 (rms 232) |
| Ni²⁺ (Ni III) | 1055.7 | 1086.2 / 3980.0 (rms 572) | 1115.4 / 4009.4 / α 89.5 (rms 10) |

The spin-allowed-separation values agree with the textbook free-ion B values for Cr³⁺ (918) and Co²⁺ (971). Ni²⁺ is often quoted as 1041, and the corresponding textbook Mn²⁺ value is 960. For consistency, use the (⁴P−⁴F)/15 column when computing nephelauxetic ratios β = B/B0. For Mn²⁺ this session could not reproduce the textbook 960 cm⁻¹. The NIST ⁴G/⁴D derivation gives 787, and the fits give 887–919. The spread comes from configuration interaction, so any quoted free-ion B for d⁵ depends on the method.

## 4. Host examples (derived here with the exact matrices)

- **Ruby, Al2O3:Cr³⁺.** Inputs: ⁴T2 18,000 cm⁻¹, ⁴T1 24,400 cm⁻¹, R1 694.3 nm. These are canonical values that were not re-fetched here. Result: Dq = 1800, B = 618.7, C = 3248, C/B = 5.25.
- **K2SiF6:Mn⁴⁺.** Inputs: ⁴T2 band at 455 nm and ²E ZPL at 620.5 nm (Sijbom et al., ECS JSS 5, R3040 (2016)). Result: Dq = 2198. B cannot be fixed without the ⁴T1 band; assuming C/B = 4.5 gives B = 768 and C = 3456.
- **Gaps:** Cr³⁺ in K2NaScF6, LiCaAlF6 and GdAlO3, and Mn²⁺ in ZnF2 and KMgF3, could not be sourced this session.

## 5. Tanabe–Sugano grids

`tm_d_ions.json → tanabe_sugano_grids` holds E/B against Dq/B = 0.1…4.0 for d²–d⁸. All grids use C/B = 4.5. For d⁴–d⁷ the ground-state crossover appears as a switch of the ground term, because energies are measured from the ground state at each Dq/B. To regenerate the grids with another C/B, run `DN(n).levels(B, C, Dq)`.
