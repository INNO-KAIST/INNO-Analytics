# Diffraction — 동역학 스폿 세기 (`src/diff_dyn.js`, `DIFF_DYN`) · 2026-10-17

정대축 SAED 의 스폿 세기를 **Bloch 파 동역학 이론 (0 차 Laue 영역, ZOLZ)** 으로 계산한다. 운동학 |F|² 가 못 그리는
이중 회절(운동학 금지 자리가 밝아짐), Gjønnes–Moodie(GM) 소멸, 두께에 따른 Pendellösung 진동을 상 찾기 점수에 쓰기 위한 것.
DOM 없음 · Worker 가능 · 큰 자료 없음 (27 kB). 불러오는 순서: `diff_sg.js → diff_scat.js → diff_core.js → diff_dyn.js`.

## 1. API

```js
const R = DIFF_DYN.zoneIntensities(phase, [u,v,w], {
  kV: 200,                 // 가속 전압
  thickness: [10,20,…,100],// nm, 수 하나 또는 배열
  gmax: 12,                // 출력 범위 (1/nm) — 관측 스폿 범위에 맞춘다
  maxBeams: 200,           // 강한 빔(정확히 푸는 빔) 상한. 같은 |g| 껍질은 자르지 않는다 (대칭 보존)
  gStrong: max(gmax,35),   // 강한 빔 반지름 상한 (1/nm)
  gBethe: 50, maxBethe: 3000, // Bethe 섭동 빔 (강한 빔 밖 ~ 50 1/nm). 0 이면 끔
  tilt: [tx,ty],           // mrad, 입사 빔 접선 성분 (2D 좌표 e1·e2 방향). 기본 정대축
  absorption: 0,           // κ: 비례 흡수 U'_g = κ U_g. > 0 이면 행렬 지수로 정확히 (느림, 0.125 nm 걸음)
  weak: true,              // g_cut < g ≤ gmax 출력 빔(큰 세포에서 강한 빔이 gmax 에 못 미칠 때) = 1 차 섭동
  Bdefault: 0.5,           // B = 0 (CIF 에 열진동 없음) 원자에 넣는 Debye–Waller B (Å²)
  dw: true, Bset: null,    // DW 끄기 / 모든 원자 B 고정 (검증용)
  gm: true, gmProbe: [17,43,89], // 정대축 GM 판정 두께 (nm)
  beams: null,             // [[h,k,l]…] 강한 빔 직접 지정 (시험용; 이때 Bethe·weak 기본 끔)
  origin: 'auto', ref: null, includeHOLZ: false (미구현 → meta.warn)
});
// R = { zone, axis, e1, e2, thickness, direct:{I:[…]},
//       beams:[{ h, g, x, y,            // 2D 좌표 (DC.zonePattern 과 같은 e1·e2)
//                I:[두께마다],          // 입사 = 1 대비 세기
//                Ikin, IkinRel,         // 운동학 |F|² (Å², DW 포함) · 출력 범위 최대 대비
//                absent, kinForbidden,  // 공간군 소멸 / (소멸 또는 |F|² < 1e-6 최대)
//                ddRoute:[[g1,g2]…], nRoute, // 운동학 금지 자리에 이르는 허용 반사 쌍 (센 것 6 개, 전체 수)
//                gmExtinct,             // 경로가 있는데 정대축에서 모든 두께 ~0 (GM 소멸)
//                dynZero,               // 정대축에서 동역학적으로도 0 (GM 이거나 ZOLZ 경로 없음)
//                pert, s, U:[re,im] (nm⁻²), xi (nm) }],
//       meta:{ lambda (nm), gamma, K, U0, nBeams, nBethe, nWeak, gCut, solver:'real'|'complex'|'expm', dim, orthoErr, ms, msSolve, warn } }

DIFF_DYN.thicknessAverage(R, [t1,t2])   // → {I:[빔마다], I0, n, tRange}  사다리꼴 평균
DIFF_DYN.zoneAverage(phase, uvw, {tRange:[10,80], tStep:2, tiltAvg:0, nTilt:6, …zoneIntensities 옵션})
                                        // 두께 모를 때: beams[i].Iavg (입사 대비), .Irel (회절 빔 최대 = 1). 상·조건별 캐시
DIFF_DYN.beamOf(R, [h,k,l])             // hkl 로 빔 찾기 (중심화 소멸 자리는 빔이 아니므로 null)
DIFF_DYN.wavelength(kV) (Å), DIFF_DYN.gammaRel(kV)
```

**두께를 모를 때의 선택**: 10–80 nm 균일 평균(2 nm 간격). SAED 영역 선택 조리개 안에는 두께 분포(쐐기·입자)가 들어가고, 시료 두께는
대개 20–100 nm 이다. 평균은 Pendellösung 진동을 지우고 "어느 두께에서나 센 / 약한 / 이중 회절로 생기는" 스폿의 차이만 남긴다.
얇은 나노입자(< 20 nm)는 `tRange:[5,30]`. 미세 기울임 평균(`tiltAvg` 1–2 mrad)은 GM 자리가 얼마나 쉽게 생기는지를 보여 주지만 7 배 느리다.

## 2. 식 (단위: 길이 nm, 역공간 1/nm, 2π 없음)

| 양 | 식 | 출처 |
|---|---|---|
| 파장·γ | λ = hc / √(E(2m₀c² + E)), γ = 1 + E/m₀c² (200 kV: 0.0250793 Å, 1.391391) | CODATA 2018 |
| 퍼텐셜 계수 | V_g = (h²/2πm₀e) · F_g* / Ω = 47.878 V·Å² · F_g*/Ω | Kirkland 2010 §5; De Graef 2003 §2.4 |
| U_g | U_g = (2m₀γe/h²) V_g = **γ F_g* / (π Ω)** [Å⁻²], F_g = Σ occ·f_e(s) e^{−Bs²} e^{+2πi g·x} (f_e: IT C 4.3.2.2, Å) | De Graef 2003 §2.4; Spence & Zuo 1992 |
| 소멸 거리 | ξ_g = K / \|U_g\| (= 1/(λ\|U_g\|)) | Hirsch et al. 1965 |
| 구조 행렬 | γʲ Cʲ = H Cʲ, H_gg = s_g, H_gh = U_{g−h}/2K, K² = K₀² + U₀ | De Graef 2003 §5.5 (고에너지 근사, 표면 ⟂ 빔) |
| 들뜸 오차 | s_g = −(2 k_t·g + g²)/(2K) (정대축 −g²/2K) | |
| 파동 | ψ_g(t) = Σ_j C₀ʲ* C_gʲ exp(2πiγʲ t), I_g = \|ψ_g\|² | |
| Bethe | H_gh ← H_gh + Σ_w H_gw H_wh / (0 − s_w)  (강한 빔 밖 g ≤ 50 1/nm) | Bethe 1928; Zuo & Weickenmeier, Ultramicroscopy 57 (1995) 375 |
| 바깥 빔 | ψ_w = Σ_j a_j W_j (e^{2πiγʲt} − e^{2πis_w t})/(γʲ − s_w), W_j = Σ_h H_wh C_hʲ (1 차 섭동) | |
| 흡수 | U → (1+iκ)U, U₀ → U₀ + iκU₀; ψ(t+δ) = exp(2πiδH)ψ(t) (축척·제곱 + Taylor 12 차) | Humphreys & Hirsch 1968 (비례 모형); Higham 2005 |

- 관례 검증: F 의 켤레(U_g ∝ F_g*) 와 전파 부호는 비중심대칭 ZnO [110] 에서 다중 슬라이스와 Friedel 짝 세기 I(002)/I(00-2) 까지 0.2 % 안에서 같다 (부호가 틀리면 두 값이 뒤바뀜).
- 고유값: 실대칭 Householder 3 대각화 + 암시적 QL (Wilkinson–Reinsch, JAMA tred2/tql2). 반전 중심이 있으면 원점을 그리로 옮겨 U 가 실수 → N 실대칭.
  없으면 Hermite H 를 2N 실대칭 [[Re, −Im],[Im, Re]] 로 풀고 v = x + iy 를 가중 ½ 로 모두 더한다 (각 복소 방향이 v·iv 로 두 번 나오고, 축퇴 부분공간에서도 Σvv† = 2P 라 고르기 불필요).
- 빔: [uvw] 평면의 모든 정수 hkl. **격자 중심화(R = I, t ≠ 0) 소멸 자리만 뺀다** — 허용 g 들은 덧셈에 닫힌 부분격자라 어떤 다중 산란으로도 못 간다.
  미끄럼·나사·구조상 0 인 자리(Si 200, 222 …)는 넣는다 (이중 회절로 들뜰 수 있음).
- GM 판정은 대칭식이 아니라 수치: 껍질 단위로 고른 대칭 빔 집합에서 세 두께(17·43·89 nm) 모두 I < 1e-9 × 최대 회절 빔이면 `dynZero`, 그중 ZOLZ 경로가 있으면 `gmExtinct`.
- 1 차 섭동인 흡수 계산(γ → γ + iq)은 쓰지 않는다: 비례 모형에서 H' 의 Bloch 기저 비대각 성분 κ·C†DC 가 고유값 간격과 비슷해 Si [011] κ = 0.1 에서 48 % 틀렸다 (scipy expm 대조). 그래서 행렬 지수 걸음.

## 3. 검증 (`node test/test_dyn.js` 79 통과 · 기록 `test/test_dyn.out.txt` · 기준 `python3 -I ref/gen_dyn_ref.py` → `test/ref_dyn.json`)

| 시험 | 결과 |
|---|---|
| 고유값 풀이 80×80 | 잔차 1.2e-14, 직교 9e-15 |
| 두 빔 Pendellösung, Si 220 s = 0 · s ≠ 0 (w = −3.4) · 흡수 κ 0.08 해석해 ½e^{−4πq₀t}[cosh(2πκt/ξ) − cos(2πt/ξ)] | 최대 차 3e-15 |
| ξ_g ↔ Hirsch et al. 표 (100 kV, B = 0; Y. Liao, *Practical Electron Microscopy and Database* 표 4148 재인용): Al 111·200, Si 111·220·400, Au 111·200·220, Fe 110·200, Cu 111·200, MgO 200·220 | 최대 5.2 % (Al 111 55.6 / 56.3 nm). Ni·Ag·Pt 는 옛 산란 인자 표와 8–11 % 차라 제외 |
| Si 400, 200 kV (B 0.46) | 168 nm ↔ 측정 194 ± 30 nm (Kossel–Möllenstedt, Acta Phys. Sin. 71 (2022), doi:10.7498/aps.71.20212271) |
| 얇은 결정 t = 0.02 nm (Bethe 없이) | I_g = (πt/ξ_g)² sinc²(πs_g t) 와 0.5 % 안, 상대 세기 = |F|² 비 0.4 % 안 (Si·Au·Y₂O₃·ZnO) |
| 독립 Bloch: gemmi(CIF 원자·c4322) + 따로 센 빔 목록 + scipy `expm` (고유값 없음, 흡수 정확) — 13 경우 (Si·Au·Fe·ZnO·루틸·β-NaYF₄·Y₂O₃, 100–300 kV, κ 0–0.1, 기울임) | 최대 차/최대 ≤ 1.3e-6 (ZnO·β-NaYF₄ 는 CIF 좌표 0.3333 처리 차이로 ≤ 1.1e-3). 자동 빔 수 = 독립 목록 수 |
| 원점 불변: Si [011] 실수 N 경로 ↔ CIF 원점 복소 2N 경로 | 6e-14 |
| 흡수 없음 Σ I | 1 ± 1e-12 (Bethe 포함 Hermite 유지) |
| **다중 슬라이스 (같은 ZOLZ 모형, 독립 풀이)**: abTEM 투영 퍼텐셜(B 0.5) → c/8 슬라이스 numpy FFT, t = 10·30·60 nm | 회절 빔 최대 차/최대: Si [011] 0.09 %, Si [001] 0.11 %, ZnO [110] 0.15 %, 루틸 [001] 0.39 %, α-Fe [111] 2.2 %, Au [001] 5.7 % → 2.3 % (Δx 0.04 → 0.025 Å; 다중 슬라이스 쪽 이산화 오차). 루틸 5.2 mrad 기울임 0.1 % |
| **3D 다중 슬라이스 (abTEM 기본, ≈1 Å 슬라이스, HOLZ 포함)** ↔ ZOLZ Bloch | Si 0.5–1.2 %, ZnO 3.5 %, 루틸 5 %, α-Fe 22 %, Au [001] 76 % (t = 60 nm) — 무거운 원소·짧은 반복 길이에서 HOLZ 결합이 ZOLZ 세기를 바꾼다 |
| 순위 상관: 3D 다중 슬라이스(3 두께 평균)를 '실측' 으로 볼 때 Spearman ρ (동역학 ZOLZ · 운동학 |F|²) | Si [011] 0.991 · 0.686, Si [001] 0.976 · 0.904, ZnO [110] 0.999 · 0.847, 루틸 [001] 0.988 · 0.790 (센 20 개: 0.883 · **−0.785**), Au [001] 0.865 · 0.856, α-Fe [111] 0.942 · 0.944 |
| 이중 회절·GM | Si [011] (200): 소멸인데 I = 1.5e-3 / 2.3e-2 / 4.5e-2 / 0.12 (10/20/30/60 nm), 경로 (11-1)+(1-11) 등. Si [011] (22-2) 구조상 0 → 1.8e-3. Si [001] (200): ZOLZ 경로 없음 (허용 hk0 은 h+k ≡ 0 mod 4) → 6e-30. 루틸 P4₂/mnm [001] (100)(010)(300)(030), ZnO P6₃mc [110] (001)(003), ZrO₂ P2₁/c [100] (001)(010)(003)(030), Y₂O₃ [110] (1-10): 경로가 있어도 GM 소멸 (≤ 1e-26). 루틸 5.2 mrad 기울임 → (100) 2.3e-2 (20 nm); 다중 슬라이스도 정대축 1.9e-7 (격자 표본화 잡음) ↔ 기울임 1.2e-2. Fe₃O₄ [110] (002) 소멸이지만 4.5e-2 |
| 수렴 (두께 10–80 nm 평균, 기본값 ↔ 강한 빔 ≤ 45 1/nm + Bethe 60) | Si [011] 0.02 %, Au [001] 0.06 %, Au [011] 0.9 %, α-Fe [111] 0.02 %, LiYF₄ [100] 0.85 %. Y₂O₃ [001] (a = 1.06 nm) 2.2 %, β-NaYF₄ [001] 3.4 % (단일 두께 최대 차는 10–20 %) |
| 바깥 빔 1 차 섭동 (Si [011] 12 < g ≤ 16) | 0.07 % |

Bethe 없이 빔만 자르면 수렴이 매우 느리다 (B = 0.5 에서 Au [001] 강한 빔 ≤ 12 1/nm: 680 %, ≤ 30: 36 %; B = 0 이면 더 나쁨 — 무거운 원자의 1s 채널링 상태가 높은 g 까지 필요).
B = 0 은 Au·Fe 가 8 Å⁻¹ 밖으로 수 % 를 흩어 Bloch·다중 슬라이스 모두 수렴하지 않으므로 기본 `Bdefault` 0.5 Å² (상온 무기물 전형값)를 넣는다.

## 4. 속도 (Node 22, 컨테이너 2 코어)

| 경우 | 시간 |
|---|---|
| Si [011] 120 빔, 실수, 10 두께 | 10 ms (풀이 6 ms) · 같은 계산 복소 2N 경로 38 ms |
| 내장 21 상 × 낮은 정대축 2 개, 기본값 (강한 ≤ 200 + Bethe ≤ 3000, GM 판정 포함) | 중앙값 85 ms, 최악 0.42 s (Y₂O₃ [100], Bethe 3000). 비중심대칭 β-NaYF₄ [001] ≈ 0.4–0.5 s |
| `zoneAverage` (36 두께) | 위와 같음 (두께는 거의 공짜) · `tiltAvg` 는 × 7 (0.2–2 s) |
| 흡수 κ > 0 (행렬 지수, N 161, 800 걸음) | 0.26 s |

## 5. 한계
- **ZOLZ 만**: HOLZ 결합 무시. 3D 다중 슬라이스와 Si 1 %, ZnO·루틸 5 %, α-Fe 20 %, Au 최대 76 % 차 (개별 빔). 순위 상관은 여전히 운동학 이상. Si [001] (200) 처럼 ZOLZ 경로가 없는 자리는 HOLZ 경유로만 약하게 생긴다.
- 고에너지·정상 입사 근사 (표면 ⟂ 빔, 후방 산란·곡률 무시), 평면파 (CBED 수렴각 없음), 굽힘·변형·모자이크·적층 결함·비탄성(TDS 흡수 계수) 없음 — 흡수는 비례 모형만.
- 부분 점유 = 평균 퍼텐셜 (가상 결정). 원자 B 없으면 0.5 Å² 가정.
- Bethe 는 t ≫ 1/|s_w| (≈ 1 nm) 에서만 맞다 — t → 0 극한(운동학)은 Bethe 를 끄면 정확.
- 큰 세포(Y₂O₃·β-NaYF₄, a ≈ 0.6–1.1 nm) 는 기본값에서 두께 평균 2–3 %, 단일 두께 10–20 % 정도. `maxBeams: 300` 이면 1 % (×2–3 느림).
- GM 판정은 대칭 빔 집합에 기대므로 `beams` 로 비대칭 목록을 주면 무의미. 정대축에서 1–2 mrad 만 기울어도 GM 자리는 회절 최대의 3–15 % 까지 생긴다 (루틸 (100): 1 mrad 평균 2.6 %, 2 mrad 8 %, 4 mrad 17 %; ZnO (001): 4.9 / 16 / 42 %).

## 6. 상 찾기에 쓰는 법 (제안 — `DI.rankSpots` / `verifyZone` 통합은 주 세션)

```js
// rankSpots(fit, phases, opt) 안, 기존 기하 FOM·운동학 점수로 정렬한 다음:
const nDyn = opt.nDyn ?? 12;                       // 동역학은 상위 후보만 (후보당 0.05–0.5 s, 캐시)
for (const c of C.slice(0, nDyn)) {
  if (!c.hasF) continue;
  const V = c.verify;                                // verifyZone 결과 (rows: 관측 스폿 h, missing …)
  const gObs = Math.max(...V.rows.map(r => 1 / r.d)) * 1.05;
  const D = DIFF_DYN.zoneAverage(c.phase, c.zone, { kV: opt.kV || 200, gmax: gObs,
              tRange: opt.tRange || [10, 80], tiltAvg: opt.tiltAvg || 0 });
  // 1) 관측 세기 ↔ 동역학 평균 세기
  const Io = V.rows.map(r => Math.max(0, r.spot.I)), tot = Io.reduce((a, b) => a + b, 0) || 1;
  const Id = V.rows.map(r => { const b = DIFF_DYN.beamOf(D, r.h); return b ? b.Irel : 0; });
  const rhoDyn = DI.spearman(Io, Id);
  // 2) Forb: '동역학적으로도 0' 인 자리의 관측 세기 몫만 벌점
  //    · 중심화 소멸(빔 없음)·ZOLZ 경로 없음 → 1.0 (HOLZ 경유 아니면 불가능)
  //    · GM 소멸 → 0.3 (1–2 mrad 기울임으로 수 % 생김)
  //    · 이중 회절로 생길 수 있는 운동학 금지 자리 → 0 (예측 세기 Id 가 ρ 에 이미 반영)
  let forbDyn = 0;
  V.rows.forEach((r, i) => { const b = DIFF_DYN.beamOf(D, r.h);
    const w = !b ? 1 : (b.dynZero && !b.gmExtinct) ? 1 : b.gmExtinct ? 0.3 : 0; forbDyn += w * Io[i]; });
  forbDyn /= tot;
  // 3) Miss: 관측 범위(영상 안·가리지 않은 곳) 안 예측 Irel ≥ 0.02 인데 관측 없음 (Irel 가중)
  const got = new Set(V.rows.map(r => r.h.join(",")));
  let zI = 0, oI = 0;
  for (const b of D.beams) { if (b.Irel < 0.02) continue; if (opt.isMasked && opt.isMasked(b.h)) continue;
    if (got.has(b.h.join(","))) oI += b.Irel; else zI += b.Irel; }
  const missDyn = zI + oI > 0 ? zI / (zI + oI) : 0;
  c.dyn = { rho: rhoDyn, forb: forbDyn, miss: missDyn, nBeams: D.meta.nBeams, ms: D.meta.avg.ms };
  c.scoreKin = c.score;
  c.score = c.fom + 5 * forbDyn + (1 - rhoDyn) + 3 * missDyn;   // 같은 가중, 운동학 대신 동역학 항
}
// 정렬: 상위 nDyn 은 동역학 점수, 나머지는 운동학 점수 그대로 (nDyn 바깥이 위로 올라오지 않게 nDyn 안끼리 다시 정렬)
const head = C.slice(0, nDyn).sort((a, b) => a.score - b.score); C.splice(0, nDyn, ...head);
```

- `verifyZone` 의 행 표시: 운동학 금지인데 `ddRoute` 가 있으면 "이중 회절 (g1+g2)", `gmExtinct` 면 "GM 소멸 — 정대축이면 0, 기울면 생김", `dynZero && !gmExtinct` 면 "ZOLZ 로 불가 (HOLZ 또는 다른 상)".
- **두께**: 기본은 평균. 화면에서 "두께 추정" 을 보여 주려면 창 [10,30]·[20,50]·[40,80] 별 ρ 를 계산해 가장 높은 창을 표시만 한다 (후보마다 자유도가 늘어 순위에 쓰면 과적합).
- **미세 기울임**: 관측에 GM 자리가 뚜렷하면 `tiltAvg: 1–2` 로 그 후보만 다시 계산 (× 7 느림).
- 운동학 ρ 는 상위 후보 바깥(동역학 계산 안 함)·원자 없는 상에 그대로 쓴다. Worker 주입 목록에 `diff_dyn.js` 를 `diff_core.js` 다음에 넣는다.

## 7. 파일
- `src/diff_dyn.js` — 본체 (DIFF_DYN). `diff_core.js` 는 고치지 않음 (DC.cell·gvec·rvec·planeBasis·reduceInt·expandAtoms·fElectron·absent·vec 만 사용).
- `ref/gen_dyn_ref.py` — 독립 Bloch (gemmi + scipy expm) · 다중 슬라이스 (abTEM 1.0.10 퍼텐셜 + numpy FFT ZOLZ 전파 / abTEM 3D). `pip install --break-system-packages abtem` (ase·numba 포함).
  abTEM `projection='infinite'` 는 Si [011] 초격자에서 평균 퍼텐셜이 13 % 작게 나와 `'finite'` 를 쓴다.
- `test/test_dyn.js`, `test/ref_dyn.json` (기준 자료 207 kB), `test/test_dyn.out.txt` (마지막 실행 기록).

참고문헌: P. Hirsch, A. Howie, R. Nicholson, D. Pashley, M. Whelan, *Electron Microscopy of Thin Crystals* (1965) · J. Gjønnes, A. F. Moodie, Acta Cryst. 19, 65 (1965) ·
C. J. Humphreys, P. B. Hirsch, Phil. Mag. 18, 115 (1968) · J. C. H. Spence, J. M. Zuo, *Electron Microdiffraction* (1992) · J. M. Zuo, A. L. Weickenmeier, Ultramicroscopy 57, 375 (1995) ·
M. De Graef, *Introduction to Conventional TEM* (2003) · E. J. Kirkland, *Advanced Computing in Electron Microscopy*, 2nd ed. (2010) · *International Tables for Crystallography* Vol. C, Table 4.3.2.2 (Peng et al. 1996) ·
J. Madsen, T. Susi, abTEM, Open Research Europe 1, 24 (2021) · N. J. Higham, SIAM J. Matrix Anal. Appl. 26, 1179 (2005).
