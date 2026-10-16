# 내장 결정상 팩 출처표 (DIFF_PACK 2026-10-17)

`ref/gen_pack.py` 가 이 폴더의 CIF 에서 `src/diff_pack.js` 와 함께 생성한다. CIF 는 COD 에서 받은 원문 그대로이며, 생성기가 SHA-256 앞 16 자리를 확인한다.
격자상수는 nm (CIF Å ÷ 10). B(Å²) 는 `_atom_site_U_iso_or_equiv`×8π² 또는 `_atom_site_B_iso_or_equiv`, 둘 다 없으면 0 (비등방 U 만 있는 CIF 는 B = 0).

| id | 이름 | 화학식 | 공간군 (CIF 설정) | 격자상수 (nm) | COD | DOI | 인용 | 고른 이유 |
|---|---|---|---|---|---|---|---|---|
| pack:Au | Gold | Au (Z = 4) | F m -3 m (#225) | a = 0.407825 | [9008463](https://www.crystallography.net/cod/9008463.html) | — (CIF 에 없음) | Wyckoff, R. W. G. (1963) Crystal Structures 1, 7-83 | Wyckoff 1963 (AMCSD) 표준 참조 구조, a = 4.07825 Å (현대 값 4.0782 Å 와 일치) |
| pack:Al | Aluminium | Al (Z = 4) | F m -3 m (#225) | a = 0.404958 | [9008460](https://www.crystallography.net/cod/9008460.html) | — (CIF 에 없음) | Wyckoff R W G (1963) Crystal Structures 1, 7-83 | Wyckoff 1963 (AMCSD) 표준 참조 구조, a = 4.04958 Å |
| pack:Si | Silicon | Si (Z = 8) | F d -3 m :1 (#227) | a = 0.543094 | [9011998](https://www.crystallography.net/cod/9011998.html) | 10.1107/S0021889875010965 | Hom, T., Kiszenick, W., Post, B. (1975) Journal of Applied Crystallography 8, 457-458 | 298.15 K 다중반사 정밀 격자상수 측정(a = 5.430941 Å), DOI 있음 |
| pack:Cu | Copper | Cu (Z = 4) | F m -3 m (#225) | a = 0.361496 | [9008468](https://www.crystallography.net/cod/9008468.html) | — (CIF 에 없음) | Wyckoff, R. W. G. (1963) Crystal Structures 1, 7-83 | Wyckoff 1963 (AMCSD) 표준 참조 구조, a = 3.61496 Å |
| pack:Ni | Nickel | Ni (Z = 4) | F m -3 m (#225) | a = 0.352387 | [9008476](https://www.crystallography.net/cod/9008476.html) | — (CIF 에 없음) | Wyckoff, R. W. G. (1963) Crystal Structures 1, 7-83 | Wyckoff 1963 (AMCSD) 표준 참조 구조, a = 3.52387 Å |
| pack:Pt | Platinum | Pt (Z = 4) | F m -3 m (#225) | a = 0.39231 | [9008480](https://www.crystallography.net/cod/9008480.html) | — (CIF 에 없음) | Wyckoff, R. W. G. (1963) Crystal Structures 1, 7-83 | Wyckoff 1963 (AMCSD) 표준 참조 구조, a = 3.9231 Å |
| pack:Ag | Silver | Ag (Z = 4) | F m -3 m (#225) | a = 0.40862 | [9008459](https://www.crystallography.net/cod/9008459.html) | — (CIF 에 없음) | Wyckoff, R. W. G. (1963) Crystal Structures 1, 7-83 | Wyckoff 1963 (AMCSD) 표준 참조 구조, a = 4.0862 Å |
| pack:aNaYF4 | α-NaYF4 (cubic) | NaYF4 (Z = 2) | F m -3 m (#225) | a = 0.547 | [1517676](https://www.crystallography.net/cod/1517676.html) | — (CIF 에 없음) | Roy, D.M., Roy, R. (1964) Journal of the Electrochemical Society 111, 421-429 | COD 의 화학량론 α-NaYF4 유일 항목(Fm-3m, Na/Y 4a 자리 0.5/0.5 무질서, F 8c). 1517677 은 Na0.39Y0.61F2.2 비화학량론 결함 구조라 제외 |
| pack:bNaYF4 | β-NaYF4 (hexagonal) | Na1.5Y1.5F6 (Z = 1) | P -6 (#174) | a = 0.59148, c = 0.3496 | [1517672](https://www.crystallography.net/cod/1517672.html) | — (CIF 에 없음) | Grzechnik, A., Bouvier, P., Mezouar, M. et al. (2002) Journal of Solid State Chemistry 165, 159-164 | 고압 연구(Grzechnik 2002)의 상압 P-6 구조 — 1a Y, 1f Y/Na 0.5/0.5, 2h Na 0.5 (CIF 점유율 그대로). 같은 논문의 1517673–75 는 고압 |
| pack:CaF2 | Fluorite CaF2 | CaF2 (Z = 4) | F m -3 m (#225) | a = 0.5462 | [1000043](https://www.crystallography.net/cod/1000043.html) | 10.1088/0022-3719/4/18/016 | Cheetham, A K, Fender, B E F, Cooper, M J (1971) Journal of Physics C 4, 3107-3121 | 중성자 회절 정밀화(R = 0.008), 화학량론 CaF2. 단 a = 5.462(30) Å 로 불확도가 크다(흔한 값 5.463 Å 와 0.02 % 차) |
| pack:LiYF4 | LiYF4 (scheelite) | LiYF4 (Z = 4) | I 41/a :2 (#88) | a = 0.5164, c = 1.0741 | [2001129](https://www.crystallography.net/cod/2001129.html) | 10.1107/S0108270193005876 | Garcia, E., Ryan, R. R. (1993) Acta Crystallographica Section C 49, 2053-2054 | 상온 단결정 정밀화(Acta C, 레이저 호스트), I41/a 원점 2 |
| pack:Y2O3 | Y2O3 (bixbyite, C-type) | Y2O3 (Z = 16) | I a -3 (#206) | a = 1.06056 | [1009013](https://www.crystallography.net/cod/1009013.html) | 10.4028/www.scientific.net/MSF.278-281.680 | Baldinozzi, G., Berar, J.-F., Calvarin, G. (1998) Materials Science Forum 278-281, 680-685 | COD 의 순수 Y2O3 Ia-3 중 위치 매개변수 정밀도가 가장 높음(R = 0.011). 1513300(a = 10.604 Å, 295 K)은 Uovl = 0.00058 Å² 로 비정상적이고 좌표 불확도가 큼 |
| pack:aFe | α-Fe (ferrite, bcc) | Fe (Z = 2) | I m -3 m (#229) | a = 0.28665 | [9008536](https://www.crystallography.net/cod/9008536.html) | — (CIF 에 없음) | Wyckoff, R. W. G. (1963) Crystal Structures 1, 7-83 | Wyckoff 1963 (AMCSD) 표준 참조 구조, 298 K, a = 2.8665 Å |
| pack:Fe3O4 | Magnetite Fe3O4 | Fe3O4 (Z = 8) | F d -3 m :2 (#227) | a = 0.83969 | [9007706](https://www.crystallography.net/cod/9007706.html) | 10.1107/S0108270184008489 | Fleet, M. E. (1984) Acta Crystallographica, Section C 40, 1491-1493 | 단결정 상온 정밀화(Fleet 1984, 시료 633 Fe3.005O4), a = 8.3969 Å 가 표준 값과 일치. 1981 년 항목(9007644)은 a = 8.3941 Å |
| pack:TiO2-rutile | Rutile TiO2 | TiO2 (Z = 2) | P 42/m n m (#136) | a = 0.45937, c = 0.29587 | [9015662](https://www.crystallography.net/cod/9015662.html) | 10.1107/S010876819100335X | Howard, C. J., Sabine, T. M., Dickson, F. (1991) Acta Crystallographica, Section B 47, 462-468 | 중성자 분말 상온 정밀화(Howard 1991), 루틸·아나타제 같은 논문 |
| pack:TiO2-anatase | Anatase TiO2 | TiO2 (Z = 4) | I 41/a m d :1 (#141) | a = 0.37845, c = 0.95143 | [9015929](https://www.crystallography.net/cod/9015929.html) | 10.1107/S010876819100335X | Howard, C. J., Sabine, T. M., Dickson, F. (1991) Acta Crystallographica, Section B 47, 462-468 | 중성자 분말 상온 정밀화(Howard 1991) |
| pack:ZnO | Zincite ZnO (wurtzite) | ZnO (Z = 2) | P 63 m c (#186) | a = 0.32494, c = 0.52054 | [2300112](https://www.crystallography.net/cod/2300112.html) | 10.1107/S0021889805042457 | Sowa, Heidrun, Ahsbahs, Hans (2006) Journal of Applied Crystallography 39, 169-175 | 293 K 단결정 정밀화(R = 0.0149, Uiso 있음) — 고압 연구의 상압(블록 zno0) 측정, a·c 가 흔한 값(3.2498, 5.2066 Å)과 0.02 % 이내 |
| pack:Al2O3 | Corundum α-Al2O3 | Al2O3 (Z = 6) | R -3 c :H (#167) | a = 0.47602, c = 1.29933 | [9007496](https://www.crystallography.net/cod/9007496.html) | 10.1107/S0567739482001478 | Lewis, J., Schwarzenbach, D., Flack, H. D. (1982) Acta Crystallographica, Section A 38, 733-739 | 정밀 구조인자 연구(Lewis 1982), a = 4.7602, c = 12.9933 Å |
| pack:NaCl | Halite NaCl | NaCl (Z = 4) | F m -3 m (#225) | a = 0.564056 | [9008678](https://www.crystallography.net/cod/9008678.html) | — (CIF 에 없음) | Wyckoff, R. W. G. (1963) Crystal Structures 1, 85-237 | Wyckoff 1963 (AMCSD) 표준 참조 구조, a = 5.64056 Å |
| pack:MgO | Periclase MgO | MgO (Z = 4) | F m -3 m (#225) | a = 0.42112 | [9008671](https://www.crystallography.net/cod/9008671.html) | — (CIF 에 없음) | Wyckoff, R. W. G. (1963) Crystal Structures 1, 85-237 | Wyckoff 1963 (AMCSD) 표준 참조 구조, a = 4.2112 Å |
| pack:graphite-2H | Graphite-2H | C (Z = 4) | P 63/m m c (#194) | a = 0.2464, c = 0.6711 | [9011577](https://www.crystallography.net/cod/9011577.html) | 10.1038/258136a0 | Trucano, P., Chen, R. (1975) Nature 258, 136-137 | 천연 결정 중성자 회절(Trucano & Chen 1975), 2H 적층 P63/mmc |

## COD 에 없는 상

| 상 | 결과 |
|---|---|
| α-NaLuF4 | 없음 — Na·Lu·F 원소 검색(3–6 원소) 결과 NaLuF4 구조 없음 (2026-10-17 COD 검색) |
| β-NaLuF4 | 없음 — 같음 (2026-10-17 COD 검색) |
| NaGdF4 (육방) | 없음 — Na·Gd·F 원소 검색(3–6 원소) 결과 NaGdF4 구조 없음 (2026-10-17 COD 검색) |

## 확인할 점

- LiYF4 (2001129): CIF 의 Uiso(Y 0.100, Li 0.20, F 0.170 Å²)는 상온의 흔한 값(0.005–0.015 Å²)보다 열 배가량 커서 입력 오류로 보고 쓰지 않았다 (B = 0). 좌표·격자는 그대로.
- Fe3O4·루틸·아나타제·커런덤·흑연: CIF 에 비등방 U 만 있어 B = 0. 루틸·아나타제(9015662, 9015929)의 비등방 U 열은 이름(U11 U22 U33 U12 U13 U23)과 값 배열(U11 U12 U13 U22 U23 U33 로 읽어야 정방정 자리 대칭에 맞음)이 어긋나 있다 — 쓰지 않았으므로 팩에는 영향 없음.
- α-NaYF4 (1517676, Roy & Roy 1964): 오래된 격자상수(a = 5.47 Å, 유효숫자 3)이고 변위 매개변수 없음. 실험 시료의 α 상 격자상수는 조성·크기에 따라 다를 수 있다.
- β-NaYF4 (1517672): CIF 의 0.6667/0.3333 좌표를 2/3·1/3 로 바꾸지 않고 그대로 썼다. Na2 (2h) z = 0.689(5) 점유 0.5.
- Y2O3 (1009013): Zr 도핑 2 상 시료의 Y2O3 상 (화학식은 순수 Y2O3 로 기록). 측정 온도 미기재(상온으로 봄).
- Wyckoff (1963) 원소·암염 항목은 AMCSD 경유이며 DOI·측정 온도가 없는 것이 많다 (α-Fe 는 298 K 명시).

## 사용 조건

- COD 는 수록 자료를 공개 영역(public domain)으로 둔다. CIF 머리말 문구는 두 가지다: (1) "All data on this site have been placed in the public domain by the contributors." — 1000043, 1009013, 1517672, 1517676, 9008460; (2) AMCSD·IUCr 제공 항목(90xxxxx, 2001129, 2300112) — "The file may be used within the scientific community so long as proper attribution is given to the journal article from which the data were obtained." 그래서 각 상의 `source.citation`·`doi` 로 원 논문을 밝힌다.
