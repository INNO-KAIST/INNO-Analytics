"""전자 산란 인자 표 — International Tables for Crystallography Vol. C, Table 4.3.2.2 (Peng, Ren, Dudarev & Whelan 1996, Acta Cryst. A52, 257):
f_e(s) = Σ_{i=1..5} a_i exp(-b_i s²) [Å], s = sinθ/λ ≤ 2 Å⁻¹ (중성 원자, 탄성). 수치는 gemmi 의 C4322 표에서 읽는다(표 값 자체는 사실 자료).
원자량은 gemmi Element.weight. 출력: src/diff_scat.js (DIFF_SCAT)"""
import gemmi, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
T = {}
for z in range(1, 119):
    e = gemmi.Element(z); c = e.c4322
    if c is None: continue
    v = [round(x, 4) for x in c.get_coefs()]
    T[e.name] = [z, round(e.weight, 4)] + v
print(len(T), 'elements; last', list(T)[-1])
js = ("/* 전자 산란 인자 (탄성, 중성 원자) — International Tables Vol. C Table 4.3.2.2 (Peng et al. 1996, Acta Cryst A52 257), gemmi 표에서 ref/gen_scat.py 로.\n"
      "   원소: [Z, 원자량, a1..a5 (Å), b1..b5 (Å²)]   f_e(s)=Σ a_i exp(-b_i s²), s=sinθ/λ=1/(2d) */\n"
      "(function(){ const G=(typeof self!=='undefined')?self:globalThis;\nG.DIFF_SCAT=" + json.dumps(T, separators=(',', ':')) + ";\n"
      "if(typeof module!=='undefined'&&module.exports) module.exports=G.DIFF_SCAT; })();\n")
open(os.path.join(HERE, '..', 'src', 'diff_scat.js'), 'w').write(js); print('size', len(js))
