"""공간군 530 설정(Hall 번호)의 대칭 연산 표 — spglib 2.x 데이터베이스(BSD-3; 내용은 International Tables A)에서 만든다.
출력: src/diff_sg.js  (window/worker 의 DIFF_SG).  사용: python3 -I ref/gen_sg.py"""
import spglib, json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
rots = {}; rot_list = []; entries = []
for hall in range(1, 531):
    t = spglib.get_spacegroup_type(hall); s = spglib.get_symmetry_from_database(hall)
    ops = []
    for R, tr in zip(s['rotations'], s['translations']):
        key = ','.join(str(int(v)) for v in R.flatten())
        if key not in rots: rots[key] = len(rot_list); rot_list.append([int(v) for v in R.flatten()])
        t24 = [int(round(v * 24)) % 24 for v in tr]
        assert all(abs(v * 24 - round(v * 24)) < 1e-6 for v in tr)
        ops.append(rots[key] * 13824 + t24[0] * 576 + t24[1] * 24 + t24[2])   # 한 수로: 회전 번호·평행이동(24 분의 1)
    entries.append([hall, t.number, t.international_full, t.international_short, t.choice, t.hall_symbol, ops])
print('rotations', len(rot_list), 'ops', sum(len(e[6]) for e in entries))
js = ("/* 공간군 530 설정 — spglib " + spglib.__version__ + " 데이터베이스(BSD-3, International Tables A) 에서 ref/gen_sg.py 로 만든 표. 손으로 고치지 않는다.\n"
      "   R: 회전 행렬(행 우선 9 수) 목록. S: [hall, 번호, 완전 기호, 짧은 기호, 설정(choice), Hall 기호, 연산들] — 연산 = 회전번호·13824 + tx·576 + ty·24 + tz (평행이동은 1/24 단위) */\n"
      "(function(){ const G=(typeof self!=='undefined')?self:globalThis;\nG.DIFF_SG={R:" + json.dumps(rot_list, separators=(',', ':')) + ",\nS:" + json.dumps(entries, separators=(',', ':'), ensure_ascii=False) + "};\n"
      "if(typeof module!=='undefined'&&module.exports) module.exports=G.DIFF_SG; })();\n")
open(os.path.join(HERE, '..', 'src', 'diff_sg.js'), 'w').write(js)
print('size', len(js))
