"""실측 Leica LIF 의 정답 — Bio-Formats 7.3.1 (ref/bf/BFDump.java) 이 읽은 series 별 화소·크기·보정.
liffile(C. Gohlke, 독립 구현) 이 읽은 배열과도 같은지 확인한다(두 구현이 같은 값일 때만 정답으로).
사용: python3 -I gen_lif.py <LIF 폴더> <bioformats_package.jar> <BFDump 컴파일 폴더>
파일: sh fetch_lif.sh <폴더>  (readlif·bioio-lif 시험 자료, GitHub)"""
import sys, os, json, glob, subprocess, tempfile, itertools
import numpy as np, liffile
D, JAR, OUT = sys.argv[1:4]
T = os.path.join(D, 'truth'); os.makedirs(T, exist_ok=True); truth = {}
def bf(path, series=None):
    td = tempfile.mkdtemp(); cmd = ['java', '-Xmx6g', '-cp', JAR + ':' + OUT, 'BFDump', path] + ([td, ','.join(map(str, series))] if series else [])
    r = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, JAVA_TOOL_OPTIONS=''))
    S = {}
    for l in r.stdout.splitlines():
        if l.startswith('S '):
            t = l.split(); d = dict(zip(t[2::2], t[3::2])); S[int(t[1])] = d
    return S, td
for p in sorted(glob.glob(os.path.join(D, '*.lif'))):
    fn = os.path.basename(p)
    S0, _ = bf(p)                      # 메타데이터만
    n = len(S0); pick = list(range(n)) if n <= 6 else [0, 1, 2, n - 1]
    S, td = bf(p, pick)
    # liffile: (요소, tile) 순서로 BF series 와 짝
    lif = liffile.LifFile(p); pairs = []
    for im in lif.images:
        m = im.sizes.get('M', 1)
        for k in range(m): pairs.append((im, k))
    assert len(pairs) == n, (fn, len(pairs), n)
    for s in pick:
        d = S[s]; X, Y, Z, C, Tt = (int(d[k]) for k in 'X Y Z C T'.split())
        a = np.fromfile(os.path.join(td, f's{s}.f32'), '<f4').reshape(Tt, Z, C, Y, X)
        im, k = pairs[s]; arr = np.asarray(im.asarray(), np.float32); dims = list(im.dims)
        if 'M' in dims: arr = np.take(arr, k, axis=dims.index('M')); dims.remove('M')
        # 영상 면 = 마지막 두 차원, 그 앞: Z(남은 공간 차원), C, 그리고 나머지는 T 로 — BF 의 순서와 맞는 배치를 찾는다
        plane = dims[-2:]; other = [x for x in dims[:-2]]
        sp = [x for x in other if x in 'ZY' and x not in plane]; zax = sp[0] if sp else None
        tax = [x for x in other if x not in ('C', zax)]
        ok = False
        for perm in itertools.permutations(tax):
            order = list(perm) + ([zax] if zax else []) + (['C'] if 'C' in dims else []) + plane
            b = np.transpose(arr, [dims.index(x) for x in order]).reshape(Tt, Z, C, Y, X)
            if np.array_equal(a, b): ok = True; break
        assert ok, (fn, s, dims, 'liffile ≠ Bio-Formats')
        key = f'{fn}#{s}'; a.tofile(os.path.join(T, key + '.truth'))
        truth[key] = dict(w=X, h=Y, nz=Z, nc=C, nt=Tt, px=float(d['px']) if d['px'] != '-' else None, py=float(d['py']) if d['py'] != '-' else None,
                          pz=float(d['pz']) if d['pz'] != '-' else None, name=d['name'], series=n, tfold=[x for x in perm] if tax else [])
    print(fn, n, 'series → 정답', len(pick), '(Bio-Formats = liffile)')
json.dump(truth, open(os.path.join(T, 'truth.json'), 'w'), ensure_ascii=False, indent=0)
print('정답', len(truth))
