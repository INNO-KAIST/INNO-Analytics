"""StackReg/TurboReg 대조 자료 생성 — pystackreg 0.2.8 를 기준으로 삼는다.

    python3 -I test/stackreg_parity.py [출력폴더]      (기본: <tempdir>/stackreg_parity)
    node test/test_stackreg.js [같은 폴더]

만드는 것
  pair  : ref, mov (float32 또는 uint8 값) → StackReg(mode).register(ref, mov) 의 행렬·점,
          StackReg.transform(mov, 그 행렬) 결과 영상, 그리고 (알면) 실제로 적용한 변환(정답 행렬).
  stack : 6 장 이동 스택 → register_stack(reference, n_frames) 행렬, transform_stack 결과.
정답 행렬 T 의 뜻은 pystackreg 와 같다: mov(T·p) = ref(p), p=(x=열, y=행, 1).
합성 영상
  blob  : 무작위 가우스 덩어리의 합(해석식) + 약한 잡음. mov 는 해석식을 T^-1 q 에서 직접 계산 (보간 없음).
  canvas: 넓은 캔버스(잡음을 가우스로 흐림)를 scipy.ndimage.map_coordinates(3차 spline) 로 옮긴 뒤 잘라 냄.
"""
import json, os, sys, tempfile, time
import numpy as np
from scipy import ndimage
from pystackreg import StackReg

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(tempfile.gettempdir(), "stackreg_parity")
os.makedirs(OUT, exist_ok=True)
MODES = {"translation": StackReg.TRANSLATION, "rigid": StackReg.RIGID_BODY, "scaled": StackReg.SCALED_ROTATION,
         "affine": StackReg.AFFINE, "bilinear": StackReg.BILINEAR}


def wr(name, arr):
    np.ascontiguousarray(arr, dtype="<f8").tofile(os.path.join(OUT, name))
    return name


def blob_field(rng, W, H, neg=False):
    n = int(min(400, max(16, W * H / 150)))
    cx = rng.uniform(-0.1 * W, 1.1 * W, n); cy = rng.uniform(-0.1 * H, 1.1 * H, n)
    smin = 2.0 + 0.004 * min(W, H); smax = 3.0 + 0.06 * min(W, H)
    s = rng.uniform(smin, smax, n)
    A = rng.uniform(0.3, 1.0, n) * (rng.choice([-1.0, 1.0], n) if neg else 1.0)
    off = -0.6 if neg else 0.2

    def f(x, y):
        out = np.full(np.shape(x), off, dtype=np.float64)
        for i in range(n):
            out += A[i] * np.exp(-((x - cx[i]) ** 2 + (y - cy[i]) ** 2) / (2.0 * s[i] ** 2))
        return out
    return f


def T_about_center(W, H, lin, t):
    c = np.array([(W - 1) / 2.0, (H - 1) / 2.0]); lin = np.asarray(lin, float)
    T = np.eye(3); T[:2, :2] = lin; T[:2, 2] = c + np.asarray(t, float) - lin @ c
    return T


def rot(deg, s=1.0):
    a = np.deg2rad(deg); return s * np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])


def true_T(mode, W, H):
    if mode == "translation": return T_about_center(W, H, np.eye(2), (0.0731 * W, -0.0523 * H))
    if mode == "rigid": return T_about_center(W, H, rot(4.3), (2.3, -1.7))
    if mode == "scaled": return T_about_center(W, H, rot(-3.1, 1.035), (1.2, 2.1))
    if mode == "affine": return T_about_center(W, H, rot(2.2) @ np.array([[1.03, 0.035], [-0.02, 0.975]]), (-1.6, 1.3))
    return None


def grid(W, H):
    y, x = np.mgrid[0:H, 0:W].astype(np.float64); return x, y


def finish(img, rng, kind):
    img = img + rng.normal(0.0, 0.003 * (img.max() - img.min()), img.shape)
    if kind == "u8":
        lo, hi = img.min(), img.max()
        return np.clip(np.round(10 + 235 * (img - lo) / (hi - lo)), 0, 255).astype(np.uint8)
    return img.astype(np.float32)


def make_pair_blob(rng, mode, W, H, kind, neg):
    f = blob_field(rng, W, H, neg); x, y = grid(W, H)
    ref = f(x, y)
    if mode == "bilinear":  # 완만한 쌍선형 왜곡: p = q + d(q)
        cx, cy = (W - 1) / 2, (H - 1) / 2
        px = x + 1.4 + 0.012 * (x - cx) + 0.00055 * (x - cx) * (y - cy)
        py = y - 0.9 - 0.008 * (y - cy) + 0.0004 * (x - cx) * (y - cy)
        T = None
    else:
        T = true_T(mode, W, H); Ti = np.linalg.inv(T)
        px = Ti[0, 0] * x + Ti[0, 1] * y + Ti[0, 2]; py = Ti[1, 0] * x + Ti[1, 1] * y + Ti[1, 2]
    mov = f(px, py)
    # 같은 범위로 양자화되도록 두 장을 함께 다룸 (u8)
    both = finish(np.stack([ref, mov]), rng, kind)
    return both[0], both[1], T


def make_pair_canvas(rng, mode, W, H, P=24):
    C = ndimage.gaussian_filter(rng.normal(size=(H + 2 * P, W + 2 * P)), 2.5) * 10.0
    ref = C[P:P + H, P:P + W]
    T = true_T(mode, W, H); Ti = np.linalg.inv(T); x, y = grid(W, H)
    px = Ti[0, 0] * x + Ti[0, 1] * y + Ti[0, 2]; py = Ti[1, 0] * x + Ti[1, 1] * y + Ti[1, 2]
    mov = ndimage.map_coordinates(C, [py + P, px + P], order=3, mode="mirror", prefilter=True)
    return ref.astype(np.float32), mov.astype(np.float32), T


def gt_err(M, T, W, H):
    """정답과의 차: 영상 네 귀퉁이·중심에서 (M p − T p) 의 최대/중심 길이, 회전각 차(°)"""
    pts = np.array([[0, 0, 1], [W - 1, 0, 1], [0, H - 1, 1], [W - 1, H - 1, 1], [(W - 1) / 2, (H - 1) / 2, 1]], float).T
    d = np.hypot(*((M @ pts)[:2] - (T @ pts)[:2]))
    ang = np.degrees(np.arctan2(M[1, 0], M[0, 0]) - np.arctan2(T[1, 0], T[0, 0]))
    return {"dmax": float(d.max()), "dcenter": float(d[-1]), "angle_deg": float(abs(ang)),
            "scale": float(abs(np.sqrt(abs(np.linalg.det(M[:2, :2]))) - np.sqrt(abs(np.linalg.det(T[:2, :2])))))}


cases = []
rng = np.random.default_rng(20261016)
tpy = {}


def add_pair(cid, mode, ref, mov, T, kind):
    H, W = ref.shape
    sr = StackReg(MODES[mode])
    t0 = time.perf_counter(); mat = sr.register(ref, mov); dt = time.perf_counter() - t0
    refpts, movpts = sr.get_points()
    out = sr.transform(mov, mat)
    c = {"id": cid, "kind": "pair", "mode": mode, "w": W, "h": H, "bit": 8 if kind == "u8" else 32,
         "ref": wr(cid + "_ref.f64", ref.astype(np.float64)), "mov": wr(cid + "_mov.f64", mov.astype(np.float64)),
         "mat": mat.tolist(), "refPts": np.asarray(refpts).tolist(), "movPts": np.asarray(movpts).tolist(),
         "out": wr(cid + "_out.f64", out), "py_ms": dt * 1000.0}
    if T is not None:
        c["T"] = T.tolist(); c["gt_py"] = gt_err(mat, T, W, H)
    cases.append(c)
    return c


# 1) 다섯 변환 × 네 크기 (blob, float32; 음수 포함 판 섞음)
SIZES = [(64, 48), (100, 100), (37, 53), (256, 256)]
for mi, mode in enumerate(MODES):
    for si, (W, H) in enumerate(SIZES):
        neg = (mi + si) % 2 == 1
        ref, mov, T = make_pair_blob(rng, mode, W, H, "f32", neg)
        add_pair(f"blob_{mode}_{W}x{H}{'_neg' if neg else ''}", mode, ref, mov, T, "f32")
# 2) 8-bit 정수 영상
for mode in ["translation", "rigid", "affine", "bilinear"]:
    ref, mov, T = make_pair_blob(rng, mode, 100, 80, "u8", False)
    add_pair(f"u8_{mode}_100x80", mode, ref, mov, T, "u8")
# 3) scipy 캔버스 (map_coordinates 3차 spline 후 자름)
for mode in ["translation", "rigid", "scaled", "affine"]:
    ref, mov, T = make_pair_canvas(rng, mode, 128, 96)
    add_pair(f"canvas_{mode}_128x96", mode, ref, mov, T, "f32")
# 4) 512×512 (속도 측정용, 대조도 함)
for mode in ["rigid", "affine"]:
    ref, mov, T = make_pair_blob(rng, mode, 512, 512, "f32", True)
    c = add_pair(f"big_{mode}_512x512", mode, ref, mov, T, "f32")
    tpy[mode] = c["py_ms"]


# 5) 스택: 6 장, 점점 움직이는 영상
def make_stack(rng, W, H, n, kind):
    f = blob_field(rng, W, H, kind == "f32"); x, y = grid(W, H); frames = []
    Ts = []
    for i in range(n):
        T = T_about_center(W, H, rot(0.9 * i), (0.83 * i, -0.57 * i + 0.1 * i * i)); Ts.append(T); Ti = np.linalg.inv(T)
        frames.append(f(Ti[0, 0] * x + Ti[0, 1] * y + Ti[0, 2], Ti[1, 0] * x + Ti[1, 1] * y + Ti[1, 2]))
    return finish(np.stack(frames), rng, kind), Ts


stk32, Ts32 = make_stack(rng, 80, 64, 6, "f32")
stk8, Ts8 = make_stack(rng, 72, 60, 6, "u8")
STACKS = [("f32", stk32, Ts32), ("u8", stk8, Ts8)]
SCASES = [("rigid", "previous", 1), ("rigid", "first", 1), ("rigid", "first", 3), ("rigid", "mean", 1),
          ("translation", "previous", 1), ("scaled", "first", 3), ("affine", "mean", 1), ("bilinear", "first", 1), ("bilinear", "mean", 1)]
for kind, stk, Ts in STACKS:
    n, H, W = stk.shape
    fname = wr(f"stack_{kind}.f64", stk.astype(np.float64))
    for mode, refr, nf in SCASES:
        if kind == "u8" and mode in ("translation", "scaled", "bilinear") and refr != "mean": continue
        sr = StackReg(MODES[mode])
        mats = sr.register_stack(stk, reference=refr, n_frames=nf, axis=0, suppress_axis_warning=True)
        if refr == "first" and nf > 1:   # 앱 규칙: 첫 n 장 평균 기준이면 첫 장도 그 기준에 맞춘다 (pystackreg 는 첫 장을 맞추지 않음)
            mats[0] = sr.register(np.mean(stk.take(range(nf), axis=0), axis=0), stk[0])
        out = sr.transform_stack(stk, tmats=mats)
        cid = f"stack_{kind}_{mode}_{refr}{nf if refr == 'first' else ''}"
        cases.append({"id": cid, "kind": "stack", "mode": mode, "w": W, "h": H, "n": n, "bit": 8 if kind == "u8" else 32,
                      "frames": fname, "reference": refr, "nFrames": nf, "mats": mats.tolist(), "out": wr(cid + "_out.f64", out)})

with open(os.path.join(OUT, "manifest.json"), "w") as fh:
    json.dump({"pystackreg": "0.2.8", "numpy": np.__version__, "cases": cases}, fh)

print(f"{len(cases)} 개 사례 → {OUT}")
print("pystackreg 자체의 정답 복원 오차 (px, °):")
for c in cases:
    if "gt_py" in c:
        g = c["gt_py"]; print(f"  {c['id']:34s} dmax {g['dmax']:.4f}  dcenter {g['dcenter']:.4f}  angle {g['angle_deg']:.4f}  scale {g['scale']:.2e}")
print("pystackreg 512×512 register 시간: " + ", ".join(f"{k} {v:.0f} ms" for k, v in tpy.items()))
