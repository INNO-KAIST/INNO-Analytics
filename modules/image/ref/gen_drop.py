#!/usr/bin/env python3
"""Drop Shape 시험 자료 — Young–Laplace 정답 윤곽(scipy, 앱과 다른 구현)으로 펜던트·세실 방울 그림을 만든다.
  python3 -I ref/gen_drop.py            → test/drop/*.png + truth.json (seed 1)

무차원 Young–Laplace (apex 곡률 반지름 b 로 나눔, s = 호 길이, φ = 접선 각, x = 반지름, z = apex 에서 높이):
  dx/ds = cos φ, dz/ds = sin φ,
  펜던트(apex 아래, z 위로):  dφ/ds = 2 − Bo·z − sin φ / x
  세실(apex 위, z 아래로):    dφ/ds = 2 + Bo·z − sin φ / x      (Bo = Δρ g b² / γ, apex 에서 sin φ / x → 1)
렌더: 4×4 부분 화소 면적 → 방울·바늘·기판 어둡게, 바탕 밝게 + 가우스 흐림 + 잡음 (8-bit).
정답: 펜던트 γ = Δρ g b² / Bo, 부피·겉넓이, Worthington Wo = Δρ g V / (π γ D_needle). 세실: 접촉각 θ, 높이·폭, 부피.
"""
import json, math, os, argparse
import numpy as np
from PIL import Image, ImageFilter
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "test", "drop")


def yl(Bo, mode, s_max=12.0, events=None):
    sgn = -1.0 if mode == "pendant" else 1.0
    def f(s, y):
        x, z, p = y
        sx = math.sin(p) / x if x > 1e-9 else 1.0
        return [math.cos(p), math.sin(p), 2.0 + sgn * Bo * z - sx]
    s0 = 1e-6  # 원점 특이점: 작은 s 에서 급수로 시작 (x≈s, z≈s²/2, φ≈s)
    y0 = [s0, s0 * s0 / 2, s0]
    return solve_ivp(f, (s0, s_max), y0, method="DOP853", rtol=1e-11, atol=1e-12, dense_output=True, events=events, max_step=0.01)


def profile_pendant(Bo, u):
    """apex 에서 바늘 끝까지. 적도(x 최대)를 지나 x 가 줄어드는 구간에서, 목(x 최소, 없으면 x → 0 직전) 과 0.8·x_max 사이
    u ∈ (0,1) 비율 자리를 바늘 끝으로 삼는다 → rn = 그 자리의 x (실제 펜던트 방울처럼 바늘 지름 = 윤곽이 붙는 폭)"""
    sol = yl(Bo, "pendant", 30.0)
    s = np.linspace(1e-6, sol.t[-1], 60000); x, z, p = sol.sol(s)
    dx = np.diff(x); i_max = int(np.argmax(dx < 0)) if (dx < 0).any() else len(x) - 1   # 첫 극대(적도)
    j = i_max
    while j + 1 < len(x) and x[j + 1] < x[j] and x[j + 1] > 0.05 and p[j] < math.pi * 0.98:
        j += 1
    x_lo, x_hi = x[j] * 1.03, 0.8 * x[i_max]
    if not x_hi > x_lo:
        raise ValueError("목이 너무 넓음")
    rn = x_hi - u * (x_hi - x_lo)
    k = i_max + int(np.argmax(x[i_max:j + 1] <= rn))
    return s[:k + 1], x[:k + 1], z[:k + 1], p[:k + 1], float(rn)


def profile_sessile(Bo, theta_deg):
    th = math.radians(theta_deg)
    def ev(s, y):
        return y[2] - th
    ev.terminal = True; ev.direction = 1
    sol = yl(Bo, "sessile", 30.0, [ev])
    s_end = sol.t_events[0][0]
    s = np.linspace(1e-6, s_end, 20000)
    x, z, p = sol.sol(s)
    return s, x, z, p


def vol_area(x, z):
    """회전체: V = ∫ π x² dz, A = ∫ 2π x ds"""
    dz = np.diff(z); xm = (x[1:] + x[:-1]) / 2; ds = np.hypot(np.diff(x), dz)
    return float(np.sum(math.pi * xm * xm * dz)), float(np.sum(2 * math.pi * xm * ds))


def render(polys, W, H, bg, fg, rng, blur, noise, ss=4):
    """polys: [(xs, ys)] 이미지 좌표 다각형(닫힘) — 안쪽 = fg. 부분 화소 ss×ss 로 면적"""
    big = Image.new("L", (W * ss, H * ss), 0)
    from PIL import ImageDraw
    d = ImageDraw.Draw(big)
    for xs, ys in polys:
        d.polygon([(float(a) * ss, float(b) * ss) for a, b in zip(xs, ys)], fill=255)
    m = np.asarray(big.resize((W, H), Image.BOX), np.float64) / 255.0
    img = bg * (1 - m) + fg * m
    img = np.asarray(Image.fromarray(np.clip(img, 0, 255).astype(np.uint8)).filter(ImageFilter.GaussianBlur(blur)), np.float64) if blur > 0 else img
    img = img + rng.normal(0, noise, img.shape)
    return np.clip(np.round(img), 0, 255).astype(np.uint8)


def make_pendant(rng, k, b_mm, Bo, u, px_per_mm, tilt_deg, out):
    s, x, z, p, rnb = profile_pendant(Bo, u); rn_mm = rnb * b_mm
    X = x * b_mm * px_per_mm; Z = z * b_mm * px_per_mm  # px
    H = int(Z[-1] + 0.35 * Z[-1] + 80); W = int(2 * X.max() + 160)
    cx = W / 2 + rng.uniform(-0.5, 0.5); ay = H - 40 + rng.uniform(-0.5, 0.5)   # apex (이미지 y 아래로)
    # 방울 윤곽 (왼쪽 아래→위, 오른쪽 위→아래) + 바늘(위로 영상 끝까지)
    xs = np.concatenate([cx - X[::-1], cx + X]); ys = np.concatenate([ay - Z[::-1], ay - Z])
    top = ay - Z[-1]
    rnp = rn_mm * px_per_mm
    xs = np.concatenate([[cx - rnp, cx - rnp], xs[1:-1] if False else xs, [cx + rnp, cx + rnp]])
    ys = np.concatenate([[-5, top], ys, [top, -5]])
    # 기울임 (카메라가 기울어짐): apex 둘레로 회전
    t = math.radians(tilt_deg); c, sn = math.cos(t), math.sin(t)
    xr = cx + (xs - cx) * c - (ys - ay) * sn; yr = ay + (xs - cx) * sn + (ys - ay) * c
    bg, fg = rng.uniform(150, 210), rng.uniform(15, 45)
    img = render([(xr, yr)], W, H, bg, fg, rng, rng.uniform(0.6, 1.2), rng.uniform(1.0, 3.0))
    fn = f"pendant_{k:02d}.png"; Image.fromarray(img).save(os.path.join(out, fn))
    V, A = vol_area(x * b_mm, z * b_mm)   # mm³, mm²
    drho, g = 1000.0, 9.80665
    gamma = drho * g * (b_mm * 1e-3) ** 2 / Bo * 1e3   # mN/m
    Wo = drho * g * (V * 1e-9) / (math.pi * gamma * 1e-3 * 2 * rn_mm * 1e-3)
    return dict(file=fn, kind="pendant", W=W, H=H, b_mm=b_mm, Bo=Bo, rn_mm=rn_mm, needle_mm=2 * rn_mm, px_per_mm=px_per_mm, tilt_deg=tilt_deg,
                apex=[cx, ay], drho=drho, g=g, gamma_mNm=gamma, V_mm3=V, A_mm2=A, Wo=Wo, Rmax_mm=float(x.max() * b_mm))


def make_sessile(rng, k, b_mm, Bo, theta, px_per_mm, tilt_deg, out):
    s, x, z, p = profile_sessile(Bo, theta)
    X = x * b_mm * px_per_mm; Z = z * b_mm * px_per_mm
    h = Z[-1]; W = int(2 * X.max() + 200); H = int(h + 260)
    cx = W / 2 + rng.uniform(-0.5, 0.5); base_y = H - 120 + rng.uniform(-0.5, 0.5); ay = base_y - h
    xs = np.concatenate([cx - X[::-1], cx + X]); ys = np.concatenate([ay + Z[::-1], ay + Z])
    # 기판: 기준선 아래 전체 (방울과 붙음)
    sub_x = np.array([-10, W + 10, W + 10, -10]); sub_y = np.array([base_y, base_y, H + 10, H + 10])
    t = math.radians(tilt_deg); c, sn = math.cos(t), math.sin(t)
    rot = lambda X_, Y_: (cx + (X_ - cx) * c - (Y_ - base_y) * sn, base_y + (X_ - cx) * sn + (Y_ - base_y) * c)
    d1 = rot(xs, ys); d2 = rot(sub_x, sub_y)
    bg, fg = rng.uniform(150, 210), rng.uniform(15, 45)
    img = render([d1, d2], W, H, bg, fg, rng, rng.uniform(0.6, 1.2), rng.uniform(1.0, 3.0))
    fn = f"sessile_{k:02d}.png"; Image.fromarray(img).save(os.path.join(out, fn))
    V, A = vol_area(x * b_mm, z * b_mm)
    cL = rot(np.array([cx - X[-1]]), np.array([base_y])); cR = rot(np.array([cx + X[-1]]), np.array([base_y]))
    return dict(file=fn, kind="sessile", W=W, H=H, b_mm=b_mm, Bo=Bo, theta=theta, px_per_mm=px_per_mm, tilt_deg=tilt_deg,
                baseline=[[float(cL[0][0]), float(cL[1][0])], [float(cR[0][0]), float(cR[1][0])]], height_mm=float(z[-1] * b_mm),
                width_mm=float(2 * x[-1] * b_mm), V_mm3=V)


def make_sessile_needle(rng, k, b_mm, Bo, theta, px_per_mm, out):
    """바늘이 방울 위로 꽂힌 세실 (needle-in-drop: 전진·후퇴각 측정 배치). 바늘 = 영상 위 끝에서 방울 높이 40 % 까지 내려온 사각형."""
    s, x, z, p = profile_sessile(Bo, theta)
    X = x * b_mm * px_per_mm; Z = z * b_mm * px_per_mm
    h = Z[-1]; W = int(2 * X.max() + 200); H = int(h + 300)
    cx = W / 2 + rng.uniform(-0.5, 0.5); base_y = H - 120 + rng.uniform(-0.5, 0.5); ay = base_y - h
    xs = np.concatenate([cx - X[::-1], cx + X]); ys = np.concatenate([ay + Z[::-1], ay + Z])
    sub_x = np.array([-10, W + 10, W + 10, -10]); sub_y = np.array([base_y, base_y, H + 10, H + 10])
    rn = 0.15 * X[-1] if theta <= 90 else 0.15 * X.max()
    nd_x = np.array([cx - rn, cx + rn, cx + rn, cx - rn]); nd_y = np.array([-10, -10, ay + 0.4 * h, ay + 0.4 * h])
    bg, fg = rng.uniform(150, 210), rng.uniform(15, 45)
    img = render([(xs, ys), (sub_x, sub_y), (nd_x, nd_y)], W, H, bg, fg, rng, rng.uniform(0.6, 1.2), rng.uniform(1.0, 3.0))
    fn = f"sessile_needle_{k:02d}.png"; Image.fromarray(img).save(os.path.join(out, fn))
    return dict(file=fn, kind="sessile_needle", W=W, H=H, b_mm=b_mm, Bo=Bo, theta=theta, px_per_mm=px_per_mm, tilt_deg=0.0,
                baseline=[[cx - float(X[-1]), base_y], [cx + float(X[-1]), base_y]], needle_px=2 * rn)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seed", type=int, default=1); a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(a.seed); truth = []
    # 펜던트: γ 25–75 mN/m 에서 b = √(Bo γ / Δρ g), 바늘 끝 = 적도와 목 사이 임의 자리, 120–260 px/mm, 기울임 ±1.5°
    k = 0
    for Bo in [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]:
        for rep in range(2):
            gam = rng.uniform(25, 75); b = math.sqrt(Bo * gam * 1e-3 / (1000.0 * 9.80665)) * 1e3   # mm
            u = rng.uniform(0.1, 0.9)
            try:
                truth.append(make_pendant(rng, k, b, Bo, u, rng.uniform(120, 260) * (1.2 / b), rng.uniform(-1.5, 1.5), OUT)); k += 1
            except ValueError:
                pass
    # 세실: θ 20–160°, Bo 0–0.8, 150–300 px/mm, 기울임 ±3°
    k = 0
    for th in [20, 40, 60, 75, 90, 105, 120, 135, 150, 160]:
        for rep in range(2):
            Bo = float(rng.choice([0.0005, 0.1, 0.3, 0.6]))
            truth.append(make_sessile(rng, k, rng.uniform(0.5, 1.4), Bo, th, rng.uniform(150, 300), rng.uniform(-3, 3), OUT)); k += 1
    # 바늘이 꽂힌 세실 (앞의 그림들의 난수 순서는 그대로 — 뒤에 덧붙임)
    for k, (th, Bo) in enumerate([(70, 0.3), (95, 0.1), (120, 0.3)]):
        truth.append(make_sessile_needle(rng, k, rng.uniform(0.8, 1.2), Bo, th, rng.uniform(180, 260), OUT))
    json.dump(truth, open(os.path.join(OUT, "truth.json"), "w"), indent=1)
    print(len(truth), "장 →", OUT)


if __name__ == "__main__":
    main()
