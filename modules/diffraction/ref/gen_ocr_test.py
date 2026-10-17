#!/usr/bin/env python3
"""스케일바 이름표 OCR 시험 그림 만들기 (Pillow) → test/ocr/*.png + truth.json

틀(gen_ocr_templates.py)과 겹치지 않게: 크기는 틀 크기(12 15 19 24 30 40 52 px)를 피하고,
'heldout' 조건은 틀에 없는 글꼴(Carlito · Inter · Noto Sans CJK 라틴 · TeX Gyre Heros)만 쓴다.
나머지 조건은 실제 GMS 처럼 Arial 치수의 Liberation Sans (틀 글꼴과 같은 가족, 크기만 다름).
위첨자 'nm⁻¹' 은 틀(0.62·0.36)과 다른 비율(0.58·0.40)로 그린다.

  python3 ref/gen_ocr_test.py                  # 시험 묶음 (seed 1) → test/ocr
  python3 ref/gen_ocr_test.py --dev OUT --seed 7   # 맞춤용 dev 묶음 (다른 seed·다른 미사용 글꼴)
"""
import argparse, io, json, os
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"
HELD_TEST = [("Carlito", "/usr/share/fonts/truetype/crosextra/Carlito-Regular.ttf", 0),
             ("Inter", "/usr/share/fonts/opentype/inter/Inter-Regular.otf", 0),
             ("NotoSansCJK", "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", 0),
             ("TeXGyreHeros", "/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreheros-regular.otf", 0)]
HELD_DEV = [("IPAPGothic", "/usr/share/fonts/opentype/ipafont-gothic/ipagp.ttf", 0),
            ("DejaVuSansCondensed", "/usr/share/fonts/truetype/dejavu/DejaVuSansCondensed.ttf", 0),
            ("LiberationSansB", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 0)]
VALUES = ["0.5", "1", "2", "5", "10", "20", "50", "100", "200", "500"]
ODD = ["0.25", "0.2", "1.5", "3", "4", "25", "40", "250", "300", "7.5"]   # 1·2·5 가 아닌 값 (값 사전이 맞는 읽기를 덮지 않는지)
UNITS = [("1/nm", "1/nm"), ("nm", "nm"), ("nm-1", "1/nm"), ("1/Å", "1/Å"), ("µm", "µm"), ("Å", "Å")]
SIZES = [14, 17, 21, 27, 34, 44, 58, 70, 83]   # 대문자 높이 ≈ 0.72 × 크기 → 10–60 px
CONDS = ["odd", "clean", "below", "side", "box_dark", "box_white", "inverted", "degraded", "heldout", "small", "large", "hollow", "distract", "gms"]
# gms: GMS 기본 이름표 꼴 — Times 꼴 명조 검은 글자 + 흰 그림자(오른쪽 아래로 밀림), 흰 속 빈 막대(속은 검정), 어두운 잡음 바탕 (REDI 예제 B3·T13 실측 모양).
#      글꼴은 Liberation Serif(Times 치수 = 틀과 같은 가족, 크기만 다름) 절반 + 틀에 없는 명조(아래) 절반
GMS_TEST = [("LiberationSerif", "/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf", 0), ("Caladea", "/usr/share/fonts/truetype/crosextra/Caladea-Regular.ttf", 0),
            ("DejaVuSerifCondensed", "/usr/share/fonts/truetype/dejavu/DejaVuSerifCondensed.ttf", 0)]
GMS_DEV = [("LiberationSerif", "/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf", 0), ("GFSBaskerville", "/usr/share/fonts/truetype/baskerville/GFSBaskerville.otf", 0),
           ("FreeSerifBold", "/usr/share/fonts/truetype/freefont/FreeSerifBold.ttf", 0)]
DISTRACT = ["SAD 1", "200 kV", "DF-2", "x 40k", "Spot 3", "img_0042.dm3", "CL 250 mm", "HAADF"]


def font(path, size, idx=0):
    return ImageFont.truetype(path, size, index=idx)


def text_layer(W, H, x, base, parts, aa=True):
    """parts = [(문자열, 글꼴, 올림 px)]: 기준선 base 에서 왼쪽 x 부터 이어 그린다 → (잉크 배열, bbox)"""
    im = Image.new("L", (W, H), 0)
    d = ImageDraw.Draw(im)
    d.fontmode = "L" if aa else "1"
    for s, f, rise in parts:
        d.text((x, base - rise), s, font=f, fill=255, anchor="ls")
        x += d.textlength(s, font=f)
    return np.asarray(im, np.float64) / 255.0


def label_parts(val, unit, fp, idx, size, sup):
    f = font(fp, size, idx)
    if unit == "nm-1" and sup:
        fs = font(fp, max(6, round(size * 0.58)), idx)
        return [(val + " nm", f, 0), ("-1", fs, round(size * 0.40))]
    return [(val + " " + unit, f, 0)]


def measure(parts):
    im = Image.new("L", (10, 10)); d = ImageDraw.Draw(im)
    w = sum(d.textlength(s, font=f) for s, f, _ in parts)
    f0 = parts[0][1]
    asc, desc = f0.getmetrics()
    return int(np.ceil(w)), asc, desc


def background(rng, W, H, level, noise, nspots):
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float64)
    bg = level + rng.uniform(-8, 8) * (xx / W) + rng.uniform(-8, 8) * (yy / H)
    # 중앙빔 번짐 비슷한 큰 흐림
    if rng.random() < 0.5:
        cx, cy, s = rng.uniform(0, W), rng.uniform(0, H), rng.uniform(60, 250)
        bg += rng.uniform(10, 60) * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * s * s))
    for _ in range(nspots):
        cx, cy, s = rng.uniform(0, W), rng.uniform(0, H), rng.uniform(1.0, 4.0)
        bg += rng.uniform(30, 260) * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * s * s))
    bg += rng.normal(0, noise, bg.shape)
    return bg


def make(rng, cond, val, unit, fname, fp, fidx, size, out_png):
    aa = rng.random() > 0.15   # 15 % 는 계단 글자(앤티에일리어싱 없음)
    sup = unit == "nm-1" and rng.random() < 0.5
    parts = label_parts(val, unit, fp, fidx, size, sup)
    tw, asc, desc = measure(parts)
    capH = font(fp, size, fidx).getbbox("0", anchor="ls")
    capH = -capH[1]
    L = int(rng.uniform(max(60, 0.5 * tw), max(120, 2.5 * tw)))
    if cond == "large":
        L = int(rng.uniform(max(200, tw * 0.8), max(260, tw * 1.6)))
    thick = max(2, int(round(capH * rng.uniform(0.18, 0.45))))
    hollow = cond in ("hollow", "gms")
    if hollow:
        thick = max(6, int(round(capH * rng.uniform(0.35, 0.6))))
    place = "above"
    if cond == "below": place = "below"
    elif cond == "side": place = rng.choice(["left", "right"])
    elif cond in ("box_dark", "box_white"): place = rng.choice(["above", "above", "below", "right"])
    gap = int(round(capH * rng.uniform(0.15, 0.7)))
    m = int(rng.uniform(16, 50))
    # 배치 (캔버스 좌표는 나중에 여백 더함)
    if place in ("above", "below"):
        bw_ = max(L, tw); bx = (bw_ - L) // 2; tx = (bw_ - tw) // 2 + int(rng.uniform(-0.2, 0.2) * max(0, L - tw))
        if place == "above":
            base = asc + 2; by = base + desc + gap
        else:
            by = 0; base = thick + gap + capH + 2
        cw, ch = bw_, max(by + thick, base + desc) + 2
    else:
        sgap = int(round(capH * rng.uniform(0.3, 1.0)))
        if place == "right":
            bx, tx = 0, L + sgap
        else:
            tx, bx = 0, tw + sgap
        base = asc + 2; by = base - capH // 2 - thick // 2
        cw, ch = L + sgap + tw, base + desc + 2
        if by < 0:
            base -= by; by = 0; ch -= 0
    extra_top = 0
    if cond == "distract":
        extra_top = int(capH * rng.uniform(3.5, 6.0)) + asc
    W = cw + 2 * m; H = ch + 2 * m + extra_top
    W = max(W, 160); H = max(H, 100)
    ox, oy = m, m + extra_top
    level = rng.uniform(10, 45)
    noise = rng.uniform(0.7, 2.2)
    if cond == "gms":   # 실측 B3·T13: 바탕 중앙값 11–13, MAD 4 안팎
        level, noise = rng.uniform(8, 20), rng.uniform(3, 6)
    img = background(rng, W, H, level, noise, int(rng.uniform(3, 12) * W * H / 1e5) + 2)
    box = None
    if cond in ("box_dark", "box_white"):
        pad = int(rng.uniform(4, 14))
        box = (ox - pad, oy - pad, ox + cw + pad, oy + ch + pad)
        img[box[1]:box[3], box[0]:box[2]] = 255.0 if cond == "box_white" else rng.uniform(0, 8)
    ink = text_layer(W, H, ox + tx, oy + base, parts, aa)
    barm = np.zeros((H, W))
    X0, Y0 = ox + bx, oy + by
    barm[Y0:Y0 + thick, X0:X0 + L] = 1
    if hollow:
        fr = int(rng.choice([2, 3]))
        barm[Y0 + fr:Y0 + thick - fr, X0 + fr:X0 + L - fr] = 0
    if cond == "gms":
        # 흰 그림자(글자를 (sx, sy) 밀어 흰색) → 그 위에 검은 글자. 막대: 흰 테두리 + 검은 속
        sx = sy = max(1, int(round(capH * rng.uniform(0.03, 0.06))))
        sh = text_layer(W, H, ox + tx + sx, oy + base + sy, parts, aa)
        bo = np.zeros((H, W)); bo[Y0:Y0 + thick, X0:X0 + L] = 1
        fr = max(2, int(round(thick * 0.12)))
        inner = np.zeros((H, W)); inner[Y0 + fr:Y0 + thick - fr, X0 + fr:X0 + L - fr] = 1
        img = img * (1 - sh) + 255.0 * sh
        img = img * (1 - ink) + 0.0 * ink
        img = img * (1 - bo) + 255.0 * bo
        img = img * (1 - inner)
    else:
        lay = np.maximum(ink, barm)
        inkv = 0.0 if cond == "box_white" else 255.0
        img = img * (1 - lay) + inkv * lay
    if cond == "distract":
        dtext = DISTRACT[int(rng.integers(len(DISTRACT)))]
        dl = text_layer(W, H, int(rng.uniform(4, max(5, W - 200))), int(rng.uniform(asc + 2, max(asc + 3, extra_top - capH))), [(dtext, font(fp, size, fidx), 0)], aa)
        img = img * (1 - dl) + 255 * dl
    if cond == "degraded":
        img = np.clip(img, 0, 255)
        pim = Image.fromarray(img.astype(np.uint8)).filter(ImageFilter.GaussianBlur(rng.uniform(0.5, 1.0)))
        img = np.asarray(pim, np.float64) + rng.normal(0, rng.uniform(3, 7), (H, W))
        img = np.clip(img, 0, 255)
        buf = io.BytesIO(); Image.fromarray(img.astype(np.uint8)).save(buf, "JPEG", quality=int(rng.uniform(55, 80)))
        img = np.asarray(Image.open(io.BytesIO(buf.getvalue())).convert("L"), np.float64)
    if cond == "inverted":
        img = 255 - np.clip(img, 0, 255)
    img = np.clip(np.round(img), 0, 255).astype(np.uint8)
    Image.fromarray(img).save(out_png, optimize=True)
    return dict(W=W, H=H, bar=dict(x0=X0, x1=X0 + L - 1, y0=Y0, y1=Y0 + thick - 1), place=place, sup=sup, aa=aa, capH=capH, box=box)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dev", default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--conds", default=None)
    a = ap.parse_args()
    out = a.dev or os.path.join(HERE, "..", "test", "ocr")
    os.makedirs(out, exist_ok=True)
    held = HELD_DEV if a.dev else HELD_TEST
    rng = np.random.default_rng(a.seed)
    truth = []
    k = 0
    for cond in (a.conds.split(",") if a.conds else CONDS):
        for rep in range(a.reps):
            for val in (ODD if cond == "odd" else VALUES):
                for unit, canon in UNITS:
                    fname, fp, fidx = "LiberationSans", LIB, 0
                    if cond == "heldout":
                        fname, fp, fidx = held[int(rng.integers(len(held)))]
                    if cond == "gms":
                        gl = GMS_DEV if a.dev else GMS_TEST
                        fname, fp, fidx = gl[int(rng.integers(len(gl)))]
                    size = int(rng.choice(SIZES))
                    if cond == "small": size = int(rng.choice([13, 14, 16, 17]))
                    if cond == "large": size = int(rng.choice([62, 70, 76, 83]))
                    if cond == "clean": size = int(rng.choice(SIZES[:7]))
                    fn = "%s_%03d.png" % (cond, k); k += 1
                    meta = make(rng, cond, val, unit, fname, fp, fidx, size, os.path.join(out, fn))
                    v = float(val)
                    inv = {"1/nm": v, "1/Å": v * 10}.get(canon)
                    truth.append(dict(file=fn, cond=cond, text=val + " " + unit, value=v, unit=canon, valueInvNm=inv, font=fname, size=size, **meta))
    json.dump(truth, open(os.path.join(out, "truth.json"), "w"), ensure_ascii=False, indent=0)
    tot = sum(os.path.getsize(os.path.join(out, t["file"])) for t in truth)
    print(len(truth), "장,", round(tot / 1e6, 2), "MB →", out)


if __name__ == "__main__":
    main()
