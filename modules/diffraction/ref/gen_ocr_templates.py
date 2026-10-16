#!/usr/bin/env python3
"""스케일바 이름표 OCR 글자 틀(template) 만들기 → src/diff_ocr.js 의 /*TPL-BEGIN*/ … /*TPL-END*/ 자리에 넣는다.

글꼴(틀용): Liberation Sans(Arial 치수 호환) · DejaVu Sans · FreeSans(Helvetica 계열), 보통·굵게.
크기(px): 12 15 19 24 30 40 52, 계단(앤티에일리어싱 없음)·부드러움 둘 다.
시험(gen_ocr_test.py)은 다른 크기와 다른 글꼴(Carlito·Inter·Noto Sans CJK 라틴·TeX Gyre Heros)을 쓴다.

틀 하나 = 클래스 1 B + 12×18 비트(27 B) + ln(가로/세로) 1 B + 높이/H 1 B + 바닥/H 1 B + 구멍 수 1 B (32 B).
H = 같은 글꼴·크기 숫자의 중앙 높이, 바닥 = (글자 아래끝+1 − 기준선)/H.
특징 뽑기(norm_bitmap·holes)는 diff_ocr.js 의 feat()·holes() 와 똑같이 맞춘다.
"""
import base64, os, re, sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
JS = os.path.join(HERE, "..", "src", "diff_ocr.js")
GW, GH = 12, 18

FONTS = [
    ("LiberationSans", "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"),
    ("LiberationSansB", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"),
    ("DejaVuSans", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
    ("DejaVuSansB", "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
    ("FreeSans", "/usr/share/fonts/truetype/freefont/FreeSans.ttf"),
    ("FreeSansB", "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf"),
]
SIZES = [12, 15, 19, 24, 30, 40, 52]
# 대상 글자 + 거절용 글자(다른 덧글씨를 이름표로 착각하지 않게)
TARGET = "0123456789./nmµuAÅ-p"
SUP = "⁻¹"            # 위첨자: 작은 '-' '1' 을 올려 그린 것 (+ 글꼴 자체 ¹)
REJECT = "BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijkloqrstvwxyz:()=,+%°'"
CLASSES = TARGET + SUP + REJECT + "μ"   # μ(그리스) 는 µ 로 합친다
SUP_SCALE, SUP_RISE = 0.62, 0.36
# 거의 같은 틀 빼기 기준 (비트 차 · ln 가로세로비 차 · 높이·바닥 차)
DUP_BITS, DUP_LA, DUP_POS = int(os.environ.get("DUP_BITS", 14)), float(os.environ.get("DUP_LA", 0.12)), float(os.environ.get("DUP_POS", 0.08))


def norm_bitmap(ink, x0, y0, x1, y1):
    """bbox [x0,x1]×[y0,y1] (포함) 의 잉크(0..1) 를 GW×GH 로 면적 평균 — JS feat() 와 같은 식"""
    w, h = x1 - x0 + 1, y1 - y0 + 1
    return wmat(GH, h) @ ink[y0:y1 + 1, x0:x1 + 1] @ wmat(GW, w).T


def wmat(n, L):
    """칸 i 가 덮는 [i·L/n, (i+1)·L/n) 와 화소 [p, p+1) 의 겹침 / 칸 폭"""
    s = L / n
    W = np.zeros((n, L))
    for i in range(n):
        a, b = i * s, (i + 1) * s
        for p in range(int(np.floor(a)), min(L, int(np.ceil(b)))):
            W[i, p] = max(0.0, min(b, p + 1) - max(a, p)) / s
    return W


def holes(binm):
    """구멍 수: bbox 안 배경(4-연결) 성분 중 테두리에 닿지 않는 것 (면적 ≥ 1)"""
    h, w = binm.shape
    lab = np.zeros((h, w), np.int32)
    n = 0
    for y in range(h):
        for x in range(w):
            if binm[y, x] or lab[y, x]:
                continue
            n += 1
            st = [(y, x)]; lab[y, x] = n; edge = False
            while st:
                cy, cx = st.pop()
                if cy == 0 or cx == 0 or cy == h - 1 or cx == w - 1:
                    edge = True
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ny, nx = cy + dy, cx + dx
                    if 0 <= ny < h and 0 <= nx < w and not binm[ny, nx] and not lab[ny, nx]:
                        lab[ny, nx] = n; st.append((ny, nx))
            if edge:
                lab[lab == n] = -1
    return len([k for k in range(1, n + 1) if (lab == k).any()])


def render(font_path, size, ch, aa):
    """글자 하나를 기준선(anchor ls) 위에 그려 잉크 배열과 기준선 y 를 돌려준다"""
    S = size * 4
    im = Image.new("L", (S, S), 0)
    d = ImageDraw.Draw(im)
    d.fontmode = "L" if aa else "1"
    base = int(size * 2.6)
    if ch in SUP:
        f = ImageFont.truetype(font_path, max(6, round(size * SUP_SCALE)))
        d.text((size, base - round(size * SUP_RISE)), "-" if ch == "⁻" else "1", font=f, fill=255, anchor="ls")
    else:
        f = ImageFont.truetype(font_path, size)
        d.text((size, base), ch, font=f, fill=255, anchor="ls")
    return np.asarray(im, np.float64) / 255.0, base


def glyph_feat(ink, base, H):
    m = ink >= 0.5
    ys, xs = np.nonzero(m)
    if not len(xs):
        return None
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
    w, h = x1 - x0 + 1, y1 - y0 + 1
    bm = norm_bitmap(ink, x0, y0, x1, y1)
    return dict(bits=(bm >= 0.5).astype(np.uint8).ravel(), la=np.log(w / h), rh=h / H, rb=(y1 + 1 - base) / H,
                ho=holes(m[y0:y1 + 1, x0:x1 + 1]))


def has_glyph(font_path, ch):
    from fontTools.ttLib import TTFont
    try:
        return ord(ch) in TTFont(font_path, fontNumber=0).getBestCmap()
    except Exception:
        return True


def main():
    recs = []
    for fname, fp in FONTS:
        for size in SIZES:
            for aa in (True, False):
                hs = []
                for dch in "0123456789":
                    ink, base = render(fp, size, dch, aa)
                    ys = np.nonzero(ink >= 0.5)[0]
                    hs.append(ys.max() - ys.min() + 1)
                H = float(np.median(hs))
                chars = list(CLASSES)
                for ch in chars:
                    if ch not in SUP and not has_glyph(fp, ch):
                        continue
                    ink, base = render(fp, size, ch, aa)
                    f = glyph_feat(ink, base, H)
                    if f is None:
                        continue
                    f["cls"] = "µ" if ch == "μ" else ch
                    recs.append(f)
                # 글꼴 자체 위첨자 ¹ (있으면)
                if has_glyph(fp, "¹"):
                    S = size * 4; im = Image.new("L", (S, S), 0); dd = ImageDraw.Draw(im); dd.fontmode = "L" if aa else "1"
                    base = int(size * 2.6); dd.text((size, base), "¹", font=ImageFont.truetype(fp, size), fill=255, anchor="ls")
                    f = glyph_feat(np.asarray(im, np.float64) / 255.0, base, H)
                    if f: f["cls"] = "¹"; recs.append(f)
    # 같은 클래스 안 거의 같은 틀 빼기
    cls_list = [c for c in CLASSES if c != "μ"]
    kept = {c: [] for c in cls_list}
    for r in recs:
        dup = False
        for k in kept[r["cls"]]:
            if (np.count_nonzero(k["bits"] != r["bits"]) <= DUP_BITS and abs(k["la"] - r["la"]) < DUP_LA and abs(k["rh"] - r["rh"]) < DUP_POS
                    and abs(k["rb"] - r["rb"]) < DUP_POS and k["ho"] == r["ho"]):
                dup = True; break
        if not dup:
            kept[r["cls"]].append(r)
    out = bytearray()
    n = 0
    for ci, c in enumerate(cls_list):
        for r in kept[c]:
            out.append(ci)
            out += np.packbits(r["bits"]).tobytes()  # 216 비트 = 27 B
            out.append(int(np.clip(round((r["la"] + 4) * 40), 0, 255)))
            out.append(int(np.clip(round(r["rh"] * 100), 0, 255)))
            out.append(int(np.clip(round((r["rb"] + 1) * 100), 0, 255)))
            out.append(int(min(r["ho"], 3)))
            n += 1
    b64 = base64.b64encode(bytes(out)).decode()
    lines = "\n".join('"' + b64[i:i + 160] + '"+' for i in range(0, len(b64), 160)) + '""'
    block = ("/*TPL-BEGIN*/\n// ref/gen_ocr_templates.py 가 만든 글자 틀 %d 개 (글꼴 %s, 크기 %s px, 부드러움 켜고/끄고)\n"
             "const TPL_CLASSES=%s;\nconst TPL_B64=\n%s;\n/*TPL-END*/") % (
        n, " · ".join(f for f, _ in FONTS), " ".join(map(str, SIZES)), _js_str("".join(cls_list)), lines)
    src = open(JS, encoding="utf8").read()
    src2, k = re.subn(r"/\*TPL-BEGIN\*/.*?/\*TPL-END\*/", lambda m: block, src, flags=re.S)
    if k != 1:
        sys.exit("diff_ocr.js 에 TPL 표시가 없다")
    open(JS, "w", encoding="utf8").write(src2)
    per = {c: len(kept[c]) for c in cls_list}
    print("틀", n, "개, ", len(b64), "글자 base64;", per)


def _js_str(s):
    return '"' + "".join(c if ord(c) < 128 and c not in '"\\' else "\\u%04x" % ord(c) for c in s) + '"'


if __name__ == "__main__":
    main()
