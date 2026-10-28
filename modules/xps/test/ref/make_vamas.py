"""Synthetic (NOT measured) VAMAS ISO 14976 files for the series / ARXPS / map readers and engines.
    python3 -I test/ref/make_vamas.py      (from modules/xps)  ->  test/data/synthetic/*
Every file has a *_truth.json with the generating parameters. Block layout follows ISO 14976 (Dench, Hazell & Seah,
Surf. Interface Anal. 13 (1988) 63) in the order the reader (src/xps_io.js vamasParse) expects.
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "synthetic")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(14976)
HV = 1486.6


def gl(E, c, w, m=0.3):
    x = (E - c) / w
    return (1 - m) * np.exp(-4 * np.log(2) * x * x) + m / (1 + 4 * x * x)


def area1(w, m=0.3):
    """∫ gl dE over ℝ = (1−m)·w·sqrt(π/(4 ln2)) + m·π w/2"""
    return (1 - m) * w * np.sqrt(np.pi / (4 * np.log(2))) + m * np.pi * w / 2


def shirley_spectrum(E, peaks, b0, k):
    """peaks [(center, width, area)], background b0 + k·∫_{Emin}^{E} P (exact Shirley form)"""
    P = sum(A / area1(w, 0) * gl(E, c, w, 0) for c, w, A in peaks)   # Gaussian peaks: the whole area lies inside the windows
    Q = np.concatenate([[0], np.cumsum((P[1:] + P[:-1]) * np.diff(E) / 2)])
    return P + b0 + k * Q


class W:
    def __init__(self): self.L = []
    def __call__(self, *v): self.L.extend(str(x) for x in v)


def header(w, mode, nreg, expvars, nblocks, comment, map_dims=None):
    w("VAMAS Surface Chemical Analysis Standard Data Transfer Format 1988 May 4")
    w("INNO synthetic", "generator make_vamas.py", "nobody", "synthetic " + mode)
    w(1, comment)
    w(mode, "REGULAR")
    if mode in ("MAP", "MAPDP", "NORM", "SDP"): w(nreg)
    if mode in ("MAP", "MAPDP"): w(*map_dims)
    w(len(expvars))
    for lab, unit in expvars: w(lab, unit)
    w(0, 0, 0, 0)            # inclusion list, manual items, future exp, future block
    w(nblocks)


def block(w, mode, bid, sample, expvals, species, trans, be, counts, takeoff=0.0, xy=None, sputter=None, dwell=0.1, scans=10):
    w(bid, sample, 2026, 10, 28, 12, 0, 0, 9, 0, "XPS")
    if mode in ("MAP", "MAPDP"): w(*xy)
    for v in expvals: w(v)
    w("Al")
    if mode in ("MAPDP", "MAPSVDP", "SDP", "SDPSV"): w(18, 1, 1)           # Ar+
    w(HV, 300, 400, 400)
    if mode in ("MAP", "MAPDP", "MAPSV", "MAPSVDP", "SEM"): w(4000, 3000)   # field of view
    w(54.7, 0, "FAT", 20, 1, 4.5, 0, 400, 400, takeoff, 0, species, trans, -1)
    ke = HV - be[::-1]                 # store on an ascending KE axis (BE descending), like most exporters
    y = counts[::-1]
    w("kinetic energy", "eV", f"{ke[0]:.6f}", f"{ke[1] - ke[0]:.6f}", 1, "counts", "d")
    w("pulse counting", dwell, scans, 0)
    if mode in ("MAPDP", "MAPSVDP", "SDP", "SDPSV"): w(*sputter)
    w(0, 0, 0, 0)                       # tilt, tilt azimuth, rotation, n additional params
    w(len(y), f"{y.min():.6g}", f"{y.max():.6g}")
    for v in y: w(f"{v:.6g}")


# ---------------------------------------------------------------- 1) sputter depth profile (SDP): SiO2 / Si
levels = [0, 30, 60, 90, 120, 150, 180, 240, 300, 360, 420, 480]
Si_S, O_S = 0.865, 2.93      # sensitivity factors used by the generator (arbitrary but fixed; truth at% uses them)
Ebe_si = np.arange(94.0, 110.0 + 1e-9, 0.05)
Ebe_o = np.arange(524.0, 540.0 + 1e-9, 0.05)
cfac = 0.1 * 10
blocks = []
truth_at = []
for t in levels:
    f = 1 / (1 + np.exp((t - 200) / 35))           # oxide fraction vs etch time
    nSi, nO = 1.0, 2.0 * f                          # SiO2 → Si
    tot = nSi + nO
    truth_at.append([100 * nSi / tot, 100 * nO / tot])
    scale = 20000
    A_si4, A_si0 = scale * Si_S * f, scale * Si_S * (1 - f)
    A_o = scale * O_S * nO
    si = shirley_spectrum(Ebe_si, [(99.3, 0.9, A_si0), (103.4, 1.5, A_si4)], 400, 0.02)
    o = shirley_spectrum(Ebe_o, [(532.8, 1.6, A_o)], 900, 0.03) if A_o > 1e-9 else np.full_like(Ebe_o, 900.0)
    blocks.append(("Si", "2p", t, Ebe_si, rng.poisson(si * cfac).astype(float)))
    blocks.append(("O", "1s", t, Ebe_o, rng.poisson(o * cfac).astype(float)))
w = W()
header(w, "SDP", 2, [("Etch Time", "s")], len(blocks), "synthetic SiO2/Si sputter depth profile")
for i, (sp, tr, t, be, cnt) in enumerate(blocks):
    block(w, "SDP", f"{sp} {tr} level {levels.index(t)}", "SiO2 on Si (synthetic)", [t], sp, tr, be, cnt,
          sputter=(3000, 1e-6, 2000, 2000, 45, 0, "continuous"))
w("end of experiment")
open(os.path.join(OUT, "depth_profile_sdp.vms"), "w").write("\n".join(w.L) + "\n")
json.dump(dict(regions=["Si 2p", "O 1s"], levels=levels, at=truth_at, sputterEnergy=3000,
               regionInfo={"Si 2p": dict(el="Si", S=Si_S, BE=100, lo=95.0, hi=108.0), "O 1s": dict(el="O", S=O_S, BE=532.8, lo=526.0, hi=539.0)}),
          open(os.path.join(OUT, "depth_profile_truth.json"), "w"), indent=1)

# multi-column text version of the Si 2p series (energy + 12 spectra)
with open(os.path.join(OUT, "depth_profile_multicol.txt"), "w") as fh:
    fh.write("Synthetic Si 2p depth profile (not measured data)\n")
    fh.write("Binding Energy (eV)\t" + "\t".join(f"Etch Time {t} s" for t in levels) + "\n")
    Y = [b[4] / cfac for b in blocks if b[0] == "Si"]
    for i in range(len(Ebe_si) - 1, -1, -1):      # descending BE as Avantage exports
        fh.write(f"{Ebe_si[i]:.3f}\t" + "\t".join(f"{y[i]:.4f}" for y in Y) + "\n")

# ---------------------------------------------------------------- 2) ARXPS angle series: SiO2 overlayer d on Si
d_true, Lo, Ls = 1.8, 3.4, 3.1           # nm; attenuation lengths of Si4+ and Si0 electrons in SiO2
Rinf = 0.9
angles = [0, 15, 30, 40, 50, 55, 60, 65]
Eb = np.arange(94.0, 110.0 + 1e-9, 0.05)
ablocks = []
truthI = []
for th in angles:
    c = np.cos(np.radians(th))
    Io = Rinf * (1 - np.exp(-d_true / (Lo * c)))
    Is = np.exp(-d_true / (Ls * c))
    G = 30000 * (0.6 + 0.4 * c)               # angle-dependent instrument factor (unknown to the analysis)
    truthI.append([G * Io, G * Is])
    sp = shirley_spectrum(Eb, [(99.3, 0.9, G * Is), (103.5, 1.5, G * Io)], 300, 0.015)
    ablocks.append((th, rng.poisson(sp * cfac).astype(float)))
w = W()
header(w, "NORM", 1, [("Angle", "degrees")], len(ablocks), "synthetic ARXPS SiO2/Si")
for th, cnt in ablocks:
    block(w, "NORM", f"Si 2p {th} deg", "SiO2 1.8 nm on Si (synthetic)", [th], "Si", "2p", Eb, cnt, takeoff=th)
w("end of experiment")
open(os.path.join(OUT, "arxps_angles.vms"), "w").write("\n".join(w.L) + "\n")
json.dump(dict(d=d_true, Lo=Lo, Ls=Ls, Rinf=Rinf, angles=angles, I=truthI, peaks=dict(Si0=[99.3, 0.9], Si4=[103.5, 1.5]), cfac=cfac),
          open(os.path.join(OUT, "arxps_truth.json"), "w"), indent=1)

# ---------------------------------------------------------------- 3) MAP mode: C 1s chemical-state map 8 x 6
nx, ny = 8, 6
Ec = np.arange(280.0, 294.0 + 1e-9, 0.1)
mblocks = []
mapA, mapB = np.zeros((ny, nx)), np.zeros((ny, nx))
for iy in range(ny):
    for ix in range(nx):
        a = 4000 * (ix + 1) / nx
        b = 3000 * (1 if (ix - 3.5) ** 2 + (iy - 2.5) ** 2 < 4 else 0.15)
        mapA[iy, ix], mapB[iy, ix] = a, b
        sp = shirley_spectrum(Ec, [(284.8, 1.1, a), (288.6, 1.3, b)], 200, 0.01)
        mblocks.append(((ix * 50.0, iy * 50.0), rng.poisson(sp * cfac).astype(float)))
w = W()
header(w, "MAP", 1, [], len(mblocks), "synthetic C 1s map", map_dims=(nx * ny, nx, ny))
for (x, y), cnt in mblocks:
    block(w, "MAP", f"C 1s x{x:g} y{y:g}", "patterned C (synthetic)", [], "C", "1s", Ec, cnt, xy=(x, y))
w("end of experiment")
open(os.path.join(OUT, "map_c1s.vms"), "w").write("\n".join(w.L) + "\n")
json.dump(dict(nx=nx, ny=ny, areaA=mapA.ravel().tolist(), areaB=mapB.ravel().tolist(), peaks=dict(A=[284.8, 1.1], B=[288.6, 1.3]), cfac=cfac,
               shirley=dict(b0=200, k=0.01)), open(os.path.join(OUT, "map_truth.json"), "w"), indent=1)
print("wrote", sorted(os.listdir(OUT)))
