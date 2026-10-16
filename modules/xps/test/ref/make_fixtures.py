"""
make_fixtures.py -- writes JSON fixtures from the independent Python reference (xps_ref.py).

usage:  python3 make_fixtures.py [lineshapes|backgrounds|fits|quant|all]
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy import integrate

import xps_ref as X

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "fixtures")
DATA = os.path.normpath(os.path.join(HERE, "..", "data", "bct600"))
PREFIX = "_D.I._on_BCT__600_C_5h___"
os.makedirs(OUT, exist_ok=True)


def dump(name, obj):
    with open(os.path.join(OUT, name), "w") as fh:
        json.dump(obj, fh, indent=1, allow_nan=False)
    print("wrote", name)


def r17(a):
    return [float(v) for v in np.asarray(a, float)]


# ---------------------------------------------------------------------------
# lineshapes
# ---------------------------------------------------------------------------
LS_CASES = [
    ("GLS", dict(F=0.5, m=0.0)), ("GLS", dict(F=1.0, m=0.3)), ("GLS", dict(F=2.0, m=1.0)),
    ("GLS", dict(F=1.2, m=0.5)), ("GLS", dict(F=3.0, m=0.1)),
    ("GLP", dict(F=0.5, m=0.0)), ("GLP", dict(F=1.0, m=0.3)), ("GLP", dict(F=1.0, m=1.0)),
    ("GLP", dict(F=2.0, m=0.5)), ("GLP", dict(F=1.5, m=0.9)), ("GLP", dict(F=0.8, m=0.05)),
    ("V", dict(wG=1.0, wL=0.0)), ("V", dict(wG=0.0, wL=1.0)), ("V", dict(wG=1.0, wL=0.3)),
    ("V", dict(wG=0.3, wL=1.0)), ("V", dict(wG=0.7, wL=0.7)), ("V", dict(wG=2.0, wL=0.1)),
    ("V", dict(wG=0.05, wL=1.5)),
    ("DS", dict(alpha=0.0, wL=0.5, wG=0.6)), ("DS", dict(alpha=0.1, wL=0.3, wG=0.0)),
    ("DS", dict(alpha=0.1, wL=0.3, wG=0.6)), ("DS", dict(alpha=0.3, wL=0.2, wG=0.5)),
    ("DS", dict(alpha=0.3, wL=0.6, wG=1.0)), ("DS", dict(alpha=0.45, wL=0.4, wG=0.4)),
    ("DS", dict(alpha=0.05, wL=0.1, wG=0.8)), ("DS", dict(alpha=0.2, wL=1.0, wG=0.3)),
    ("LA", dict(F=1.0, alpha=1.0, beta=1.0, wG=0.0)), ("LA", dict(F=1.0, alpha=1.0, beta=1.0, wG=0.5)),
    ("LA", dict(F=1.2, alpha=2.5, beta=1.2, wG=0.0)), ("LA", dict(F=1.2, alpha=2.5, beta=1.2, wG=0.4)),
    ("LA", dict(F=0.8, alpha=1.1, beta=2.2, wG=0.6)), ("LA", dict(F=0.6, alpha=0.6, beta=5.0, wG=0.3)),
    ("LA", dict(F=2.0, alpha=20.0, beta=0.8, wG=1.0)), ("LA", dict(F=1.0, alpha=0.55, beta=0.55, wG=0.2)),
    ("GLT", dict(F=1.0, m=0.3, T=0.0, tau=0.5)), ("GLT", dict(F=1.0, m=0.3, T=0.3, tau=0.5)),
    ("GLT", dict(F=1.0, m=0.3, T=1.0, tau=0.5)), ("GLT", dict(F=1.2, m=0.0, T=0.3, tau=1.5)),
    ("GLT", dict(F=0.8, m=1.0, T=0.3, tau=0.2)), ("GLT", dict(F=1.5, m=0.5, T=1.0, tau=3.0)),
    ("GLT", dict(F=1.0, m=0.2, T=0.5, tau=0.05)),
]


def width_scale(shape, p):
    if shape == "V":
        return X.voigt_fwhm_OL(p["wG"], p["wL"])
    if shape == "DS":
        return X.voigt_fwhm_OL(p["wG"], p["wL"])
    if shape == "LA":
        return math.hypot(p["F"], p["wG"])
    return p["F"]


def make_lineshapes():
    cases = []
    for i, (shape, p) in enumerate(LS_CASES):
        W = width_scale(shape, p)
        region = None
        if shape in ("GLS", "GLP", "V"):
            x = np.linspace(-4 * W, 4 * W, 201)
        elif shape == "DS":
            region = [-5 * W, 10 * W]
            x = np.linspace(region[0], region[1], 201)
        elif shape == "LA":
            x = np.linspace(-6 * W, 6 * W, 201)
        else:  # GLT
            x = np.linspace(-4 * W, 4 * W + 6 * p["tau"] * (p["T"] > 0), 201)
        meth = "quad"
        if shape == "DS":
            f = X.ds(x, p["alpha"], p["wL"], p["wG"], region, method="quad")
        elif shape == "LA":
            f = X.la(x, p["F"], p["alpha"], p["beta"], p["wG"], method="quad")
        else:
            f = X.shape_eval(shape, x, p)
        # scalar function for FWHM / area (fine-grid conv is accurate to <1e-9 vs quad)
        if shape == "DS":
            fs = lambda v: X.ds(np.array([v]), p["alpha"], p["wL"], p["wG"], region, method="fine")[0]
        else:
            fs = lambda v: X.shape_eval(shape, np.array([v]), p, region)[0]
        fw, xm, fm = X.numeric_fwhm(fs, 0.0, W)
        area = None
        if not (shape == "DS" and p["alpha"] > 0):
            if shape == "LA" and p["wG"] > 0:
                area = 1.0  # analytic (normalised b/N convolved with unit Gaussian); checked in validate.py
                aq, _ = integrate.quad(fs, -np.inf, np.inf, limit=800, epsabs=0, epsrel=1e-10)
                area = aq
            else:
                area, _ = integrate.quad(fs, -np.inf, np.inf, limit=800, epsabs=0, epsrel=1e-12)
        c = dict(id=f"ls{i:02d}_{shape}", shape=shape, params=p, E0=0.0, x=r17(x), f=r17(f),
                 area_R=area, fwhm_true=fw, x_max=xm, f_max=fm, region=region)
        if shape == "DS":
            c["norm"] = "R (alpha=0, literal reading)" if p["alpha"] == 0 else "region"
            if p["alpha"] == 0:
                c["f_region_norm"] = r17(X.ds(x, 0.0, p["wL"], p["wG"], region, method="quad",
                                              alpha0_norm="region"))
            c["region_integral_unnormalised"] = X.ds_region_integral(region[0], region[1], p["alpha"],
                                                                     p["wL"], p["wG"], method="quad")
        if shape == "V":
            c["fwhm_OL"] = X.voigt_fwhm_OL(p["wG"], p["wL"])
        if shape == "GLP":
            c["J"] = X.glp_J(p["m"])
        if shape == "LA":
            c["N"] = X.la_N(p["F"], p["alpha"], p["beta"])
        cases.append(c)
    dump("lineshapes.json", dict(
        spec="SPEC.md v1 sec.1", generator="xps_ref.py (scipy wofz / quad / closed forms)",
        note="x is E-E0 (E0=0). f is unit-area over R, except DS: alpha>0 normalised over region "
             "(region given relative to E0); DS alpha=0 normalised over R (literal reading), "
             "with f_region_norm giving the other reading.",
        cases=cases))


# ---------------------------------------------------------------------------
# backgrounds
# ---------------------------------------------------------------------------
BG_REGIONS = [
    ("C1s", 281.0, 292.5), ("C1s", 280.0, 296.0),
    ("O1s", 526.0, 535.5),
    ("Co2p", 772.0, 808.0), ("Co2p", 773.0, 789.0),
    ("Co3p", 55.5, 66.0),
    ("Ta4f", 21.0, 30.5),
    ("Ba4d", 84.5, 94.5),
]
BG_TYPES = [
    ("linear", {}), ("shirley", {}), ("smart", {}),
    ("tougaard", dict(kind="U2", C=1643.0)), ("tougaard", dict(kind="U3", C=551.0, D=436.0)),
]


def load(core):
    return X.load_txt(os.path.join(DATA, PREFIX + core + ".txt"))


def bg_case(cid, file, E, y, typ, prm, n_end, expected=None, note=None):
    spec = dict(type=typ, n_end=n_end, params=prm)
    info = {}
    if typ in ("shirley", "smart"):
        B, info = X.bg_shirley(E, y, n_end, smart=(typ == "smart"), return_info=True)
    elif typ == "tougaard":
        B, info = X.bg_tougaard(E, y, n_end, prm["kind"], prm["C"], prm.get("D", 0.0), return_info=True)
    else:
        B = X.background(E, y, spec)
    c = dict(id=cid, file=file, type=typ, n_end=n_end, params=prm, E=r17(E), y=r17(y), B=r17(B),
             endpoints=X.endpoints(E, y, n_end), info=info)
    if typ == "smart":
        c["ys"] = r17(X.smooth_quadratic(E, y))
        c["smooth_h"] = X.smooth_h(E)
    if expected is not None:
        c["B_expected"] = r17(expected)
        c["max_abs_dev_from_expected"] = float(np.max(np.abs(B - expected)))
    if note:
        c["note"] = note
    return c


def make_backgrounds():
    cases = []
    for core, a, b in BG_REGIONS:
        E, y = load(core)
        Er, yr = X.region_slice(E, y, a, b)
        for typ, prm in BG_TYPES:
            tag = typ if typ != "tougaard" else prm["kind"]
            cases.append(bg_case(f"{core}_{a:g}-{b:g}_{tag}", PREFIX + core + ".txt", Er, yr, typ, prm, 5))
    # synthetic 1: flat + symmetric peak, equal ends -> all backgrounds constant = y_lo
    E = np.round(np.arange(0.0, 20.0001, 0.1), 10)
    P = 4000 * X.gauss(E - 10.0, 1.5)
    y = 500.0 + P
    for typ, prm in BG_TYPES:
        tag = typ if typ != "tougaard" else prm["kind"]
        cases.append(bg_case(f"syn_flat_{tag}", None, E, y, typ, prm, 5,
                             expected=np.full_like(E, 500.0 + P[0] * 0 + np.mean(P[:5])),
                             note="y_lo == y_hi -> B constant = y_lo (shirley rule; linear trivially; "
                                  "tougaard: s=0; smart: min(y_lo, ys))"))
    # synthetic 2: Gaussian on a self-consistent (discrete) Shirley background, n_end=1
    E = np.round(np.arange(0.0, 20.0001, 0.05), 10)
    P = 10000 * X.gauss(E - 9.0, 1.2)
    Bt = 300.0 + 0.08 * X.cumtrapz(P, E)
    y = Bt + P
    cases.append(bg_case("syn_shirley_gauss", None, E, y, "shirley", {}, 1, expected=Bt,
                         note="y = 300 + 0.08*cumtrapz(P) + P: fixed point of the discrete Shirley "
                              "iteration (n_end=1, P(ends)~1e-30)"))
    # synthetic 3: linear ramp + Gaussian -> linear recovers the ramp (n_end=5)
    E = np.round(np.arange(100.0, 115.0001, 0.1), 10)
    ramp = 1000.0 + 25.0 * (E - 100.0)
    y = ramp + 5000 * X.gauss(E - 107.5, 1.0)
    cases.append(bg_case("syn_linear_ramp", None, E, y, "linear", {}, 5, expected=ramp,
                         note="mean of 5 collinear points lies on the line; Gaussian ~0 at ends"))
    dump("backgrounds.json", dict(spec="SPEC.md v1 sec.2", cases=cases))


# ---------------------------------------------------------------------------
# fits
# ---------------------------------------------------------------------------

def F_(v, lo=None, hi=None):
    d = dict(mode="free", value=v)
    if lo is not None:
        d["lo"] = lo
    if hi is not None:
        d["hi"] = hi
    return d


def FX(v):
    return dict(mode="fixed", value=v)


def LK(ref, k=1.0, b=0.0):
    return dict(mode="link", ref=ref, k=k, b=b)


def synth(model, truth_theta, E, c, seed, B_truth=None):
    """Poisson data: counts ~ Poisson(c * mu_true); y = counts / c (CPS)."""
    m = X.Model(model, E, np.ones_like(E), c=c)
    mu, B, parts, vals = m.mu(np.asarray(truth_theta, float), return_parts=True)
    if B_truth is not None:
        mu = mu - B + B_truth
        B = B_truth
    rng = np.random.default_rng(seed)
    y = rng.poisson(c * mu) / c
    return y, mu, B, vals


def truth_dict(vals):
    out = {}
    for (own, pn), v in vals.items():
        out.setdefault(own, {})[pn] = v
    return out


def run_fit(cid, kind, model, E, y, c, file=None, truth=None, n_starts=24, extra=None):
    t0 = time.time()
    m = X.Model(model, E, y, c=c)
    res = X.fit(m, n_starts=n_starts)
    res["seconds"] = time.time() - t0
    res["kappa_noise"] = X.noise_kappa(E, y)
    case = dict(id=cid, kind=kind, file=file, region=[float(E[0]), float(E[-1])], c=c,
                E=r17(E), y=r17(y), model=model, result=res)
    if truth is not None:
        case["truth"] = truth
    if extra:
        case.update(extra)
    print(f"  {cid}: chi2={res['chi2']:.3f} phi={res['phi']:.4f} N={res['N']} p={res['p']} "
          f"rho1={res['rho1']:.3f} uniq={res['uniqueness']:.2f} at_bound={res['at_bound']} "
          f"({res['seconds']:.1f}s)")
    return case


def synthetic_fits():
    cases = []
    # S1 single GLS on linear background --------------------------------------
    E = np.round(np.arange(280.0, 290.0001, 0.1), 10)
    lin = 200.0 + 5.0 * (E - 280.0)
    model = dict(background=dict(type="linear", n_end=5), components=[
        dict(name="p1", shape="GLS", E0=F_(285.0, 283.0, 287.0), A=F_(3000.0),
             F=F_(1.0, 0.3, 3.0), m=F_(0.5, 0.0, 1.0))])
    gen = dict(background=dict(type="none"), components=[
        dict(name="p1", shape="GLS", E0=FX(285.1), A=FX(5000.0), F=FX(1.2), m=FX(0.3))])
    y, mu, B, vals = synth(gen, [], E, 1.0, 1, B_truth=lin)
    cases.append(run_fit("syn1_GLS_linear", "synthetic", model, E, y, 1.0,
                         truth=dict(params=truth_dict(vals), background="200+5*(E-280)", mu=r17(mu), seed=1)))
    # S2 C1s-like 3 GLS on Shirley (generated as active-Shirley truth), fitted with static shirley
    E = np.round(np.arange(280.0, 294.0001, 0.05), 10)
    gen = dict(background=dict(type="active", c0=FX(2000.0), k=FX(0.02), s=FX(0.0)), components=[
        dict(name="CC", shape="GLS", E0=FX(284.8), A=FX(12000.0), F=FX(1.1), m=FX(0.3)),
        dict(name="CO", shape="GLS", E0=FX(286.3), A=FX(4000.0), F=FX(1.1), m=FX(0.3)),
        dict(name="OCO", shape="GLS", E0=FX(288.6), A=FX(2500.0), F=FX(1.1), m=FX(0.3))])
    c = 2.0
    y, mu, B, vals = synth(gen, [], E, c, 2)
    model = dict(background=dict(type="shirley", n_end=5), components=[
        dict(name="CC", shape="GLS", E0=F_(284.7, 284.0, 285.5), A=F_(10000.0), F=F_(1.0, 0.5, 2.5),
             m=FX(0.3)),
        dict(name="CO", shape="GLS", E0=F_(286.2, 285.6, 287.2), A=F_(3000.0), F=LK("CC"), m=FX(0.3)),
        dict(name="OCO", shape="GLS", E0=F_(288.5, 287.6, 289.5), A=F_(2000.0), F=LK("CC"), m=FX(0.3))])
    cases.append(run_fit("syn2_C1s_3GLS_shirley", "synthetic", model, E, y, c,
                         truth=dict(params=truth_dict(vals), background="active-Shirley c0=2000,k=0.02",
                                    B=r17(B), mu=r17(mu), seed=2)))
    # S3 Voigt doublet Co2p-like, delta=15 fixed, r=0.5 fixed, q free --------------------
    E = np.round(np.arange(770.0, 805.0001, 0.1), 10)
    gen = dict(background=dict(type="active", c0=FX(30000.0), k=FX(0.15), s=FX(0.0)), components=[
        dict(name="Co2p", shape="V", E0=FX(780.2), A=FX(60000.0), wG=FX(1.0), wL=FX(1.4),
             doublet=dict(delta=FX(15.0), r=FX(0.5), q=FX(1.15)))])
    c = 1.0
    y, mu, B, vals = synth(gen, [], E, c, 3)
    model = dict(background=dict(type="active", c0=F_(29000.0), k=F_(0.1, 0.0, 2.0), s=FX(0.0)), components=[
        dict(name="Co2p", shape="V", E0=F_(780.0, 779.0, 781.5), A=F_(50000.0), wG=F_(1.2, 0.2, 3.0),
             wL=F_(1.0, 0.05, 4.0), doublet=dict(delta=FX(15.0), r=FX(0.5), q=F_(1.0, 0.7, 1.6)))])
    cases.append(run_fit("syn3_Co2p_Voigt_doublet", "synthetic", model, E, y, c,
                         truth=dict(params=truth_dict(vals), B=r17(B), mu=r17(mu), seed=3)))
    # S4 DS metal + GLS oxide on Shirley -------------------------------------------------
    E = np.round(np.arange(66.0, 80.0001, 0.05), 10)
    gen = dict(background=dict(type="active", c0=FX(1500.0), k=FX(0.05), s=FX(0.0)), components=[
        dict(name="metal", shape="DS", E0=FX(71.0), A=FX(8000.0), alpha=FX(0.12), wL=FX(0.35), wG=FX(0.6)),
        dict(name="oxide", shape="GLS", E0=FX(73.1), A=FX(3000.0), F=FX(1.5), m=FX(0.3))])
    c = 1.0
    y, mu, B, vals = synth(gen, [], E, c, 4)
    model = dict(background=dict(type="active", c0=F_(1400.0), k=F_(0.03, 0.0, 2.0), s=FX(0.0)), components=[
        dict(name="metal", shape="DS", E0=F_(71.0, 70.3, 71.7), A=F_(7000.0), alpha=F_(0.1, 0.0, 0.45),
             wL=F_(0.4, 0.05, 1.5), wG=F_(0.6, 0.1, 1.5)),
        dict(name="oxide", shape="GLS", E0=F_(73.0, 72.2, 74.0), A=F_(2500.0), F=F_(1.4, 0.6, 3.0), m=FX(0.3))])
    cases.append(run_fit("syn4_DS_metal_GLS_oxide", "synthetic", model, E, y, c,
                         truth=dict(params=truth_dict(vals), B=r17(B), mu=r17(mu), seed=4,
                                    note="DS area is the area inside the fit region [E[0],E[-1]]"),
                         n_starts=12))
    # S5 LA asymmetric on linear ---------------------------------------------------------
    E = np.round(np.arange(845.0, 860.0001, 0.05), 10)
    lin = 4000.0 + 30.0 * (E - 845.0)
    gen = dict(background=dict(type="none"), components=[
        dict(name="Ni", shape="LA", E0=FX(852.6), A=FX(20000.0), F=FX(0.9), alpha=FX(2.5), beta=FX(1.15),
             wG=FX(0.5))])
    y, mu, B, vals = synth(gen, [], E, 1.0, 5, B_truth=lin)
    model = dict(background=dict(type="linear", n_end=5), components=[
        dict(name="Ni", shape="LA", E0=F_(852.5, 851.5, 853.5), A=F_(15000.0), F=F_(1.0, 0.3, 2.5),
             alpha=F_(2.0, 0.55, 20.0), beta=F_(1.3, 0.55, 20.0), wG=F_(0.4, 0.05, 1.5))])
    cases.append(run_fit("syn5_LA_asym_linear", "synthetic", model, E, y, 1.0,
                         truth=dict(params=truth_dict(vals), background="4000+30*(E-845)", mu=r17(mu), seed=5),
                         n_starts=12))
    # S6 active Shirley with 2 peaks -----------------------------------------------------
    E = np.round(np.arange(525.0, 538.0001, 0.05), 10)
    gen = dict(background=dict(type="active", c0=FX(5000.0), k=FX(0.06), s=FX(0.0)), components=[
        dict(name="lattice", shape="GLS", E0=FX(529.6), A=FX(15000.0), F=FX(1.2), m=FX(0.25)),
        dict(name="OH", shape="GLS", E0=FX(531.3), A=FX(9000.0), F=FX(1.6), m=FX(0.25))])
    c = 0.5
    y, mu, B, vals = synth(gen, [], E, c, 6)
    model = dict(background=dict(type="active", c0=F_(4800.0), k=F_(0.04, 0.0, 1.0), s=FX(0.0)), components=[
        dict(name="lattice", shape="GLS", E0=F_(529.5, 528.8, 530.4), A=F_(12000.0), F=F_(1.0, 0.5, 2.5),
             m=FX(0.25)),
        dict(name="OH", shape="GLS", E0=F_(531.2, 530.5, 532.5), A=F_(8000.0), F=F_(1.4, 0.5, 3.0), m=FX(0.25))])
    cases.append(run_fit("syn6_active_shirley_2GLS", "synthetic", model, E, y, c,
                         truth=dict(params=truth_dict(vals), B=r17(B), mu=r17(mu), seed=6)))
    # S7 GLP + GLT on smart background (Tougaard-free check of the other shapes) ---------
    E = np.round(np.arange(20.0, 32.0001, 0.05), 10)
    gen = dict(background=dict(type="active", c0=FX(800.0), k=FX(0.03), s=FX(0.0)), components=[
        dict(name="a", shape="GLP", E0=FX(24.0), A=FX(5000.0), F=FX(0.9), m=FX(0.4)),
        dict(name="b", shape="GLT", E0=FX(26.5), A=FX(4000.0), F=FX(1.0), m=FX(0.2), T=FX(0.35), tau=FX(1.2))])
    y, mu, B, vals = synth(gen, [], E, 1.0, 7)
    model = dict(background=dict(type="smart", n_end=5), components=[
        dict(name="a", shape="GLP", E0=F_(24.1, 23.3, 24.8), A=F_(4000.0), F=F_(1.0, 0.4, 2.0), m=F_(0.3, 0.0, 1.0)),
        dict(name="b", shape="GLT", E0=F_(26.4, 25.7, 27.2), A=F_(3000.0), F=F_(1.1, 0.4, 2.0), m=FX(0.2),
             T=F_(0.3, 0.0, 1.0), tau=F_(1.0, 0.1, 5.0))])
    cases.append(run_fit("syn7_GLP_GLT_smart", "synthetic", model, E, y, 1.0,
                         truth=dict(params=truth_dict(vals), B=r17(B), mu=r17(mu), seed=7)))
    return cases


def real_fits():
    cases = []
    # R1 C1s ------------------------------------------------------------------------
    E, y = load("C1s")
    Er, yr = X.region_slice(E, y, 280.5, 293.0)
    model = dict(background=dict(type="smart", n_end=5), components=[
        dict(name="CC", shape="GLS", E0=F_(284.6, 284.0, 285.2), A=F_(15000.0), F=F_(1.2, 0.6, 2.5), m=FX(0.3)),
        dict(name="CO", shape="GLS", E0=F_(286.0, 285.4, 286.8), A=F_(3000.0), F=LK("CC"), m=FX(0.3)),
        dict(name="CdO", shape="GLS", E0=F_(287.8, 287.0, 288.3), A=F_(1500.0), F=LK("CC"), m=FX(0.3)),
        dict(name="CO3", shape="GLS", E0=F_(288.9, 288.4, 289.6), A=F_(5000.0), F=F_(1.2, 0.6, 2.5), m=FX(0.3))])
    cases.append(run_fit("real1_C1s_4GLS_smart", "real", model, Er, yr, 1.0, file=PREFIX + "C1s.txt"))
    # R2 O1s ------------------------------------------------------------------------
    E, y = load("O1s")
    Er, yr = X.region_slice(E, y, 525.5, 536.0)
    model = dict(background=dict(type="shirley", n_end=5), components=[
        dict(name="O_lat", shape="GLS", E0=F_(528.6, 527.8, 529.4), A=F_(100000.0), F=F_(1.3, 0.6, 2.5), m=FX(0.1)),
        dict(name="O_2", shape="GLS", E0=F_(530.8, 530.0, 531.5), A=F_(80000.0), F=F_(1.5, 0.6, 3.0), m=FX(0.1)),
        dict(name="O_3", shape="GLS", E0=F_(532.3, 531.6, 533.5), A=F_(20000.0), F=F_(1.6, 0.6, 3.0), m=FX(0.1))])
    cases.append(run_fit("real2_O1s_3GLS_shirley", "real", model, Er, yr, 1.0, file=PREFIX + "O1s.txt"))
    # R3 Ta4f -----------------------------------------------------------------------
    E, y = load("Ta4f")
    Er, yr = X.region_slice(E, y, 21.5, 30.5)
    model = dict(background=dict(type="shirley", n_end=5), components=[
        dict(name="Ta4f", shape="GLS", E0=F_(24.2, 23.6, 24.8), A=F_(40000.0), F=F_(1.0, 0.4, 2.0),
             m=F_(0.3, 0.0, 1.0), doublet=dict(delta=F_(1.9, 1.75, 2.05), r=FX(0.75), q=FX(1.0))),
        dict(name="Ta4f_b", shape="GLS", E0=F_(25.6, 25.0, 26.6), A=F_(3000.0), F=LK("Ta4f"), m=LK("Ta4f"),
             doublet=dict(delta=LK("Ta4f"), r=FX(0.75), q=FX(1.0))),
        dict(name="X29", shape="GLS", E0=F_(29.3, 28.3, 30.3), A=F_(2000.0), F=F_(2.0, 0.8, 4.0), m=FX(0.3))])
    cases.append(run_fit("real3_Ta4f_2doublets_GLS_shirley", "real", model, Er, yr, 1.0, file=PREFIX + "Ta4f.txt"))
    return cases


def make_fits(which="all"):
    cases = []
    if which in ("all", "syn"):
        cases += synthetic_fits()
    if which in ("all", "real"):
        cases += real_fits()
    dump("fits.json" if which == "all" else f"fits_{which}.json",
         dict(spec="SPEC.md v1 sec.3", note="see README for schema", cases=cases))


# ---------------------------------------------------------------------------
# quant
# ---------------------------------------------------------------------------
# TPP-2M check materials: Tanuma, Powell, Penn, Surf. Interface Anal. 35 (2003) 268 (IMFP VII),
# Table 1 (M, rho, Nv, Ep, Eg) and Table 4 (IMFPs from optical data, Angstrom).
TPP_VII = [
    # name, M, rho, Nv, Ep_printed, Eg, optical IMFP at 100, 500, 1000, 2000 eV
    ("Li", 6.941, 0.534, 1, 7.99, 0, [7.0, 22.7, 39.7, 71.0]),
    ("Be", 9.0122, 1.848, 2, 18.44, 0, [4.3, 11.7, 19.9, 34.9]),
    ("Diamond", 12.011, 3.515, 4, 31.16, 5.5, [4.7, 9.9, 16.4, 28.1]),
    ("Graphite", 12.011, 2.25, 4, 24.93, 0, [4.4, 10.8, 18.2, 31.6]),
    ("Na", 22.9898, 0.971, 1, 5.92, 0, [7.6, 23.1, 40.0, 71.0]),
    ("Al", 26.9815, 2.699, 3, 15.78, 0, [4.6, 12.7, 21.3, 37.2]),
    ("K", 39.0983, 0.862, 1, 4.28, 0, [10.0, 31.4, 54.8, 97.1]),
    ("Sc", 44.9559, 2.989, 3, 12.86, 0, [4.9, 11.9, 20.2, 35.1]),
    ("Ge", 72.59, 5.32, 4, 15.59, 0.67, [5.0, 12.9, 21.2, 36.5]),
    ("In", 114.818, 7.31, 3, 12.59, 0, [5.4, 12.3, 20.6, 35.3]),
    ("Sn", 118.71, 7.31, 4, 14.29, 0, [6.4, 13.8, 22.8, 38.9]),
    ("Cs", 132.905, 1.88, 1, 3.43, 0, [9.3, 29.3, 50.8, 90.2]),
    ("Gd", 157.25, 8.23, 9, 19.77, 0, [4.6, 11.8, 19.5, 33.7]),
    ("Tb", 158.9253, 8.25, 9, 19.69, 0, [4.2, 10.2, 16.9, 29.2]),
    ("Dy", 162.5, 8.78, 9, 20.08, 0, [4.6, 10.6, 17.5, 30.0]),
]

# a few XPS-relevant materials (inputs are illustrative; Nv convention = valence electrons per formula unit)
TPP_EXTRA = [
    ("Au", 196.967, 19.32, 11, 0.0), ("Si", 28.086, 2.33, 4, 1.1), ("SiO2", 60.08, 2.19, 16, 9.0),
    ("Ta2O5", 441.89, 8.2, 40, 4.0), ("Co3O4", 240.80, 6.11, 44, 1.6), ("BaTiO3", 233.19, 6.02, 24, 3.2),
    ("PMMA(C5H8O2)", 100.12, 1.19, 40, 5.0),
]


def make_quant():
    energies = [100.0, 200.0, 500.0, 1000.0, 1200.0, 1486.6 - 284.8, 2000.0]
    mats = []
    for name, M, rho, Nv, Eg in TPP_EXTRA:
        lam, prm = X.tpp2m(np.array(energies), Nv, rho, M, Eg)
        mats.append(dict(name=name, M=M, rho=rho, Nv=Nv, Eg=Eg, E=energies, lambda_A=r17(lam), derived=prm))
    lit = []
    rows = []
    for name, M, rho, Nv, Ep_pr, Eg, opt in TPP_VII:
        E4 = np.array([100.0, 500.0, 1000.0, 2000.0])
        lam, prm = X.tpp2m(E4, Nv, rho, M, Eg)
        dev = (lam - np.array(opt)) / np.array(opt)
        lit.append(dict(name=name, M=M, rho=rho, Nv=Nv, Eg=Eg, Ep_printed=Ep_pr, Ep_ref=prm["Ep"],
                        E=r17(E4), lambda_TPP2M_ref=r17(lam), lambda_optical_TPP2003=opt, rel_dev=r17(dev)))
        rows.append(dev)
    rows = np.array(rows)
    ep_dev = [abs(l["Ep_ref"] - l["Ep_printed"]) for l in lit]
    rms_by_mat = np.sqrt(np.mean(rows ** 2, axis=1))
    litcheck = dict(
        source="S. Tanuma, C. J. Powell, D. R. Penn, Surf. Interface Anal. 35 (2003) 268-275, "
               "'Calculation of electron IMFPs. VII. Reliability of the TPP-2M IMFP predictive equation' "
               "(NIMS MDR copy). Eqs. (4),(6) identical to SPEC 4; Table 1 parameters; Table 4 optical IMFPs.",
        materials=lit,
        max_abs_Ep_dev_eV=float(max(ep_dev)),
        rms_rel_dev_per_material=dict(zip([l["name"] for l in lit], r17(rms_by_mat))),
        mean_rms_rel_dev=float(np.mean(rms_by_mat)),
        note="Optical-data IMFPs are not TPP-2M outputs; agreement expected only at the ~10-20% level "
             "(the paper reports RMS deviations of this order, large for Li/Na/K/Cs/diamond/graphite). "
             "Exact agreement is expected for Ep (printed to 0.01 eV).")
    # secondary check: published TPP-2M values (nm) quoted in Jaber et al., Jordan J. Phys. 16 (2023) 381,
    # Table 1 column 'Ref. [22]' (Tanuma, Powell, Penn TPP-2M), Si/Cu/Au at 1, 3, 5 keV.
    pub = [("Si", 28.086, 2.33, 4, 1.1, [2.39, 5.79, 8.87]),
           ("Cu", 63.546, 8.96, 11, 0.0, [1.5, 3.6, 5.5]),
           ("Au", 196.967, 19.32, 11, 0.0, [1.23, 2.8, 4.3])]
    sec = []
    for name, M, rho, Nv, Eg, vals in pub:
        lam, _ = X.tpp2m(np.array([1000.0, 3000.0, 5000.0]), Nv, rho, M, Eg)
        sec.append(dict(name=name, M=M, rho=rho, Nv=Nv, Eg=Eg, E=[1000.0, 3000.0, 5000.0],
                        lambda_nm_ref=r17(lam / 10), lambda_nm_published=vals,
                        rel_dev=r17((lam / 10 - np.array(vals)) / np.array(vals))))
    litcheck["published_tpp2m_values"] = dict(
        source="A. M. D. Jaber et al., Jordan J. Phys. 16(4) (2023) 381-393, Table 1, column citing "
               "Tanuma-Powell-Penn (TPP-2M); inputs (Nv, rho, M, Eg) are the standard TPP values assumed "
               "by us (Si Nv=4, Cu/Au Nv=11) since the source does not list them.",
        cases=sec)
    # toy at%
    peaks = [dict(name="C1s", BE=284.8, A=12000.0, sigma_A=300.0, S=1.0),
             dict(name="O1s", BE=530.0, A=52000.0, sigma_A=800.0, S=2.93),
             dict(name="Ta4f", BE=24.5, A=61000.0, sigma_A=900.0, S=8.62),
             dict(name="Co2p3/2", BE=780.0, A=45000.0, sigma_A=1500.0, S=12.6)]
    PE = 50.0
    T_opts = [dict(type="none"), dict(type="power", x=0.5),
              dict(type="thermo", a=[0.0, -0.45, 0.05, -0.004], note="toy coefficients"),
              dict(type="table", KE=[100.0, 400.0, 800.0, 1200.0, 1600.0], T=[1.8, 0.62, 0.36, 0.27, 0.22])]
    L_opts = [dict(type="none"), dict(type="power", n=0.6),
              dict(type="tpp2m", Nv=40, rho=8.2, M=441.89, Eg=4.0, note="Ta2O5-like matrix for all peaks")]
    toy = []
    for To in T_opts:
        for Lo in L_opts:
            res = X.quantify(peaks, PE, To, Lo)
            toy.append(dict(T=To, lam=Lo, PE=PE, hnu=X.HNU_ALKA, result=res))
    dump("quant.json", dict(spec="SPEC.md v1 sec.4", hnu=X.HNU_ALKA, tpp2m_materials=mats,
                            tpp2m_literature_check=litcheck, toy_peaks=peaks, toy_cases=toy))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "lineshapes"):
        make_lineshapes()
    if what in ("all", "backgrounds"):
        make_backgrounds()
    if what in ("all", "quant"):
        make_quant()
    if what in ("all", "fits"):
        make_fits()
    if what in ("syn", "real"):
        make_fits(what)
