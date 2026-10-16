"""
validate.py -- self-validation of xps_ref.py and comparison with third-party code.
Writes fixtures/validation.json and prints a summary.

Third-party (optional, offline pip packages):
  * lmfitxps (shirley_calculate, tougaard_calculate)    pip install lmfitxps
  * pynxtools-xps (models.backgrounds.Shirley)          pip install pynxtools-xps
"""
import json
import math
import os
import io
import contextlib

import numpy as np
from scipy import integrate, optimize

import xps_ref as X

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
R = {}


def rec(key, val, ok=None):
    R[key] = dict(value=val, ok=ok)
    print(f"{'PASS' if ok else ('FAIL' if ok is False else 'INFO')}  {key}: {val}")


# ---------------------------------------------------------------- lineshapes
ls = json.load(open(os.path.join(FIX, "lineshapes.json")))["cases"]
areas = [abs(c["area_R"] - 1) for c in ls if c["area_R"] is not None]
rec("lineshape unit area max|A-1| (quad over R)", max(areas), max(areas) < 1e-9)
fw = [abs(c["fwhm_true"] - c["params"]["F"]) for c in ls if c["shape"] == "GLS"]
fw += [abs(c["fwhm_true"] - c["params"]["F"]) for c in ls if c["shape"] == "GLP" and c["params"]["m"] in (0.0, 1.0)]
fw += [abs(c["fwhm_true"] - c["params"]["F"]) for c in ls if c["shape"] == "GLT" and c["params"]["T"] == 0.0]
rec("FWHM == F for GLS, GLP(m=0,1), GLT(T=0)", max(fw), max(fw) < 1e-9)
ol = [abs(c["fwhm_true"] / c["fwhm_OL"] - 1) for c in ls if c["shape"] == "V"]
rec("Voigt numeric FWHM vs Olivero-Longbothum (rel; spec claims <=2e-4, true worst case ~2.37e-4)",
    max(ol), max(ol) <= 2.4e-4)

x = np.linspace(-5, 5, 401)
lim = {}
lim["V(wG=1e-7)->L"] = np.max(np.abs(X.voigt(x, 1e-7, 0.8) - X.lorentz(x, 0.8))) / np.max(X.lorentz(x, 0.8))
lim["V(wL=1e-7)->G"] = np.max(np.abs(X.voigt(x, 0.8, 1e-7) - X.gauss(x, 0.8))) / np.max(X.gauss(x, 0.8))
lim["LA(1,1,0)->L"] = np.max(np.abs(X.la(x, 0.9, 1, 1, 0) - X.lorentz(x, 0.9))) / np.max(X.lorentz(x, 0.9))
lim["LA(1,1,wG)->V(wG,F)"] = np.max(np.abs(X.la(x, 0.9, 1, 1, 0.7, method="quad") - X.voigt(x, 0.7, 0.9))) / np.max(X.voigt(x, 0.7, 0.9))
lim["DS(a=0)/pi->V"] = np.max(np.abs(X.ds_conv(x, 0, 0.6, 0.5, method="quad") / math.pi - X.voigt(x, 0.5, 0.6))) / np.max(X.voigt(x, 0.5, 0.6))
lim["DS(a=0,wG=0)/pi->L"] = np.max(np.abs(X.ds_raw(x, 0, 0.6) / math.pi - X.lorentz(x, 0.6))) / np.max(X.lorentz(x, 0.6))
lim["GLP(m=0)->G"] = np.max(np.abs(X.glp(x, 1.1, 0.0) - X.gauss(x, 1.1))) / np.max(X.gauss(x, 1.1))
lim["GLP(m=1)->L"] = np.max(np.abs(X.glp(x, 1.1, 1.0) - X.lorentz(x, 1.1))) / np.max(X.lorentz(x, 1.1))
lim["GLT(T=0)->GLS"] = np.max(np.abs(X.glt(x, 1.1, 0.3, 0.0, 1.0) - X.gls(x, 1.1, 0.3)))
rec("limits (max rel dev)", {k: float(v) for k, v in lim.items()}, max(lim.values()) < 1e-6)

# convolution accuracy: fine grid vs adaptive quad; closed-form GLT vs quad
xx = np.linspace(-4, 8, 25)
dev = {}
for a, wL, wG in [(0.1, 0.3, 0.6), (0.3, 0.2, 0.5), (0.45, 0.4, 0.4), (0.05, 0.1, 0.8)]:
    q = X.ds_conv(xx, a, wL, wG, method="quad"); f = X.ds_conv(xx, a, wL, wG)
    dev[f"DS{a},{wL},{wG}"] = float(np.max(np.abs(f - q)) / np.max(q))
for F, a, b, wG in [(1.2, 2.5, 1.2, 0.4), (0.6, 0.6, 5, 0.3), (2, 20, 0.8, 1.0)]:
    q = X.la(xx, F, a, b, wG, method="quad"); f = X.la(xx, F, a, b, wG)
    dev[f"LA{F},{a},{b},{wG}"] = float(np.max(np.abs(f - q)) / np.max(q))
for F, m, tau in [(1, 0.3, 0.5), (0.8, 1.0, 0.2), (1.5, 0.5, 3.0), (1.0, 0.2, 0.05)]:
    q = X.glt_tail_quad(xx, F, m, tau)
    f = (1 - m) * X.gauss_exp_conv(xx, F, tau) + m * X.lorentz_exp_conv(xx, F, tau)
    dev[f"GLTtail{F},{m},{tau}"] = float(np.max(np.abs(f - q)) / np.max(q))
rec("convolutions: fine-grid / closed form vs adaptive quad (rel to max)", dev, max(dev.values()) < 1e-7)
ri = {}
for a, wL, wG, reg in [(0.1, 0.3, 0.6, (-3, 6)), (0.3, 0.2, 0.5, (-2, 9)), (0.45, 0.4, 0.0, (-3, 5))]:
    an = X.ds_region_integral(reg[0], reg[1], a, wL, wG, method="quad")
    br, _ = integrate.quad(lambda t: X.ds_conv(np.array([t]), a, wL, wG, method="quad")[0], reg[0], reg[1],
                           points=[0.0], limit=400, epsrel=1e-11)
    ri[f"{a},{wL},{wG},{reg}"] = abs(an / br - 1)
rec("DS region integral: analytic-inner vs nested quad (rel)", ri, max(ri.values()) < 1e-9)

# ---------------------------------------------------------------- backgrounds
bg = json.load(open(os.path.join(FIX, "backgrounds.json")))["cases"]
it = {c["id"]: c["info"]["iterations"] for c in bg if c["type"] in ("shirley", "smart")}
conv = all(c["info"]["converged"] for c in bg if c["type"] in ("shirley", "smart"))
rec("Shirley/smart converged (iterations)", it, conv)
# fixed-point residual of converged Shirley
fp = {}
for c in bg:
    if c["type"] != "shirley":
        continue
    E = np.array(c["E"]); y = np.array(c["y"]); B = np.array(c["B"]); ep = c["endpoints"]
    Q = X.cumtrapz(y - B, E)
    Bn = ep["y_lo"] + (ep["y_hi"] - ep["y_lo"]) * Q / Q[-1]
    fp[c["id"]] = float(np.max(np.abs(Bn - B)) / (abs(ep["y_hi"] - ep["y_lo"]) + 1e-12))
rec("Shirley fixed-point residual max|T(B)-B|/|dy|", max(fp.values()), max(fp.values()) < 1e-6)
syn = {c["id"]: c["max_abs_dev_from_expected"] for c in bg if "max_abs_dev_from_expected" in c}
rec("synthetic backgrounds vs analytic expectation (abs)", syn, max(syn.values()) < 1e-5)

# third-party comparisons (n_end = 1, the only endpoint rule they implement)
tp = {}
try:
    from lmfitxps.backgrounds import shirley_calculate, tougaard_calculate
    have_lmfitxps = True
except Exception:
    have_lmfitxps = False
try:
    from pynxtools_xps.models.backgrounds import Shirley as PnxShirley
    have_pnx = True
except Exception:
    have_pnx = False
regions = sorted({(c["file"], c["E"][0], c["E"][-1]) for c in bg if c["file"]})
tp_info = {}
for c in bg:
    if not c["file"] or c["type"] != "linear":
        continue
    E = np.array(c["E"]); y = np.array(c["y"])
    Ed = E[::-1]; yd = y[::-1]   # descending BE like the instrument file
    tag = c["id"].rsplit("_", 1)[0]
    Bme = X.bg_shirley(E, y, 1, tol=1e-13, maxit=5000)
    scale = abs(y[-1] - y[0]) + 1e-12
    if have_lmfitxps:
        Bl = shirley_calculate(Ed, yd, tol=1e-20, maxit=5000)[::-1]
        tp[f"{tag} shirley vs lmfitxps"] = float(np.max(np.abs(Bme - Bl)) / scale)
        for kind, C, D, tcd, td in (("U2", 1643.0, 0.0, 1.0, 0.0), ("U3", 551.0, 436.0, -1.0, 436.0)):
            Bt, info = X.bg_tougaard(E, y, 1, kind, C, D, return_info=True)
            if info["s"] == 0.0 and y[-1] < y[0]:
                # SPEC v1.1 clamps s < 0 to 0 (unphysical negative loss); lmfitxps does not -> not comparable
                tp_info[f"{tag} tougaard {kind}"] = "s<0 clamped (SPEC 2 v1.1), not compared"
                continue
            with contextlib.redirect_stdout(io.StringIO()):
                Bl, tb = tougaard_calculate(Ed, yd, tb=2866, tc=C, tcd=tcd, td=td, maxit=2000)
            Bl = np.asarray(Bl)[::-1]
            tp[f"{tag} tougaard {kind} vs lmfitxps"] = float(np.max(np.abs(Bt - Bl)) / scale)
    if have_pnx:
        with np.errstate(all="ignore"):
            try:
                Bp = PnxShirley().calc_background(Ed, yd, tol=1e-12, maxit=5000)[::-1]
                # pynxtools-xps integrates y - B without removing the end-point offset y_end (its B starts at 0),
                # so its k and B differ from Proctor-Sherwood; reported for information only.
                tp_info[f"{tag} shirley vs pynxtools-xps (different algorithm)"] = float(np.max(np.abs(Bme - Bp)) / scale)
            except Exception as e:
                tp[f"{tag} shirley vs pynxtools-xps"] = repr(e)
num = [v for v in tp.values() if isinstance(v, float)]
rec("third-party background agreement (max|dB|/|y_end-y_start|, n_end=1)", tp,
    (max(num) < 1e-5) if num else None)
rec("third-party background notes", tp_info)

# ---------------------------------------------------------------- fits
fits = json.load(open(os.path.join(FIX, "fits.json")))["cases"]
opt = {}
for c in fits:
    E = np.array(c["E"]); y = np.array(c["y"]); r = c["result"]
    m = X.Model(c["model"], E, y, c["c"])
    v = np.array(r["v"]); th = np.array(r["theta"])
    sw = 1 / np.sqrt(v)
    chi = lambda t: float(np.sum(((y - m.mu(t)) * sw) ** 2))
    c0 = chi(th)
    # independent polish with a different algorithm (bounded Powell) from the reported optimum
    bnds = list(zip(np.where(np.isfinite(m.lo), m.lo, None), np.where(np.isfinite(m.hi), m.hi, None)))
    res = optimize.minimize(chi, th, method="Powell", bounds=bnds,
                            options=dict(xtol=1e-10, ftol=1e-14, maxfev=4000))
    # gradient of chi2 at optimum for interior params (scaled by sigma)
    J = X.jacobian_mu(m, th) * sw[:, None]
    g = -2 * J.T @ ((y - m.mu(th)) * sw)
    sig = np.array([s if s is not None else np.nan for s in r["sigma_vec"]])
    gs = np.nanmax(np.abs(g * sig)) if np.any(np.isfinite(sig)) else None
    extra = {}
    if any(cp["shape"] == "DS" for cp in c["model"]["components"]):
        mq = X.Model(c["model"], E, y, c["c"], method="quad")
        extra["DS mu fine-vs-quad rel"] = float(np.max(np.abs(mq.mu(th) - m.mu(th))) / np.max(m.mu(th)))
    opt[c["id"]] = dict(chi2=c0, powell_improvement=c0 - res.fun, max_abs_grad_x_sigma=float(gs), **extra)
worst = max(o["powell_improvement"] for o in opt.values())
rec("fit optimality (Powell polish from reported optimum; |dchi2|, |grad|*sigma)", opt, worst < 1e-6)

# heavier global search: 100 LHS starts for each fixture, pass 1 only
glob_ = {}
for c in fits:
    E = np.array(c["E"]); y = np.array(c["y"])
    if any(cp["shape"] in ("DS",) for cp in c["model"]["components"]):
        nst = 24
    elif any(cp["shape"] in ("LA",) for cp in c["model"]["components"]):
        nst = 40
    else:
        nst = 100
    m = X.Model(c["model"], E, y, c["c"])
    with contextlib.redirect_stdout(io.StringIO()):
        r2 = X.fit(m, n_starts=nst, seed=999)
    glob_[c["id"]] = dict(n_starts=nst, chi2_fixture=c["result"]["chi2"], chi2_recheck=r2["chi2"],
                          diff=r2["chi2"] - c["result"]["chi2"])
rec("global-minimum recheck (different seed, more starts): chi2_recheck - chi2_fixture", glob_,
    min(g["diff"] for g in glob_.values()) > -1e-6)

json.dump(R, open(os.path.join(FIX, "validation.json"), "w"), indent=1, default=float)
print("wrote fixtures/validation.json")
