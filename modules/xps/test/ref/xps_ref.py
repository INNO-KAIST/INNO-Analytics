"""
xps_ref.py -- independent Python reference implementation of SPEC.md (XPS module v1).

Written only from /home/claude/ia/modules/xps/SPEC.md (no look at the JS engine).
Numerics are deliberately higher-accuracy than the spec's discretisation (sec. 1.8):
  * Voigt via scipy.special.wofz (Faddeeva)
  * GLS (x) exp tail (GLT) in closed form (ex-Gaussian via erfcx, Lorentz (x) exp via complex E1)
  * DS (x) G and LA (x) G: either adaptive quad (fixtures, `*_quad`) or a very fine
    trapezoid convolution that is exponentially accurate for smooth integrands (fits).
  * DS region normalisation: inner integral done analytically with the DS antiderivative.
Backgrounds (sec. 2) are implemented *literally* on the data grid (trapezoid), because the
spec defines them as discrete operations.

Conventions: E is binding energy in ascending order; x = E - E0 (x>0 = higher BE).
"""
from __future__ import annotations

import math
import numpy as np
from scipy import special, integrate, optimize
from scipy.stats import qmc

LN2 = math.log(2.0)
SQ2 = math.sqrt(2.0)
FWHM2SIG = 1.0 / (2.0 * math.sqrt(2.0 * LN2))   # sigma = FWHM * FWHM2SIG

# ---------------------------------------------------------------------------
# 1.1 basic functions
# ---------------------------------------------------------------------------

def gauss(x, F):
    x = np.asarray(x, float)
    return (2.0 / F) * math.sqrt(LN2 / math.pi) * np.exp(-4.0 * LN2 * x * x / (F * F))


def lorentz(x, F):
    x = np.asarray(x, float)
    return (2.0 / (math.pi * F)) / (1.0 + 4.0 * x * x / (F * F))


# ---------------------------------------------------------------------------
# 1.2 GLS, 1.3 GLP, 1.4 Voigt
# ---------------------------------------------------------------------------

def gls(x, F, m):
    return (1.0 - m) * gauss(x, F) + m * lorentz(x, F)


def glp_J(m):
    a = 4.0 * LN2 * (1.0 - m)
    b = 4.0 * m
    if b == 0.0:
        return math.sqrt(math.pi / a)
    return (math.pi / math.sqrt(b)) * float(special.erfcx(math.sqrt(a / b)))


def glp(x, F, m):
    x = np.asarray(x, float)
    t2 = x * x / (F * F)
    g = np.exp(-4.0 * LN2 * (1.0 - m) * t2) / (1.0 + 4.0 * m * t2)
    return g / (F * glp_J(m))


def voigt(x, wG, wL):
    x = np.asarray(x, float)
    if wG <= 0.0:
        return lorentz(x, wL)
    if wL <= 0.0:
        return gauss(x, wG)
    sig = wG * FWHM2SIG
    gam = wL / 2.0
    z = (x + 1j * gam) / (sig * SQ2)
    return np.real(special.wofz(z)) / (sig * math.sqrt(2.0 * math.pi))


def voigt_fwhm_OL(wG, wL):
    """Olivero-Longbothum approximation (display value in spec 1.4)."""
    return 0.5346 * wL + math.sqrt(0.2166 * wL * wL + wG * wG)


# ---------------------------------------------------------------------------
# 1.5 Doniach-Sunjic
# ---------------------------------------------------------------------------

def ds_raw(x, alpha, wL):
    """Un-normalised DS kernel, tail to x>0 (higher BE)."""
    x = np.asarray(x, float)
    g = wL / 2.0
    return np.cos(math.pi * alpha / 2.0 + (1.0 - alpha) * np.arctan(-x / g)) / \
        (x * x + g * g) ** ((1.0 - alpha) / 2.0)


def ds_antideriv(x, alpha, wL):
    """Antiderivative of ds_raw.  ds = Re[e^{-i pi a/2} (g - i x)^(a-1)],
    so  int ds dx = Re[ i e^{-i pi a/2} (g - i x)^a / a ]  (a>0),  atan(x/g) (a=0)."""
    x = np.asarray(x, float)
    g = wL / 2.0
    if alpha == 0.0:
        return np.arctan(x / g)
    z = g - 1j * x
    return np.real(1j * np.exp(-1j * math.pi * alpha / 2.0) * z ** alpha / alpha)


def _gauss_grid(sig, h, nsig=10.0):
    K = int(math.ceil(nsig * sig / h))
    t = np.arange(-K, K + 1) * h
    w = np.exp(-0.5 * (t / sig) ** 2) / (sig * math.sqrt(2 * math.pi)) * h
    return t, w


def _conv_fine(func, x, sig, h):
    """p(x) = int func(x - t) G_sigma(t) dt by a fine trapezoid on +-10 sigma (exponentially
    accurate for analytic integrands; func must be vectorised)."""
    x = np.asarray(x, float)
    t, w = _gauss_grid(sig, h)
    out = np.empty_like(x)
    # chunk to bound memory
    step = max(1, int(4_000_000 // len(t)))
    for i in range(0, len(x), step):
        xi = x[i:i + step]
        out[i:i + step] = func(xi[:, None] - t[None, :]) @ w
    return out


def _conv_quad(func, x, sig, brk=0.0):
    """Adaptive-quadrature Gaussian convolution, used for fixtures/validation."""
    out = []
    for xi in np.atleast_1d(np.asarray(x, float)):
        lo, hi = -12 * sig, 12 * sig
        pts = [p for p in (xi - brk,) if lo < p < hi]
        val, err = integrate.quad(lambda t: func(np.array([xi - t]))[0] *
                                  math.exp(-0.5 * (t / sig) ** 2), lo, hi,
                                  points=pts or None, limit=400, epsabs=0.0, epsrel=1e-12)
        out.append(val / (sig * math.sqrt(2 * math.pi)))
    return np.array(out)


def ds_conv(x, alpha, wL, wG, method="fine", h=None):
    """p = ds (x) G(wG).  (unnormalised)"""
    if wG <= 0.0:
        return ds_raw(x, alpha, wL)
    sig = wG * FWHM2SIG
    f = lambda u: ds_raw(u, alpha, wL)
    if method == "quad":
        return _conv_quad(f, x, sig)
    if h is None:
        h = min(sig, wL / 2.0) / 40.0
    return _conv_fine(f, x, sig, h)


def ds_region_integral(a, b, alpha, wL, wG, method="fine", h=None):
    """int_a^b p(x) dx  with x relative to E0 (inner integral analytic)."""
    if wG <= 0.0:
        return float(ds_antideriv(b, alpha, wL) - ds_antideriv(a, alpha, wL))
    sig = wG * FWHM2SIG
    if method == "quad":
        f = lambda t: (ds_antideriv(b - t, alpha, wL) - ds_antideriv(a - t, alpha, wL)) * \
            math.exp(-0.5 * (t / sig) ** 2) / (sig * math.sqrt(2 * math.pi))
        pts = [p for p in (a, b) if -12 * sig < p < 12 * sig]
        v, _ = integrate.quad(f, -12 * sig, 12 * sig, points=pts or None, limit=400,
                              epsabs=0.0, epsrel=1e-13)
        return v
    if h is None:
        h = min(sig, wL / 2.0) / 40.0
    t, w = _gauss_grid(sig, h)
    return float(((ds_antideriv(b - t, alpha, wL) - ds_antideriv(a - t, alpha, wL))) @ w)


def ds(x, alpha, wL, wG, region, method="fine", alpha0_norm="R"):
    """Normalised DS (spec 1.5).  region = (xmin, xmax) RELATIVE TO E0.
    alpha>0: normalised over region.  alpha==0: literal reading -> normalised over R
    (= Voigt).  alpha0_norm='region' gives the other reading."""
    p = ds_conv(x, alpha, wL, wG, method=method)
    if alpha == 0.0 and alpha0_norm == "R":
        return p / math.pi
    return p / ds_region_integral(region[0], region[1], alpha, wL, wG, method=method)


# ---------------------------------------------------------------------------
# 1.6 LA
# ---------------------------------------------------------------------------

def la_N(F, alpha, beta):
    sp = math.sqrt(math.pi)
    return (F / 4.0) * sp * (math.exp(special.gammaln(alpha - 0.5) - special.gammaln(alpha)) +
                             math.exp(special.gammaln(beta - 0.5) - special.gammaln(beta)))


def la_raw(x, F, alpha, beta):
    x = np.asarray(x, float)
    u = 1.0 + 4.0 * x * x / (F * F)
    return np.where(x >= 0.0, u ** (-alpha), u ** (-beta))   # SPEC 1.6 v1.1: alpha = high-BE side


def la(x, F, alpha, beta, wG, method="fine", h=None):
    N = la_N(F, alpha, beta)
    f = lambda u: la_raw(u, F, alpha, beta) / N
    if wG <= 0.0:
        return f(np.asarray(x, float))
    sig = wG * FWHM2SIG
    if method == "quad":
        return _conv_quad(f, x, sig)
    if h is None:
        h = min(sig, F / 2.0) / 60.0
    return _conv_fine(f, x, sig, h)


# ---------------------------------------------------------------------------
# 1.7 GLT  (closed form)
# ---------------------------------------------------------------------------

def gauss_exp_conv(x, F, tau):
    """(G(.;F) (x) e)(x), e(t)=exp(-t/tau)/tau for t>=0  (ex-Gaussian)."""
    x = np.asarray(x, float)
    s = F * FWHM2SIG
    z = (s / tau - x / s) / SQ2
    out = np.empty_like(x)
    pos = z >= 0
    out[pos] = np.exp(-0.5 * (x[pos] / s) ** 2) * special.erfcx(z[pos])
    # z<0: erfcx(z) = 2 exp(z^2) - erfcx(-z);  z^2 - x^2/2s^2 = s^2/2tau^2 - x/tau
    xn = x[~pos]
    out[~pos] = 2.0 * np.exp(0.5 * (s / tau) ** 2 - xn / tau) - \
        np.exp(-0.5 * (xn / s) ** 2) * special.erfcx(-z[~pos])
    return (0.5 / tau) * out


def lorentz_exp_conv(x, F, tau):
    """(L(.;F) (x) e)(x) = (1/(pi tau)) Im[ e^{zz} E1(zz) ],  zz = -(x+i g)/tau."""
    x = np.asarray(x, float)
    g = F / 2.0
    zz = -(x + 1j * g) / tau
    # e^{z}E1(z) for large |z| loses nothing with scipy complex exp1 at these sizes;
    # use asymptotic series only where |z|>~700 to avoid overflow.
    out = np.empty(x.shape)
    big = np.abs(zz) > 600
    zs = zz[~big]
    out[~big] = np.imag(np.exp(zs) * special.exp1(zs)) / (math.pi * tau)
    if np.any(big):
        zb = zz[big]
        s = np.zeros_like(zb)
        term = np.ones_like(zb)
        for k in range(12):
            s += term
            term = term * (-(k + 1)) / zb
        out[big] = np.imag(s / zb) / (math.pi * tau)
    return out


def glt(x, F, m, T, tau):
    base = gls(x, F, m)
    if T == 0.0:
        return base
    tail = (1.0 - m) * gauss_exp_conv(x, F, tau) + m * lorentz_exp_conv(x, F, tau)
    return (1.0 - T) * base + T * tail


def glt_tail_quad(x, F, m, tau):
    """quad check of (GLS (x) e)."""
    out = []
    for xi in np.atleast_1d(np.asarray(x, float)):
        f = lambda t: gls(np.array([xi - t]), F, m)[0] * math.exp(-t / tau) / tau
        pts = [xi] if xi > 0 else None
        v1, _ = integrate.quad(f, 0.0, max(xi, 0.0) + 60 * tau if pts else 60 * tau,
                               points=pts, limit=500, epsabs=0, epsrel=1e-12)
        out.append(v1)
    return np.array(out)


# ---------------------------------------------------------------------------
# generic evaluation
# ---------------------------------------------------------------------------
SHAPE_PARAMS = {
    "GLS": ["F", "m"], "GLP": ["F", "m"], "V": ["wG", "wL"],
    "DS": ["alpha", "wL", "wG"], "LA": ["F", "alpha", "beta", "wG"],
    "GLT": ["F", "m", "T", "tau"],
}
WIDTH_PARAMS = {"F", "wL"}   # SPEC 1.9 v1.1: q scales the lifetime width only (wG common)


def shape_eval(shape, x, p, region=None, method="fine"):
    """Unit-area (or region-normalised for DS) lineshape at x = E - E0."""
    if shape == "GLS":
        return gls(x, p["F"], p["m"])
    if shape == "GLP":
        return glp(x, p["F"], p["m"])
    if shape == "V":
        return voigt(x, p["wG"], p["wL"])
    if shape == "DS":
        return ds(x, p["alpha"], p["wL"], p["wG"], region, method=method)
    if shape == "LA":
        return la(x, p["F"], p["alpha"], p["beta"], p["wG"], method=method)
    if shape == "GLT":
        return glt(x, p["F"], p["m"], p["T"], p["tau"])
    raise ValueError(shape)


def doublet_partner_params(p, q):
    pp = dict(p)
    for k in WIDTH_PARAMS:
        if k in pp:
            pp[k] = pp[k] * q
    return pp


def numeric_fwhm(fun, x0_guess, scale):
    """Numerical FWHM of a single-maximum function fun(x) (scalar->scalar)."""
    xs = np.linspace(x0_guess - 6 * scale, x0_guess + 6 * scale, 6001)
    ys = np.array([fun(v) for v in xs])
    i = int(np.argmax(ys))
    res = optimize.minimize_scalar(lambda v: -fun(v), bracket=(xs[max(i - 1, 0)], xs[i], xs[min(i + 1, len(xs) - 1)]),
                                   tol=1e-12)
    xm, fm = res.x, -res.fun
    half = fm / 2.0
    g = lambda v: fun(v) - half
    # walk outwards
    def find(direction):
        st = scale / 50.0
        a = xm
        b = xm + direction * st
        while g(b) > 0:
            a = b
            st *= 1.5
            b = b + direction * st
            if abs(b - xm) > 1e4 * scale:
                return None
        return optimize.brentq(g, min(a, b), max(a, b), xtol=1e-14, rtol=1e-14, maxiter=500)
    l = find(-1)
    r = find(+1)
    return (r - l), xm, fm


# ---------------------------------------------------------------------------
# 2. Backgrounds
# ---------------------------------------------------------------------------

def cumtrapz(y, E):
    y = np.asarray(y, float)
    E = np.asarray(E, float)
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(E))
    return out


def endpoints(E, y, n=5):
    E = np.asarray(E, float)
    y = np.asarray(y, float)
    return dict(y_lo=float(np.mean(y[:n])), y_hi=float(np.mean(y[-n:])),
                E_lo=float(np.mean(E[:n])), E_hi=float(np.mean(E[-n:])))


def round_half_up(v):
    return int(math.floor(v + 0.5))


def smooth_h(E):
    dE = float(np.median(np.abs(np.diff(E))))
    return max(2, round_half_up(0.25 / dE + 1e-9))


def smooth_quadratic(E, y, h=None):
    """Local quadratic LSQ smoothing, centred window 2h+1, asymmetric (truncated) at edges."""
    E = np.asarray(E, float)
    y = np.asarray(y, float)
    if h is None:
        h = smooth_h(E)
    N = len(y)
    ys = np.empty(N)
    for i in range(N):
        a, b = max(0, i - h), min(N - 1, i + h)
        xx = E[a:b + 1] - E[i]
        V = np.vander(xx, 3, increasing=True)
        coef, *_ = np.linalg.lstsq(V, y[a:b + 1], rcond=None)
        ys[i] = coef[0]
    return ys


def bg_linear(E, y, n=5):
    ep = endpoints(E, y, n)
    E = np.asarray(E, float)
    return ep["y_lo"] + (ep["y_hi"] - ep["y_lo"]) * (E - ep["E_lo"]) / (ep["E_hi"] - ep["E_lo"])


def bg_shirley(E, y, n=5, smart=False, ys=None, tol=1e-7, maxit=200, return_info=False):
    E = np.asarray(E, float)
    y = np.asarray(y, float)
    ep = endpoints(E, y, n)
    ylo, yhi = ep["y_lo"], ep["y_hi"]
    if yhi - ylo == 0.0:
        B = np.full_like(y, ylo)
        return (B, dict(iterations=0, converged=True)) if return_info else B
    if smart and ys is None:
        ys = smooth_quadratic(E, y)
    B = np.full_like(y, ylo)
    thr = tol * (abs(yhi - ylo) + 1e-12)
    it = 0
    conv = False
    for it in range(1, maxit + 1):
        Q = cumtrapz(y - B, E)
        Bn = ylo + (yhi - ylo) * Q / Q[-1]
        if smart:
            Bn = np.minimum(Bn, ys)
        d = np.max(np.abs(Bn - B))
        B = Bn
        if d < thr:
            conv = True
            break
    return (B, dict(iterations=it, converged=conv)) if return_info else B


def tougaard_K(T, kind="U2", C=1643.0, D=0.0):
    T = np.asarray(T, float)
    if kind == "U2":
        return T / (C + T * T) ** 2
    if kind == "U3":
        return T / ((C - T * T) ** 2 + D * T * T)
    raise ValueError(kind)


def bg_tougaard(E, y, n=5, kind="U2", C=1643.0, D=0.0, return_info=False):
    E = np.asarray(E, float)
    y = np.asarray(y, float)
    ep = endpoints(E, y, n)
    ylo, yhi = ep["y_lo"], ep["y_hi"]
    N = len(E)
    I = np.zeros(N)
    yy = y - ylo
    for k in range(1, N):
        Ek = E[k]
        integrand = tougaard_K(Ek - E[:k + 1], kind, C, D) * yy[:k + 1]
        I[k] = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(E[:k + 1]))
    I_hi = float(np.interp(ep["E_hi"], E, I))
    s = (yhi - ylo) / I_hi if I_hi != 0 else 0.0
    s = max(0.0, s)   # SPEC 2 v1.1: no negative loss cross-section
    B = ylo + s * I
    return (B, dict(s=s, I_hi=I_hi)) if return_info else B


def background(E, y, spec):
    t = spec["type"]
    n = spec.get("n_end", 5)
    if t == "linear":
        return bg_linear(E, y, n)
    if t == "shirley":
        return bg_shirley(E, y, n)
    if t == "smart":
        return bg_shirley(E, y, n, smart=True)
    if t == "tougaard":
        pr = spec.get("params", {})
        return bg_tougaard(E, y, n, pr.get("kind", "U2"), pr.get("C", 1643.0), pr.get("D", 0.0))
    if t == "none":
        return np.zeros(len(E))
    raise ValueError(t)


# ---------------------------------------------------------------------------
# 3. Fitting
# ---------------------------------------------------------------------------
class Model:
    """Model definition (see README for the JSON schema).

    model = {"background": {...}, "components": [ {name, shape, E0:P, A:P, <shape params>:P,
             "doublet": None | {"delta":P, "r":P, "q":P}} ]}
    P = {"mode":"free","value":v,"lo":lo,"hi":hi} | {"mode":"fixed","value":v}
        | {"mode":"link","ref":compname,"k":1.0,"b":0.0}
    Background "active": {"type":"active","c0":P,"k":P,"s":P}
    """

    def __init__(self, model, E, y, c=1.0, method="fine"):
        self.m = model
        self.E = np.asarray(E, float)
        self.y = np.asarray(y, float)
        self.c = float(c)
        self.method = method
        self.Emin, self.Emax = float(self.E[0]), float(self.E[-1])
        self.comps = model["components"]
        self.cnames = [cp["name"] for cp in self.comps]
        bg = model["background"]
        self.active = bg["type"] == "active"
        self.ys = smooth_quadratic(self.E, self.y)
        self.Bstatic = None if self.active else background(self.E, self.y, bg)
        # enumerate parameters
        self.keys = []      # (owner, pname)
        self.pdef = {}
        for cp in self.comps:
            names = ["E0", "A"] + SHAPE_PARAMS[cp["shape"]]
            for pn in names:
                self.keys.append((cp["name"], pn))
                self.pdef[(cp["name"], pn)] = cp[pn]
            if cp.get("doublet"):
                for pn in ("delta", "r", "q"):
                    self.keys.append((cp["name"], pn))
                    self.pdef[(cp["name"], pn)] = cp["doublet"][pn]
        if self.active:
            for pn in ("c0", "k", "s"):
                d = bg.get(pn, {"mode": "fixed", "value": 0.0})
                self.keys.append(("bg", pn))
                self.pdef[("bg", pn)] = d
        self.free = [k for k in self.keys if self.pdef[k]["mode"] == "free"]
        self.lo = []
        self.hi = []
        for k in self.free:
            d = self.pdef[k]
            lo = d.get("lo", -np.inf)
            hi = d.get("hi", np.inf)
            if k[1] == "A":
                lo = max(0.0 if lo is None else lo, 0.0)
            if k[1] == "k" and k[0] == "bg":
                lo = max(0.0 if lo is None else lo, 0.0)
            self.lo.append(-np.inf if lo is None else lo)
            self.hi.append(np.inf if hi is None else hi)
        self.lo = np.array(self.lo, float)
        self.hi = np.array(self.hi, float)

    def start(self):
        return np.array([self.pdef[k]["value"] for k in self.free], float)

    def full(self, theta):
        vals = {}
        for k, v in zip(self.free, theta):
            vals[k] = float(v)
        for k in self.keys:
            if self.pdef[k]["mode"] == "fixed":
                vals[k] = float(self.pdef[k]["value"])
        # resolve links (allow chains)
        for _ in range(10):
            for k in self.keys:
                d = self.pdef[k]
                if d["mode"] == "link":
                    ref = (d["ref"], k[1] if "param" not in d else d["param"])
                    if ref in vals:
                        vals[k] = d.get("k", 1.0) * vals[ref] + d.get("b", 0.0)
        return vals

    def components(self, theta):
        vals = self.full(theta)
        region_abs = (self.Emin, self.Emax)
        out = {}
        for cp in self.comps:
            nm = cp["name"]
            sp = {pn: vals[(nm, pn)] for pn in SHAPE_PARAMS[cp["shape"]]}
            E0 = vals[(nm, "E0")]
            A = vals[(nm, "A")]
            if cp.get("doublet"):
                d, r, q = vals[(nm, "delta")], vals[(nm, "r")], vals[(nm, "q")]
                A1 = A / (1 + r)
                A2 = A * r / (1 + r)
                f1 = shape_eval(cp["shape"], self.E - E0, sp,
                                (region_abs[0] - E0, region_abs[1] - E0), self.method)
                sp2 = doublet_partner_params(sp, q)
                E02 = E0 + d
                f2 = shape_eval(cp["shape"], self.E - E02, sp2,
                                (region_abs[0] - E02, region_abs[1] - E02), self.method)
                out[nm] = A1 * f1
                out[nm + "#2"] = A2 * f2
            else:
                out[nm] = A * shape_eval(cp["shape"], self.E - E0, sp,
                                         (region_abs[0] - E0, region_abs[1] - E0), self.method)
        return out, vals

    def mu(self, theta, return_parts=False):
        parts, vals = self.components(theta)
        P = np.sum(list(parts.values()), axis=0)
        if self.active:
            B = vals[("bg", "c0")] + vals[("bg", "k")] * cumtrapz(P, self.E) + \
                vals[("bg", "s")] * (self.E - self.Emin)
        else:
            B = self.Bstatic
        m = B + P
        if return_parts:
            return m, B, parts, vals
        return m


def variance(yhat, c):
    vmin = max(1.0 / c, 1e-6 * float(np.max(yhat)))
    return np.maximum(yhat, vmin) / c


def nnls_areas(model, theta):
    """Initial areas by NNLS with nonlinear params fixed (spec 3.3)."""
    idxA = [i for i, k in enumerate(model.free) if k[1] == "A"]
    if not idxA:
        return theta
    th = theta.copy()
    cols = []
    for i in idxA:
        t1 = th.copy()
        for j in idxA:
            t1[j] = 0.0
        t1[i] = 1.0
        base = th.copy()
        for j in idxA:
            base[j] = 0.0
        cols.append(model.mu(t1) - model.mu(base))
    base = th.copy()
    for j in idxA:
        base[j] = 0.0
    r0 = model.y - model.mu(base)
    w = 1.0 / np.sqrt(variance(model.ys, model.c))
    M = np.array(cols).T * w[:, None]
    a, _ = optimize.nnls(M, r0 * w)
    for i, v in zip(idxA, a):
        th[i] = min(max(v, model.lo[i]), model.hi[i]) if np.isfinite(model.hi[i]) else max(v, model.lo[i])
        if th[i] <= 0:
            th[i] = 1e-3 * max(1.0, float(np.max(np.abs(r0))))
    return th


def _solve(model, v, theta0):
    sw = 1.0 / np.sqrt(v)
    fun = lambda th: (model.y - model.mu(th)) * sw
    with np.errstate(invalid="ignore"):
        lo = np.where(np.isfinite(model.lo), model.lo + 1e-9 * (np.abs(model.lo) + 1), -np.inf)
        hi = np.where(np.isfinite(model.hi), model.hi - 1e-9 * (np.abs(model.hi) + 1), np.inf)
    th0 = np.clip(theta0, lo, hi)
    res = optimize.least_squares(fun, th0, bounds=(model.lo, model.hi), method="trf",
                                 x_scale="jac", ftol=1e-15, xtol=1e-15, gtol=1e-15,
                                 max_nfev=20000, diff_step=1e-7)
    return res.x, float(np.sum(res.fun ** 2))


def jacobian_mu(model, theta):
    """Central differences with one Richardson step (O(h^4))."""
    J = np.zeros((len(model.E), len(theta)))
    for j in range(len(theta)):
        h = 1e-4 * max(abs(theta[j]), 1e-2)
        if model.free[j][1] in ("E0", "delta"):
            h = 1e-4
        def d(hh):
            tp = theta.copy(); tm = theta.copy()
            tp[j] += hh; tm[j] -= hh
            return (model.mu(tp) - model.mu(tm)) / (2 * hh)
        J[:, j] = (4 * d(h / 2) - d(h)) / 3.0
    return J


def runs_z(r):
    s = np.sign(r)
    s = s[s != 0]
    n1 = np.sum(s > 0); n2 = np.sum(s < 0)
    runs = 1 + np.sum(s[1:] != s[:-1])
    n = n1 + n2
    mu = 2 * n1 * n2 / n + 1
    var = 2 * n1 * n2 * (2 * n1 * n2 - n) / (n * n * (n - 1))
    return float((runs - mu) / math.sqrt(var)), int(runs)


def fit(model, n_starts=24, seed=12345, extra_starts=()):
    """Spec 3: pass 1 with v from smoothed data, pass 2 with v from mu_1 (IRLS 2 passes)."""
    v1 = variance(model.ys, model.c)
    starts = [nnls_areas(model, model.start())]
    for s in extra_starts:
        starts.append(nnls_areas(model, np.asarray(s, float)))
    nl = [i for i, k in enumerate(model.free) if k[1] != "A"]
    if n_starts > 1 and nl:
        sampler = qmc.LatinHypercube(d=len(nl), seed=seed)
        U = sampler.random(n_starts)
        for u in U:
            th = model.start().copy()
            for jj, i in enumerate(nl):
                lo, hi = model.lo[i], model.hi[i]
                if not np.isfinite(lo) or not np.isfinite(hi):
                    c0 = th[i]
                    sc = max(abs(c0) * 0.5, 1.0)
                    lo = c0 - sc if not np.isfinite(lo) else lo
                    hi = c0 + sc if not np.isfinite(hi) else hi
                th[i] = lo + u[jj] * (hi - lo)
            starts.append(nnls_areas(model, th))
    sols = []
    for th0 in starts:
        try:
            th, c2 = _solve(model, v1, th0)
            sols.append((c2, th))
        except Exception as e:  # pragma: no cover
            print("   start failed:", repr(e))
    sols.sort(key=lambda t: t[0])
    chi1, th1 = sols[0]
    uniq = float(np.mean([abs(c - chi1) < 1.0 for c, _ in sols]))
    mu1 = model.mu(th1)
    v2 = variance(mu1, model.c)
    # pass 2: start from pass-1 best and from the 3 best distinct pass-1 solutions
    cand = []
    for c, th in sols[:4]:
        cand.append(_solve(model, v2, th))
    cand.sort(key=lambda t: t[1])
    th2, chi2 = cand[0]
    return summarize(model, th2, v2, extra=dict(chi2_pass1=chi1, theta_pass1=th1.tolist(),
                                                 uniqueness=uniq, n_starts=len(sols)))


def summarize(model, theta, v, extra=None):
    mu, B, parts, vals = model.mu(theta, return_parts=True)
    r = (model.y - mu) / np.sqrt(v)
    chi2 = float(np.sum(r * r))
    N = len(model.y)
    p = len(theta)
    phi = chi2 / (N - p)
    # bound detection
    rng = np.where(np.isfinite(model.hi - model.lo), model.hi - model.lo, np.abs(theta) + 1.0)
    at_lo = np.abs(theta - model.lo) <= 1e-6 * rng
    at_hi = np.abs(theta - model.hi) <= 1e-6 * rng
    atb = at_lo | at_hi
    J = jacobian_mu(model, theta) / np.sqrt(v)[:, None]
    keep = ~atb
    sig = np.full(p, np.nan)
    singular = None
    cond = None
    cov_full = np.full((p, p), np.nan)
    if np.any(keep):
        Jk = J[:, keep]
        H = Jk.T @ Jk
        try:
            Ck = np.linalg.inv(H) * phi
            singular = False
        except np.linalg.LinAlgError:
            Ck = np.linalg.pinv(H) * phi
            singular = True
        cond = float(np.linalg.cond(H))
        sig[keep] = np.sqrt(np.diag(Ck))
        idx = np.where(keep)[0]
        cov_full[np.ix_(idx, idx)] = Ck
    rho1 = float(np.sum(r[1:] * r[:-1]) / np.sum(r * r))
    ar1 = math.sqrt((1 + rho1) / (1 - rho1)) if rho1 > 0 else 1.0
    # SPEC 3.5 v1.1: AICc/BIC use chi2_d, the chi2 under the data-only (pass-1) variance
    vd = variance(smooth_quadratic(model.E, model.y), model.c)
    chi2_d = float(np.sum((model.y - mu) ** 2 / vd))
    AICc = N * math.log(chi2_d / N) + 2 * p + 2 * p * (p + 1) / (N - p - 1)
    BIC = N * math.log(chi2_d / N) + p * math.log(N)
    DW = float(np.sum(np.diff(r) ** 2) / np.sum(r * r))
    rz, nruns = runs_z(r)
    # chi2 re-evaluated with v = mu_final (diagnostic, Pearson with the final model)
    vf = variance(mu, model.c)
    chi2_mufinal = float(np.sum((model.y - mu) ** 2 / vf))
    params = {}
    sigmas = {}
    for (own, pn), val in vals.items():
        params.setdefault(own, {})[pn] = val
    for (own, pn), s in zip(model.free, sig):
        sigmas.setdefault(own, {})[pn] = None if not np.isfinite(s) else float(s)
    out = dict(params=params, sigma=sigmas,
               free=[f"{a}.{b}" for a, b in model.free],
               theta=theta.tolist(), sigma_vec=[None if not np.isfinite(s) else float(s) for s in sig],
               at_bound=[f"{a}.{b}" for (a, b), t in zip(model.free, atb) if t],
               chi2=chi2, chi2_d=chi2_d, phi=phi, N=N, p=p, rho1=rho1, ar1_factor=ar1,
               AICc=AICc, BIC=BIC, DW=DW, runs_z=rz, n_runs=nruns,
               chi2_with_v_mu_final=chi2_mufinal, hessian_singular=singular, hessian_cond=cond,
               mu=mu.tolist(), B=np.asarray(B).tolist(), v=v.tolist(),
               cov=[[None if not np.isfinite(c) else float(c) for c in row] for row in cov_full])
    if extra:
        out.update(extra)
    return out


# ---------------------------------------------------------------------------
# 4. Quantification
# ---------------------------------------------------------------------------
HNU_ALKA = 1486.6


def tpp2m(E, Nv, rho, M, Eg):
    """TPP-2M IMFP in Angstrom (E = kinetic energy in eV)."""
    E = np.asarray(E, float)
    Ep = 28.8 * math.sqrt(Nv * rho / M)
    beta = -0.10 + 0.944 / math.sqrt(Ep * Ep + Eg * Eg) + 0.069 * rho ** 0.1
    gamma = 0.191 * rho ** -0.5
    U = Nv * rho / M
    C = 1.97 - 0.91 * U
    D = 53.4 - 20.8 * U
    return E / (Ep * Ep * (beta * np.log(gamma * E) - C / E + D / E ** 2)), \
        dict(Ep=Ep, beta=beta, gamma=gamma, U=U, C=C, D=D)


def transmission(KE, PE, opt):
    t = opt["type"]
    if t == "none":
        return 1.0
    if t == "power":
        return (PE / KE) ** opt["x"]
    if t == "thermo":
        L = math.log(KE / PE)
        return math.exp(sum(a * L ** k for k, a in enumerate(opt["a"])))
    if t == "table":
        ke = np.log(np.asarray(opt["KE"], float))
        T = np.log(np.asarray(opt["T"], float))
        return float(np.exp(np.interp(math.log(KE), ke, T)))
    raise ValueError(t)


def imfp(KE, opt):
    t = opt["type"]
    if t == "none":
        return 1.0
    if t == "power":
        return KE ** opt.get("n", 0.6)
    if t == "tpp2m":
        return float(tpp2m(KE, opt["Nv"], opt["rho"], opt["M"], opt["Eg"])[0])
    raise ValueError(t)


def quantify(peaks, PE, Topt, Lopt, hnu=HNU_ALKA):
    """peaks: [{name, A, sigma_A, S, BE}] -> at% with delta-method sigma (areas independent)."""
    n = []
    sn = []
    rows = []
    for pk in peaks:
        KE = hnu - pk["BE"]
        T = transmission(KE, PE, Topt)
        lam = imfp(KE, Lopt)
        den = pk["S"] * T * lam
        n.append(pk["A"] / den)
        sn.append(pk.get("sigma_A", 0.0) / den)
        rows.append(dict(name=pk["name"], KE=KE, T=T, lam=lam, n=pk["A"] / den))
    n = np.array(n); sn = np.array(sn)
    tot = n.sum()
    at = 100 * n / tot
    # d at_i / d n_j = 100 (delta_ij tot - n_i)/tot^2
    Jm = 100 * (np.eye(len(n)) * tot - n[:, None]) / tot ** 2
    cov = Jm @ np.diag(sn ** 2) @ Jm.T
    for i, rw in enumerate(rows):
        rw["at"] = float(at[i]); rw["sigma_at"] = float(math.sqrt(cov[i, i]))
    return rows


# ---------------------------------------------------------------------------
# 6. Noise estimate
# ---------------------------------------------------------------------------

def noise_kappa(E, y):
    y = np.asarray(y, float)
    ys = smooth_quadratic(E, y)
    d = (y[1:-1] - 0.5 * (y[:-2] + y[2:])) / math.sqrt(1.5)
    return float(np.median(d * d / ys[1:-1]) / 0.4549)


def load_txt(path):
    d = np.loadtxt(path)
    o = np.argsort(d[:, 0])
    return d[o, 0], d[o, 1]


def region_slice(E, y, emin, emax):
    m = (E >= emin - 1e-9) & (E <= emax + 1e-9)
    return E[m], y[m]
