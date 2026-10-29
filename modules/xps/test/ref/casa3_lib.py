"""Shared simulator for the CasaXPS 2.3.27 Round-3 readings (SPEC_casa_v3.md).
Region i0 = 19, i1 = 361 (BE 280.03 / 297.13); Linear background, endpoints = mean of y over (i0+1)±1, (i1−1)±1;
unweighted LSQ over channels i0+1 … i1−1; RMS = sqrt(Σ r²).  Shapes written from Cookbook 2026 pp. 65–67 (no CasaXPS code).
Coordinates: t = (E − pos)/F, t > 0 = HIGHER BE (Cookbook x_c = −t).  Generic shape = base(t) [⊗ Gaussian FWHM s·F]."""
import os
import numpy as np
from scipy import integrate, special, interpolate
from scipy.optimize import least_squares
from scipy import signal

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', 'data', 'casaval')
LN2 = np.log(2); RHO = 1 / 144.87; I0, I1 = 19, 361
E = 279.08 + 0.05 * np.arange(381)
II = np.arange(I0 + 1, I1)


def load(name):
    L = open(os.path.join(D, name + '.vms')).read().split('\n'); k = L.index('762')
    return np.array([float(s) for s in L[k + 5:k + 5 + 762:2]])[::-1]


def linbg(y, w=1):
    a, b = I0 + 1, I1 - 1; lo = y[a - w:a + w + 1].mean(); hi = y[b - w:b + w + 1].mean()
    return lo + (hi - lo) * (E - E[a]) / (E[b] - E[a])


# ---- base functions (unit F), t>0 = high BE ----
def lor(t): return 1 / (1 + 4 * t * t)
def la(a, b): return lambda t: np.where(t >= 0, lor(t) ** a, lor(t) ** b)      # a: high-BE side (Cookbook x_c<=0)
def gl(m):  # Cookbook GL(x,m) product, m in 0..1
    return lambda t: np.exp(-4 * LN2 * (1 - m) * t * t) / (1 + 4 * m * t * t)
def sgl(m): return lambda t: (1 - m) * np.exp(-4 * LN2 * t * t) + m * lor(t)
def tla(a, mu, sgn=+1):  # T(x_c) = l^a [pi/2 - atan(2 x_c) + pi/mu], x_c = -sgn*t
    return lambda t: lor(t) ** a * (np.pi / 2 + np.arctan(2 * sgn * t) + np.pi / mu)
def ds(al, sgn=+1):      # P(x_c) = cos(pi a/2 + (1-a) atan(2 x_c)) / (1+4x_c^2)^((1-a)/2), x_c = -sgn*t
    return lambda t: np.cos(np.pi * al / 2 + (1 - al) * np.arctan(-2 * sgn * t)) / (1 + 4 * t * t) ** ((1 - al) / 2)


def quadR(f):
    return sum(integrate.quad(f, a, b, limit=1000, epsabs=0, epsrel=1e-12)[0] for a, b in [(-np.inf, -50), (-50, 0), (0, 50), (50, np.inf)])


def gconv(v, s, dt):
    if s <= 0: return v
    sg = s / (2 * np.sqrt(2 * LN2)); K = int(np.ceil(9 * sg / dt)); k = np.arange(-K, K + 1) * dt
    w = np.exp(-0.5 * (k / sg) ** 2); w /= w.sum(); return signal.fftconvolve(v, w, mode='same')


class Grid:
    """base sampled on a fine t-grid, optional Gaussian conv (FWHM s in F units), optional ST step; cubic interpolation.
    Area over R from quad of the base (convolution keeps it); ST area = A^2 (1/2 + 1/mu)."""
    def __init__(self, base, s=0.0, T=200.0, dt=0.0025, st_mu=None, st_s=0.0, A=None):
        t = np.arange(-T, T + dt / 2, dt); v = base(t)
        A0 = quadR(lambda u: float(base(np.array([u]))[0])) if A is None else A
        v = gconv(v, s, dt)
        if st_mu is not None:
            left = integrate.quad(lambda u: float(base(np.array([u]))[0]), -np.inf, -T, limit=500)[0]
            cum = left + np.concatenate([[0], np.cumsum(0.5 * (v[1:] + v[:-1]) * dt)])
            v = v * (cum + A0 / st_mu); A0 = A0 * A0 * (0.5 + 1 / st_mu)
            v = gconv(v, st_s, dt)
        self.t, self.v, self.A = t, v, A0
        self.f = interpolate.CubicSpline(t, v)

    def unit(self, x, pos, F):   # area-1 (over R) shape in eV
        return self.f((x - pos) / F) / (self.A * F)


def fwhm_t(g):   # width at half max of a Grid (t units)
    t, v = g.t, g.v; im = int(np.argmax(v)); h = v[im] / 2
    i = im
    while v[i] > h: i -= 1
    a = t[i] + (h - v[i]) / (v[i + 1] - v[i]) * (t[i + 1] - t[i])
    j = im
    while v[j] > h: j += 1
    b = t[j - 1] + (v[j - 1] - h) / (v[j - 1] - v[j]) * (t[j] - t[j - 1])
    return b - a, t[im], a, b


def fit(name, model, p0, chans=II):
    y = load(name); B = linbg(y); x = E[chans]; d = y[chans]; b = B[chans]
    r = least_squares(lambda p: b + model(x, *p) - d, p0, x_scale='jac', xtol=1e-14, ftol=1e-14, gtol=1e-14)
    res = d - b - model(x, *r.x)
    return r.x, float(np.sqrt((res ** 2).sum())), len(chans) - len(p0)


def fit_grid(name, g, F0=1.0, pos0=285.0, area0=15000.0):
    model = lambda x, pos, area, F: area * g.unit(x, pos, F)
    p, rms, dof = fit(name, model, [pos0 + 0.003, area0 * 0.97, F0 * 1.02])
    return dict(pos=p[0], area=p[1], F=p[2], rms=rms, dof=dof)
