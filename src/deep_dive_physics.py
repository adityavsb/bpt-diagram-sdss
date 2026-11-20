"""
Physics helpers for the LEDA 4572013 deep dive: dust correction, direct-Te
("electron temperature") abundances, and a small Gaussian line fitter.

Everything is plain numpy/scipy so it can be unit-tested (tests/test_deep_dive.py).

Equations and where they come from
  * Cardelli, Clayton & Mathis (1989) optical extinction polynomial, R_V = 3.1.
  * T_e([O III]) from R = I(4959+5007)/I(4363): Osterbrock (1989) Eq. 5.4,
        R = 7.90 exp(3.29e4 / T) / (1 + 4.5e-4 n_e / sqrt(T)).
    Cross-check: the empirical fit of Perez-Montero (2017), t3 = 0.7840 - 0.0001357 R + 48.44 / R.
  * Ionic abundances O++/H+, O+/H+ and N+/O+ with the T_e(O II)-T_e(O III) relation:
    Izotov, Stasinska, Meynet, Guseva & Thuan (2006, A&A 448, 955). The coefficients were
    written from memory of that paper; the notebook validates them empirically (recovers the
    N/O plateau and the MPA-JHU O/H on ~800 control galaxies) instead of assuming they are right.
"""
import numpy as np
from scipy.optimize import brentq, curve_fit

C_KMS = 299792.458

# rest wavelengths (Angstrom, air)
LINES = {
    "oii_3727": 3727.4,   # centroid of the 3726+3729 doublet
    "hgamma": 4340.47,
    "oiii_4363": 4363.21,
    "hbeta": 4861.33,
    "oiii_4959": 4958.91,
    "oiii_5007": 5006.84,
    "halpha": 6562.80,
    "nii_6584": 6583.45,
}


# ------------------------------------------------------------------ dust
def ccm89_k(wave_angstrom, r_v: float = 3.1):
    """A(lambda)/E(B-V) from the Cardelli+89 optical polynomial (valid 3030-9091 A)."""
    x = 1e4 / np.asarray(wave_angstrom, dtype=float)
    if np.any((x < 1.1) | (x > 3.3)):
        raise ValueError("wavelength outside the CCM89 optical range (3030-9091 A)")
    y = x - 1.82
    a = (1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4
         + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7)
    b = (1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4
         - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7)
    return r_v * (a + b / r_v)


def ebv_from_balmer(f_ha, f_hb, intrinsic: float = 2.86, r_v: float = 3.1):
    """E(B-V) from the Halpha/Hbeta decrement (case B, T=1e4 K -> 2.86). Floored at 0."""
    kdiff = ccm89_k(LINES["hbeta"], r_v) - ccm89_k(LINES["halpha"], r_v)
    ebv = 2.5 / kdiff * np.log10((np.asarray(f_ha) / np.asarray(f_hb)) / intrinsic)
    return np.maximum(ebv, 0.0)


def deredden(flux, wave_angstrom, ebv, r_v: float = 3.1):
    return np.asarray(flux) * 10 ** (0.4 * ccm89_k(wave_angstrom, r_v) * ebv)


# ------------------------------------------------------------------ temperatures
def oiii_ratio_R(f4959, f5007, f4363):
    return (np.asarray(f4959) + np.asarray(f5007)) / np.asarray(f4363)


def te_oiii(R: float, ne: float = 100.0) -> float:
    """Invert Osterbrock's relation for T_e([O III]) in Kelvin. NaN outside 5e3-5e4 K."""
    def f(T):
        return 7.90 * np.exp(3.29e4 / T) / (1 + 4.5e-4 * ne / np.sqrt(T)) - R
    try:
        return brentq(f, 5e3, 5e4)
    except ValueError:
        return np.nan


def t3_empirical(R):
    """Perez-Montero (2017) empirical t3 = T_e([O III]) / 1e4 K."""
    R = np.asarray(R, dtype=float)
    return 0.7840 - 0.0001357 * R + 48.44 / R


def t2_from_t3(t3):
    """T_e([O II]) / 1e4 K from T_e([O III]) / 1e4 K (Izotov+06 relation for low-metallicity H II regions)."""
    t3 = np.asarray(t3, dtype=float)
    return -0.744 + t3 * (2.338 - 0.610 * t3)


# ------------------------------------------------------------------ abundances
def o_plus2_over_h(i4959_plus_5007_over_hb, t3):
    """12 + log(O++/H+)"""
    t3 = np.asarray(t3, dtype=float)
    return (np.log10(i4959_plus_5007_over_hb) + 6.1868 + 1.2491 / t3 - 0.5816 * np.log10(t3))


def o_plus_over_h(i3727_over_hb, t2, ne: float = 100.0):
    """12 + log(O+/H+)"""
    t2 = np.asarray(t2, dtype=float)
    return (np.log10(i3727_over_hb) + 5.887 + 1.641 / t2 - 0.543 * np.log10(t2)
            + 0.000114 * ne * t2)


def n_plus_over_o_plus(i6548_plus_6584_over_i3727, t2):
    """log(N+/O+); adopted as log(N/O) (ICF = 1)."""
    t2 = np.asarray(t2, dtype=float)
    return (np.log10(i6548_plus_6584_over_i3727) + 0.307 - 0.02 * np.log10(t2) - 0.726 / t2)


def direct_abundances(f, ebv=None, ne: float = 100.0):
    """
    Direct-Te abundances from a dict of observed (reddened) line fluxes with keys
    oii_3727 (3726+3729 summed), hbeta, halpha, oiii_4363, oiii_4959, oiii_5007, nii_6584.
    If ebv is None it is derived from Halpha/Hbeta. Returns a dict (scalars).
    """
    if ebv is None:
        ebv = float(ebv_from_balmer(f["halpha"], f["hbeta"]))
    d = {k: float(deredden(f[k], LINES[k], ebv)) for k in
         ["oii_3727", "hbeta", "oiii_4363", "oiii_4959", "oiii_5007", "nii_6584"]}
    R = oiii_ratio_R(d["oiii_4959"], d["oiii_5007"], d["oiii_4363"])
    T3 = te_oiii(R, ne)
    t3 = T3 / 1e4
    t2 = float(t2_from_t3(t3))
    opp = float(o_plus2_over_h((d["oiii_4959"] + d["oiii_5007"]) / d["hbeta"], t3))
    op = float(o_plus_over_h(d["oii_3727"] / d["hbeta"], t2, ne))
    oh = 12 + np.log10(10 ** (opp - 12) + 10 ** (op - 12)) if np.isfinite(opp + op) else np.nan
    # [N II] 6548 = 6584 / 3.05 (atomic physics), so 6548+6584 = 6584 * (1 + 1/3.05)
    no = float(n_plus_over_o_plus(d["nii_6584"] * (1 + 1 / 3.05) / d["oii_3727"], t2))
    return {"ebv": float(ebv), "R": float(R), "Te_OIII": float(T3), "Te_OII": t2 * 1e4,
            "t3_empirical_PM17": float(t3_empirical(R)),
            "OPP_H": opp, "OP_H": op, "OH": float(oh), "logNO": no,
            "O2_frac": float(10 ** (op - opp) / (1 + 10 ** (op - opp)))}


def direct_abundances_mc(f, f_err, n=2000, err_inflate=1.0, ne=100.0, seed=0):
    """Monte-Carlo the direct abundances over line-flux uncertainties; returns (median, 16th, 84th) per quantity."""
    rng = np.random.default_rng(seed)
    keys = [k for k in f if k in f_err]
    out = {}
    for _ in range(n):
        draw = {k: max(rng.normal(f[k], err_inflate * f_err[k]), 1e-6) for k in keys}
        r = direct_abundances(draw, ne=ne)
        for q in ["Te_OIII", "OH", "logNO", "ebv"]:
            out.setdefault(q, []).append(r[q])
    return {q: tuple(np.nanpercentile(v, [50, 16, 84])) for q, v in out.items()}


# ------------------------------------------------------------------ line fitting
def gauss(x, amp, mu, sig):
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def fit_multiline(wave, flux, ivar, centers, window, tie=None, sig0=1.3, sig_bounds=(0.6, 4.0)):
    """
    Joint fit of several blended/neighbouring lines in one window (rest frame, Angstrom).

    centers : list of rest wavelengths
    window  : (lo, hi) rest-wavelength window
    tie     : {i: (parent_index, ratio)} makes amplitude_i = ratio * amplitude_parent
              (e.g. [N II] 6548 = 6584 / 3.05, [O III] 4959 = 5007 / 2.98)
    All lines share one velocity dispersion and one velocity offset, and sit on a linear
    continuum. Returns a list of dicts (flux, flux_err, snr) aligned with `centers`, plus
    the fit's reduced chi2 under key 'chi2' of the first dict.
    The error bars come from the covariance matrix, scaled by the reduced chi2 if it exceeds 1.
    """
    tie = tie or {}
    m = (wave > window[0]) & (wave < window[1]) & (ivar > 0)
    x, y, s = wave[m], flux[m], 1 / np.sqrt(ivar[m])
    n = len(centers)
    free = [i for i in range(n) if i not in tie]
    pivot = float(np.mean(window))

    def amps_from(p):
        a = np.zeros(n)
        for j, i in enumerate(free):
            a[i] = p[4 + j]
        for i, (par, ratio) in tie.items():
            a[i] = ratio * a[par]
        return a

    def model(x, *p):
        sig, dmu, c0, c1 = p[:4]
        a = amps_from(p)
        out = c0 + c1 * (x - pivot)
        for i, c in enumerate(centers):
            out = out + gauss(x, a[i], c + dmu, sig)
        return out

    c0 = float(np.median(y))
    p0 = [sig0, 0.0, c0, 0.0] + [max(float(np.max(y) - c0) * 0.3, 1e-3)] * len(free)
    lo = [sig_bounds[0], -3.0, -np.inf, -np.inf] + [0.0] * len(free)
    hi = [sig_bounds[1], 3.0, np.inf, np.inf] + [np.inf] * len(free)
    try:
        p, cov = curve_fit(model, x, y, p0=p0, sigma=s, absolute_sigma=True, bounds=(lo, hi), maxfev=40000)
    except Exception:
        return [{"flux": np.nan, "flux_err": np.nan, "snr": np.nan, "ok": False} for _ in centers]
    chi2 = float(np.sum(((y - model(x, *p)) / s) ** 2) / max(len(x) - len(p), 1))
    cov = cov * max(chi2, 1.0)
    sig = p[0]
    k = sig * np.sqrt(2 * np.pi)
    out = []
    for i in range(n):
        # flux_i = c_i * amp_parent * sig * sqrt(2 pi); gradient w.r.t. (sig, amp_parent)
        par, ratio = tie.get(i, (i, 1.0))
        j = 4 + free.index(par)
        amp = ratio * p[j]
        J = np.zeros(len(p))
        J[0] = amp * np.sqrt(2 * np.pi)
        J[j] = ratio * k
        err = float(np.sqrt(max(J @ cov @ J, 0.0)))
        flx = float(amp * k)
        out.append({"flux": flx, "flux_err": err, "snr": flx / err if err > 0 else np.nan,
                    "ok": True, "sigma": float(sig), "dmu": float(p[1])})
    out[0]["chi2"] = chi2
    out[0]["x"], out[0]["yfit"] = x, model(x, *p)   # for plotting the fit over the data
    return out


def fit_line(wave, flux, ivar, center, window=30.0, sig0=2.5, sig_bounds=(0.8, 8.0)):
    """
    Fit one emission line with a Gaussian plus linear continuum in +-window Angstrom of `center`.
    Returns dict(flux, flux_err, mu, sigma, snr, ok). Flux is the Gaussian integral.
    Errors come from the covariance matrix scaled by reduced chi2 (>= 1).
    """
    m = (wave > center - window) & (wave < center + window) & (ivar > 0)
    if m.sum() < 12:
        return {"flux": np.nan, "flux_err": np.nan, "mu": np.nan, "sigma": np.nan, "snr": np.nan, "ok": False}
    x, y, s = wave[m], flux[m], 1 / np.sqrt(ivar[m])

    def model(x, amp, mu, sig, c0, c1):
        return gauss(x, amp, mu, sig) + c0 + c1 * (x - center)

    c0 = float(np.median(y[np.abs(x - center) > 12])) if (np.abs(x - center) > 12).sum() > 3 else float(np.median(y))
    p0 = [max(y[np.argmin(np.abs(x - center))] - c0, 1e-3), center, sig0, c0, 0.0]
    lo = [0, center - 4, sig_bounds[0], -np.inf, -np.inf]
    hi = [np.inf, center + 4, sig_bounds[1], np.inf, np.inf]
    try:
        p, cov = curve_fit(model, x, y, p0=p0, sigma=s, absolute_sigma=True, bounds=(lo, hi), maxfev=20000)
    except Exception:
        return {"flux": np.nan, "flux_err": np.nan, "mu": np.nan, "sigma": np.nan, "snr": np.nan, "ok": False}
    chi2 = np.sum(((y - model(x, *p)) / s) ** 2) / max(len(x) - 5, 1)
    scale = max(chi2, 1.0)
    amp, mu, sig = p[:3]
    flx = amp * sig * np.sqrt(2 * np.pi)
    # propagate amp,sig covariance into the integral
    J = np.array([sig * np.sqrt(2 * np.pi), 0.0, amp * np.sqrt(2 * np.pi)])
    C = cov[np.ix_([0, 1, 2], [0, 1, 2])] * scale
    err = float(np.sqrt(max(J @ C @ J, 0.0)))
    return {"flux": float(flx), "flux_err": err, "mu": float(mu), "sigma": float(sig),
            "snr": float(flx / err) if err > 0 else np.nan, "ok": True, "chi2": float(chi2)}
