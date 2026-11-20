import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from deep_dive_physics import (
    LINES, ccm89_k, ebv_from_balmer, deredden, te_oiii, t3_empirical,
    direct_abundances, fit_line, fit_multiline, gauss,
)


def test_ccm89_known_values():
    # standard R_V = 3.1 values: k(Halpha) ~ 2.53, k(Hbeta) ~ 3.61
    assert np.isclose(ccm89_k(LINES["halpha"]), 2.53, atol=0.02)
    assert np.isclose(ccm89_k(LINES["hbeta"]), 3.61, atol=0.02)


def test_dereddening_restores_case_b():
    ebv = float(ebv_from_balmer(3.4, 1.0))
    assert ebv > 0
    ha = deredden(3.4, LINES["halpha"], ebv)
    hb = deredden(1.0, LINES["hbeta"], ebv)
    assert np.isclose(ha / hb, 2.86, rtol=1e-6)
    assert float(ebv_from_balmer(2.5, 1.0)) == 0.0  # never negative


def test_te_oiii_round_trip():
    T, ne = 12000.0, 100.0
    R = 7.90 * np.exp(3.29e4 / T) / (1 + 4.5e-4 * ne / np.sqrt(T))
    assert np.isclose(te_oiii(R, ne), T, rtol=1e-6)
    # the empirical fit agrees with the physical relation to a few percent near 1e4 K
    assert abs(float(t3_empirical(R)) - T / 1e4) / (T / 1e4) < 0.05


def test_fit_line_recovers_flux():
    rng = np.random.default_rng(1)
    wave = np.arange(4800.0, 4920.0, 1.0)
    sigma_noise = 0.5
    flux = gauss(wave, 20.0, 4861.3, 2.4) + 3.0 + rng.normal(0, sigma_noise, wave.size)
    ivar = np.full(wave.size, 1 / sigma_noise**2)
    r = fit_line(wave, flux, ivar, 4861.3)
    true = 20.0 * 2.4 * np.sqrt(2 * np.pi)
    assert r["ok"] and abs(r["flux"] - true) < 4 * r["flux_err"]


def test_fit_multiline_deblends_and_respects_ties():
    rng = np.random.default_rng(2)
    wave = np.arange(4300.0, 4400.0, 0.7)
    noise = 0.4
    # Hgamma (strong) next to a weak [O III] 4363, as in the real spectrum
    flux = (gauss(wave, 12.0, 4340.47, 1.3) + gauss(wave, 2.0, 4363.21, 1.3) + 5.0
            + rng.normal(0, noise, wave.size))
    ivar = np.full(wave.size, 1 / noise**2)
    r = fit_multiline(wave, flux, ivar, [4340.47, 4363.21], (4310, 4390))
    true = 2.0 * 1.3 * np.sqrt(2 * np.pi)
    assert abs(r[1]["flux"] - true) < 4 * r[1]["flux_err"]

    # tied doublet: 4959 = 5007 / 2.98
    w2 = np.arange(4930.0, 5040.0, 0.7)
    f2 = gauss(w2, 30.0 / 2.98, 4958.91, 1.3) + gauss(w2, 30.0, 5006.84, 1.3) + 2.0 + rng.normal(0, 0.3, w2.size)
    r2 = fit_multiline(w2, f2, np.full(w2.size, 1 / 0.3**2), [5006.84, 4958.91], (4935, 5035),
                       tie={1: (0, 1 / 2.98)})
    assert np.isclose(r2[0]["flux"] / r2[1]["flux"], 2.98, rtol=1e-6)


def test_target_direct_abundance_regression():
    # MPA-JHU fluxes of LEDA 4572013 (specObjID 2237320683103217664), 1e-17 cgs
    f = {"oii_3727": 559.3137 + 635.8513, "hbeta": 582.2193, "halpha": 1811.066,
         "oiii_4363": 31.37065, "oiii_4959": 822.5381, "oiii_5007": 2372.394, "nii_6584": 47.86558}
    r = direct_abundances(f)
    assert 7.8 < r["OH"] < 8.2
    assert 11000 < r["Te_OIII"] < 14500
    assert -1.8 < r["logNO"] < -1.2
