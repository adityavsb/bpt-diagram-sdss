import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from bpt_math import kauffmann_boundary, kewley_boundary, classify_bpt
from glitch_filter import compute_error_propagation, compute_boundary_sigma_distance, K95_2D
from anomaly_detector import to_percentile_score


def test_boundaries_known_values():
    assert np.isclose(kauffmann_boundary(np.array([-0.5]))[0], 0.61 / (-0.55) + 1.30)
    assert np.isclose(kewley_boundary(np.array([0.0]))[0], 0.61 / (-0.47) + 1.19)
    assert np.isnan(kauffmann_boundary(np.array([0.05]))[0])
    assert np.isnan(kewley_boundary(np.array([0.47]))[0])


def test_classify_regions():
    x = np.array([-1.0, -0.3, 0.5, -0.2])
    y = np.array([-0.5, 0.0, 0.0, 0.3])
    assert list(classify_bpt(x, y)) == ["Star-Forming", "Composite", "AGN", "AGN"]


def _four_lines(flux=100.0, err=10.0):
    d = {}
    for l in ["h_alpha", "h_beta", "oiii_5007", "nii_6584"]:
        d[f"{l}_flux"] = [flux]
        d[f"{l}_flux_err"] = [err]
    d["log_nii_ha"] = [0.0]
    d["log_oiii_hb"] = [0.0]
    return pd.DataFrame(d)


def test_error_propagation_and_ellipse():
    e = compute_error_propagation(_four_lines())
    assert np.isclose(e["sigma_log_nii_ha"][0], np.sqrt(0.02) / np.log(10), atol=1e-6)
    assert np.isclose(K95_2D, 2.4477, atol=1e-3)
    assert np.isclose(e["ellipse_95_semi_a"][0], K95_2D * e["sigma_log_nii_ha"][0])


def test_boundary_sigma_distance():
    df = pd.DataFrame({"log_nii_ha": [-1.0], "log_oiii_hb": [-0.5]})
    err = pd.DataFrame({"sigma_log_nii_ha": [0.02], "sigma_log_oiii_hb": [0.02]})
    r = compute_boundary_sigma_distance(df, err)
    # vertical distance to Kauffmann at x=-1 is |-0.5 - 0.719| = 1.219 -> ~61 sigma; the
    # minimum over the curve is <= that and > 30 sigma.
    assert 30 < r["n_sigma_kauffmann"][0] <= 1.219 / 0.02 + 1e-6
    assert r["class_robust_3sigma"][0]


def test_percentile_monotonic():
    s = to_percentile_score(np.array([3.0, 1.0, 2.0]))
    assert list(np.argsort(s)) == [1, 2, 0] and s.max() == 100.0
