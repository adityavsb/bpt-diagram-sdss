"""
Error propagation, confidence intervals, and pipeline glitch vs. physical anomaly
scoring for SDSS DR17 emission-line galaxy diagnostics.
All flux errors are the MPA-JHU rescaled errors (see query_data.MPA_JHU_ERR_SCALE).
"""
import numpy as np
import pandas as pd
from scipy.stats import chi2

# 1-D marginal 95% half-width factor and 2-D 95% ellipse semi-axis factor.
K95_1D = 1.959964
K95_2D = float(np.sqrt(chi2.ppf(0.95, df=2)))  # = 2.4477

def compute_line_snr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes signal-to-noise ratio for all four diagnostic lines and derives
    summary SNR metrics (min, mean, max).
    """
    res = pd.DataFrame(index=df.index)
    # Errors are already MPA-JHU-rescaled by query_data.query_sdss_catalog.
    # Zero/negative errors are excluded by the SQL selection, so no clamping.
    res["snr_ha"] = df["h_alpha_flux"] / df["h_alpha_flux_err"]
    res["snr_hb"] = df["h_beta_flux"] / df["h_beta_flux_err"]
    res["snr_oiii"] = df["oiii_5007_flux"] / df["oiii_5007_flux_err"]
    res["snr_nii"] = df["nii_6584_flux"] / df["nii_6584_flux_err"]
    
    snr_matrix = res[["snr_ha", "snr_hb", "snr_oiii", "snr_nii"]].values
    res["snr_min"] = np.min(snr_matrix, axis=1)
    res["snr_mean"] = np.mean(snr_matrix, axis=1)
    res["snr_max"] = np.max(snr_matrix, axis=1)
    return res

def compute_error_propagation(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes 1st-order analytical Gaussian error propagation for logarithmic BPT ratios:
      x = log10([N II] 6584 / H-alpha)
      y = log10([O III] 5007 / H-beta)
      
    sigma_log10(A/B) = (1 / ln(10)) * sqrt( (sigma_A/A)^2 + (sigma_B/B)^2 )
    
    Derives 1-D marginal 95% half-widths (1.96 sigma) and 2-D 95% ellipse semi-axes (2.448 sigma, chi2_2).
    """
    ln10 = np.log(10.0)
    
    # Fluxes are > 0 by selection; errors are MPA-JHU-rescaled (see query_data).
    rel_err_nii = df["nii_6584_flux_err"] / df["nii_6584_flux"]
    rel_err_ha = df["h_alpha_flux_err"] / df["h_alpha_flux"]
    rel_err_oiii = df["oiii_5007_flux_err"] / df["oiii_5007_flux"]
    rel_err_hb = df["h_beta_flux_err"] / df["h_beta_flux"]
    
    sigma_x = (1.0 / ln10) * np.sqrt(rel_err_nii**2 + rel_err_ha**2)
    sigma_y = (1.0 / ln10) * np.sqrt(rel_err_oiii**2 + rel_err_hb**2)
    
    # 1-D marginal 95% half-widths (for the printed [lo, hi] strings)
    delta_x_95 = K95_1D * sigma_x
    delta_y_95 = K95_1D * sigma_y
    # 2-D 95% ellipse semi-axes (chi-square with 2 dof); errors treated as uncorrelated
    a_95 = K95_2D * sigma_x
    b_95 = K95_2D * sigma_y
    
    x = df["log_nii_ha"]
    y = df["log_oiii_hb"]
    
    res = pd.DataFrame(index=df.index)
    res["sigma_log_nii_ha"] = sigma_x
    res["sigma_log_oiii_hb"] = sigma_y
    res["delta_x_95"] = delta_x_95
    res["delta_y_95"] = delta_y_95
    res["ellipse_95_semi_a"] = a_95
    res["ellipse_95_semi_b"] = b_95
    
    res["ci_95_nii_ha_lower"] = x - delta_x_95
    res["ci_95_nii_ha_upper"] = x + delta_x_95
    res["ci_95_oiii_hb_lower"] = y - delta_y_95
    res["ci_95_oiii_hb_upper"] = y + delta_y_95
    
    # Format readable CI ranges
    res["ci_95_nii_ha_str"] = [
        f"[{lo:.2f}, {hi:.2f}]" for lo, hi in zip(res["ci_95_nii_ha_lower"], res["ci_95_nii_ha_upper"])
    ]
    res["ci_95_oiii_hb_str"] = [
        f"[{lo:.2f}, {hi:.2f}]" for lo, hi in zip(res["ci_95_oiii_hb_lower"], res["ci_95_oiii_hb_upper"])
    ]
    
    # 2-D 95% confidence ellipse area = pi * a * b
    res["ci_95_ellipse_area"] = np.pi * a_95 * b_95
    return res

def compute_fragility_and_outlier_scores(
    df: pd.DataFrame,
    consensus_anomaly_score: np.ndarray,
    snr_df: pd.DataFrame,
    err_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Both outputs are heuristics, not calibrated probabilities.
    measurement_fragility in [0,1] rises toward 1 as the minimum line S/N approaches
    the selection floor, as the 2-D error ellipse grows relative to the sample's 95th percentile,
    and as Hbeta/[O III] S/N falls.
    quality_weighted_outlier_score = consensus percentile x (1 - fragility).
    """
    snr_min = snr_df["snr_min"].values
    ellipse_area = err_df["ci_95_ellipse_area"].values
    
    # data-driven reference: the 95th percentile ellipse area in the sample
    AREA_REF = float(np.nanpercentile(ellipse_area, 95))
    
    # --- 1. MEASUREMENT FRAGILITY MODEL ---
    # a) SNR vulnerability factor: penalty near survey boundary S/N ~ 3.0
    # Sigmoidal-like dropoff from S/N = 3 (high risk) to S/N >= 12 (low risk)
    snr_risk = np.clip(1.0 - (snr_min - 3.0) / 10.0, 0.0, 1.0)
    
    # b) Ellipse uncertainty factor: large area means position on BPT diagram is blurry
    area_risk = np.clip(ellipse_area / AREA_REF, 0.0, 1.0)
    
    # c) Faint line vulnerability (H-beta and [O III] are commonly susceptible to sky subtraction jitter)
    faint_line_snr = np.minimum(snr_df["snr_hb"].values, snr_df["snr_oiii"].values)
    faint_risk = np.clip(1.0 - (faint_line_snr - 3.0) / 8.0, 0.0, 1.0)
    
    # Weighted composite fragility (0 - 1)
    raw_glitch_risk = 0.45 * snr_risk + 0.35 * area_risk + 0.20 * faint_risk
    
    # --- 2. QUALITY-WEIGHTED OUTLIER SCORE ---
    # The consensus percentile, discounted by measurement fragility. No further
    # S/N factors: fragility already encodes S/N and ellipse size, and brightness
    # must not be rewarded twice.
    measurement_fragility = raw_glitch_risk  # 0..1
    quality_weighted = consensus_anomaly_score * (1.0 - measurement_fragility)

    res = pd.DataFrame(index=df.index)
    res["measurement_fragility"] = np.round(measurement_fragility, 3)        # 0..1, heuristic
    res["quality_weighted_outlier_score"] = np.round(quality_weighted, 2)     # 0..100, rank-based
    return res

def compute_boundary_sigma_distance(df: pd.DataFrame, err_df: pd.DataFrame, n_grid: int = 4000) -> pd.DataFrame:
    """
    Minimum Mahalanobis-like distance (in sigma units, uncorrelated errors) from
    each object to the Kauffmann+03 and Kewley+01 demarcation curves:
        d^2 = ((x - xc)/sigma_x)^2 + ((y - yc)/sigma_y)^2, minimised over curve points.
    An object with n_sigma_nearest_boundary >= 3 cannot be moved across the
    nearest classification boundary by a 3-sigma fluctuation of its line ratios.
    This is the quantity the README calls "detachment from standard tracks".
    """
    from bpt_math import kauffmann_boundary, kewley_boundary

    xg_k03 = np.linspace(-3.0, 0.05 - 1e-4, n_grid)
    yg_k03 = kauffmann_boundary(xg_k03)
    xg_k01 = np.linspace(-3.0, 0.47 - 1e-4, n_grid)
    yg_k01 = kewley_boundary(xg_k01)
    # restrict to the plotted dynamic range so the asymptote does not dominate
    keep03 = np.isfinite(yg_k03) & (yg_k03 > -3) & (yg_k03 < 3)
    keep01 = np.isfinite(yg_k01) & (yg_k01 > -3) & (yg_k01 < 3)
    xg_k03, yg_k03 = xg_k03[keep03], yg_k03[keep03]
    xg_k01, yg_k01 = xg_k01[keep01], yg_k01[keep01]

    x = df["log_nii_ha"].values[:, None]
    y = df["log_oiii_hb"].values[:, None]
    sx = err_df["sigma_log_nii_ha"].values[:, None]
    sy = err_df["sigma_log_oiii_hb"].values[:, None]

    chunk_size = 2000
    d03_list, d01_list = [], []
    for i in range(0, len(df), chunk_size):
        xc = x[i : i + chunk_size]
        yc = y[i : i + chunk_size]
        sxc = sx[i : i + chunk_size]
        syc = sy[i : i + chunk_size]
        d03_chunk = np.sqrt(np.min(((xc - xg_k03) / sxc) ** 2 + ((yc - yg_k03) / syc) ** 2, axis=1))
        d01_chunk = np.sqrt(np.min(((xc - xg_k01) / sxc) ** 2 + ((yc - yg_k01) / syc) ** 2, axis=1))
        d03_list.append(d03_chunk)
        d01_list.append(d01_chunk)

    d03 = np.concatenate(d03_list)
    d01 = np.concatenate(d01_list)

    res = pd.DataFrame(index=df.index)
    res["n_sigma_kauffmann"] = d03
    res["n_sigma_kewley"] = d01
    res["n_sigma_nearest_boundary"] = np.minimum(d03, d01)
    res["class_robust_3sigma"] = res["n_sigma_nearest_boundary"] >= 3.0
    return res

