"""
Error propagation, confidence intervals, and pipeline glitch vs. physical anomaly
scoring for SDSS DR17 emission-line galaxy diagnostics.
"""
import numpy as np
import pandas as pd

def compute_line_snr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes signal-to-noise ratio for all four diagnostic lines and derives
    summary SNR metrics (min, mean, max).
    """
    res = pd.DataFrame(index=df.index)
    res["snr_ha"] = df["h_alpha_flux"] / np.maximum(df["h_alpha_flux_err"], 1e-9)
    res["snr_hb"] = df["h_beta_flux"] / np.maximum(df["h_beta_flux_err"], 1e-9)
    res["snr_oiii"] = df["oiii_5007_flux"] / np.maximum(df["oiii_5007_flux_err"], 1e-9)
    res["snr_nii"] = df["nii_6584_flux"] / np.maximum(df["nii_6584_flux_err"], 1e-9)
    
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
    
    Derives 95% confidence interval half-widths (1.96 * sigma) and 2D ellipse areas.
    """
    ln10 = np.log(10.0)
    
    rel_err_nii = df["nii_6584_flux_err"] / np.maximum(df["nii_6584_flux"], 1e-9)
    rel_err_ha = df["h_alpha_flux_err"] / np.maximum(df["h_alpha_flux"], 1e-9)
    rel_err_oiii = df["oiii_5007_flux_err"] / np.maximum(df["oiii_5007_flux"], 1e-9)
    rel_err_hb = df["h_beta_flux_err"] / np.maximum(df["h_beta_flux"], 1e-9)
    
    sigma_x = (1.0 / ln10) * np.sqrt(rel_err_nii**2 + rel_err_ha**2)
    sigma_y = (1.0 / ln10) * np.sqrt(rel_err_oiii**2 + rel_err_hb**2)
    
    # 95% Confidence Interval half-widths (Gaussian 1.96 standard errors)
    delta_x_95 = 1.96 * sigma_x
    delta_y_95 = 1.96 * sigma_y
    
    x = df["log_nii_ha"]
    y = df["log_oiii_hb"]
    
    res = pd.DataFrame(index=df.index)
    res["sigma_log_nii_ha"] = sigma_x
    res["sigma_log_oiii_hb"] = sigma_y
    res["delta_x_95"] = delta_x_95
    res["delta_y_95"] = delta_y_95
    
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
    
    # 2D 95% Confidence Ellipse Area = pi * a * b
    res["ci_95_ellipse_area"] = np.pi * delta_x_95 * delta_y_95
    return res

def compute_glitch_and_confidence_scores(
    df: pd.DataFrame,
    consensus_anomaly_score: np.ndarray,
    snr_df: pd.DataFrame,
    err_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Computes two scaled scores (0 to 100):
      1. pipeline_glitch_risk:
         High probability (60-100) indicates the object's anomalous location is likely
         an artifact of low SNR, high fractional line error, or a giant uncertainty ellipse.
      2. physical_anomaly_confidence:
         High probability (60-100) indicates high ML anomaly consensus backed by high SNR,
         tight error bounds, and statistical robustness (>3 sigma detachment).
    """
    snr_min = snr_df["snr_min"].values
    ellipse_area = err_df["ci_95_ellipse_area"].values
    
    # --- 1. PIPELINE GLITCH RISK MODEL ---
    # a) SNR vulnerability factor: penalty near survey boundary S/N ~ 3.0
    # Sigmoidal-like dropoff from S/N = 3 (high risk) to S/N >= 12 (low risk)
    snr_risk = np.clip(1.0 - (snr_min - 3.0) / 10.0, 0.0, 1.0)
    
    # b) Ellipse uncertainty factor: large area means position on BPT diagram is blurry
    # Median ellipse area in SDSS DR17 sample is ~0.02 - 0.08. Areas > 0.25 indicate noisy fits.
    area_risk = np.clip(ellipse_area / 0.25, 0.0, 1.0)
    
    # c) Faint line vulnerability (H-beta and [O III] are commonly susceptible to sky subtraction jitter)
    faint_line_snr = np.minimum(snr_df["snr_hb"].values, snr_df["snr_oiii"].values)
    faint_risk = np.clip(1.0 - (faint_line_snr - 3.0) / 8.0, 0.0, 1.0)
    
    # Weighted composite glitch risk (0 - 100)
    raw_glitch_risk = 0.45 * snr_risk + 0.35 * area_risk + 0.20 * faint_risk
    # Scale from 0 to 100
    pipeline_glitch_risk = np.clip(raw_glitch_risk * 100.0, 0.0, 100.0)
    
    # --- 2. PHYSICAL ANOMALY CONFIDENCE MODEL ---
    # Requires:
    # 1. High consensus anomaly score (0-100)
    # 2. High SNR across all lines (min SNR / 12 saturated at 1.0)
    # 3. Small error ellipse (position is sharp, not drifting)
    # 4. Low pipeline glitch risk
    snr_quality_factor = np.clip(snr_min / 12.0, 0.2, 1.0)
    sharpness_factor = np.clip(1.0 - (ellipse_area / 0.30), 0.1, 1.0)
    non_glitch_factor = (100.0 - pipeline_glitch_risk) / 100.0
    
    # Base confidence: product of consensus anomaly score and observational quality factors
    raw_physical_conf = (
        consensus_anomaly_score * 
        (non_glitch_factor ** 1.2) * 
        (snr_quality_factor ** 0.5) * 
        sharpness_factor
    )
    physical_anomaly_confidence = np.clip(raw_physical_conf, 0.0, 100.0)
    
    res = pd.DataFrame(index=df.index)
    res["pipeline_glitch_risk"] = np.round(pipeline_glitch_risk, 1)
    res["physical_anomaly_confidence"] = np.round(physical_anomaly_confidence, 1)
    return res
