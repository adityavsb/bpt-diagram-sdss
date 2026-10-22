"""
Tri-algorithmic anomaly detection and voting ensemble for SDSS DR17 galaxies.
Combines:
1. Isolation Forest (global path-length isolation)
2. Local Outlier Factor (local neighborhood density ratio)
3. PCA Reconstruction Error (spectral manifold reconstruction residual)
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler

def extract_anomaly_features(df: pd.DataFrame) -> np.ndarray:
    """
    Constructs a standardized feature representation from emission line diagnostics:
    - BPT coordinates: log([N II]/Hα), log([O III]/Hβ)
    - Dust extinction tracker: log(Hα / Hβ) [Balmer decrement]
    - Forbidden line ratio: log([O III] / [N II])
    - Log line luminosities / fluxes for dynamic range
    """
    ha = np.maximum(df["h_alpha_flux"].values, 1e-4)
    hb = np.maximum(df["h_beta_flux"].values, 1e-4)
    oiii = np.maximum(df["oiii_5007_flux"].values, 1e-4)
    nii = np.maximum(df["nii_6584_flux"].values, 1e-4)
    
    log_nii_ha = df["log_nii_ha"].values
    log_oiii_hb = df["log_oiii_hb"].values
    log_balmer = np.log10(ha / hb)
    log_oiii_nii = np.log10(oiii / nii)
    log_ha = np.log10(ha)
    log_oiii = np.log10(oiii)
    
    features = np.column_stack([
        log_nii_ha,
        log_oiii_hb,
        log_balmer,
        log_oiii_nii,
        log_ha,
        log_oiii
    ])
    
    scaler = RobustScaler()
    return scaler.fit_transform(features)

def to_percentile_score(raw_scores: np.ndarray) -> np.ndarray:
    """
    Transforms raw anomaly scores into a clean empirical percentile scale [0.0, 100.0].
    Higher value = more anomalous.
    """
    ranks = rankdata(raw_scores, method="average")
    return np.round((ranks / len(raw_scores)) * 100.0, 2)

def detect_anomalies_tri_ensemble(df: pd.DataFrame, random_state: int = 42) -> pd.DataFrame:
    """
    Runs the Tri-Algorithmic Voting Ensemble on SDSS DR17 galaxies:
    1. Isolation Forest
    2. Local Outlier Factor (LOF)
    3. PCA Reconstruction Error
    
    Returns an enriched DataFrame with individual calibrated scores (0-100),
    equal-weighted consensus anomaly score, algorithm vote count (0 to 3),
    and astrophysical archetype classification.
    """
    X = extract_anomaly_features(df)
    n_samples = len(df)
    
    # --- Algorithm 1: Isolation Forest ---
    iso = IsolationForest(
        n_estimators=250,
        contamination=0.03,
        max_samples="auto",
        random_state=random_state,
        n_jobs=-1
    )
    iso.fit(X)
    # decision_function gives negative values for anomalies -> invert
    raw_iforest = -iso.decision_function(X)
    score_iforest = to_percentile_score(raw_iforest)
    
    # --- Algorithm 2: Local Outlier Factor ---
    lof = LocalOutlierFactor(
        n_neighbors=35,
        contamination=0.03,
        n_jobs=-1
    )
    lof.fit_predict(X)
    # negative_outlier_factor_ gives large negative numbers for outliers -> invert
    raw_lof = -lof.negative_outlier_factor_
    score_lof = to_percentile_score(raw_lof)
    
    # --- Algorithm 3: PCA Reconstruction Error ---
    # Low-rank projection to capture 90%+ normal spectral manifold
    pca = PCA(n_components=2, random_state=random_state)
    X_pca = pca.fit_transform(X)
    X_recon = pca.inverse_transform(X_pca)
    raw_pca_err = np.linalg.norm(X - X_recon, axis=1)
    score_pca = to_percentile_score(raw_pca_err)
    
    # --- Equal-Weighted Consensus Scoring ---
    consensus_score = np.round((score_iforest + score_lof + score_pca) / 3.0, 2)
    
    # --- Voting Agreement (Algorithm Votes in Top 5%) ---
    vote_iforest = (score_iforest >= 95.0).astype(int)
    vote_lof = (score_lof >= 95.0).astype(int)
    vote_pca = (score_pca >= 95.0).astype(int)
    total_votes = vote_iforest + vote_lof + vote_pca
    
    results = pd.DataFrame(index=df.index)
    results["score_iforest"] = score_iforest
    results["score_lof"] = score_lof
    results["score_pca"] = score_pca
    results["consensus_anomaly_score"] = consensus_score
    results["algorithm_votes"] = total_votes
    
    return results

def assign_physical_archetype(
    log_nii_ha: np.ndarray,
    log_oiii_hb: np.ndarray,
    glitch_risk: np.ndarray,
    phys_conf: np.ndarray
) -> list:
    """
    Assigns an astrophysical archetype or pipeline glitch label based on
    diagnostic position, physical confidence, and glitch risk.
    """
    archetypes = []
    for x, y, gr, pc in zip(log_nii_ha, log_oiii_hb, glitch_risk, phys_conf):
        if gr >= 75.0 and pc < 35.0:
            archetypes.append("Suspected Pipeline Glitch / Measurement Artifact")
        elif y >= 0.70 and x <= -0.70:
            archetypes.append("Extreme Green Pea / EELG Candidate")
        elif y >= 0.80 and x > 0.0:
            archetypes.append("Extreme Seyfert / Fast Accretion Outflow")
        elif y <= -0.45 and x >= 0.0:
            archetypes.append("Under-Ionized LINER / Shock Front")
        elif x >= 0.35:
            archetypes.append("Nitrogen-Enriched Super-Solar System")
        elif -0.30 <= x <= 0.10 and -0.10 <= y <= 0.45:
            archetypes.append("Transition Composite / Intermediate Outlier")
        elif y >= 0.90:
            archetypes.append("Ultra-High Ionization Outlier")
        elif x <= -1.10:
            archetypes.append("Extremely Metal-Poor Galaxy (XMPG)")
        else:
            archetypes.append("Unconventional Emission Outlier")
    return archetypes
