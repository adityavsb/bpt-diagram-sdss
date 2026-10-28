"""
Tri-algorithmic anomaly detection and voting ensemble for SDSS DR17 galaxies.
Combines:
1. Isolation Forest (global path-length isolation)
2. Local Outlier Factor (local neighborhood density ratio)
3. PCA Reconstruction Error (spectral manifold reconstruction residual)
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler

FEATURE_NAMES = ["log_nii_ha", "log_oiii_hb", "log_balmer_decrement", "log_ew_ha"]

def extract_anomaly_features(df: pd.DataFrame):
    """
    Distance- and luminosity-independent feature set (4 features, full rank):
      1. log10([N II] 6584 / Halpha)        - BPT x
      2. log10([O III] 5007 / Hbeta)        - BPT y
      3. log10(Halpha / Hbeta)              - Balmer decrement (dust / case-B deviation)
      4. log10(EW(Halpha) / Angstrom)       - WHAN-style ionisation strength
                                              (MPA-JHU stores emission EW as negative)
    Absolute fluxes are deliberately excluded: in a flux-limited survey they
    track distance and luminosity, not emission-line physics.

    Returns (X_scaled, X_raw, scaler).
    """
    ha = df["h_alpha_flux"].values
    hb = df["h_beta_flux"].values
    ew_ha = -df["h_alpha_eqw"].values  # emission is negative in MPA-JHU
    n_clipped = int((ew_ha <= 0.1).sum())
    if n_clipped:
        print(f"[Features] {n_clipped} objects have EW(Halpha) <= 0.1 A; clipped to 0.1 A.")
    ew_ha = np.maximum(ew_ha, 0.1)

    X_raw = np.column_stack([
        df["log_nii_ha"].values,
        df["log_oiii_hb"].values,
        np.log10(ha / hb),
        np.log10(ew_ha),
    ])
    scaler = RobustScaler()
    X = scaler.fit_transform(X_raw)
    assert np.linalg.matrix_rank(X) == X.shape[1], "feature matrix is rank-deficient"
    return X, X_raw, scaler

def to_percentile_score(raw_scores: np.ndarray) -> np.ndarray:
    """
    Transforms raw anomaly scores into a clean empirical percentile scale [0.0, 100.0].
    Higher value = more anomalous.
    """
    ranks = rankdata(raw_scores, method="average")
    return np.round((ranks / len(raw_scores)) * 100.0, 2)

def detect_anomalies_tri_ensemble(
    df: pd.DataFrame,
    random_state: int = 42,
    lof_k: int = 35,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Runs the Tri-Algorithmic Voting Ensemble on SDSS DR17 galaxies:
    1. Isolation Forest
    2. Local Outlier Factor (LOF)
    3. PCA Reconstruction Error
    
    Returns an enriched DataFrame with individual calibrated scores (0-100),
    equal-weighted consensus anomaly score, algorithm vote count (0 to 3),
    and astrophysical archetype classification.
    """
    X, X_raw, scaler = extract_anomaly_features(df)
    n_samples = len(df)
    
    # --- Algorithm 1: Isolation Forest ---
    iso = IsolationForest(
        n_estimators=250,
        contamination="auto",  # contamination only affects the discarded binary labels, not the scores
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
        n_neighbors=lof_k,
        contamination="auto",  # contamination only affects the discarded binary labels, not the scores
        n_jobs=-1
    )
    lof.fit_predict(X)
    # negative_outlier_factor_ gives large negative numbers for outliers -> invert
    raw_lof = -lof.negative_outlier_factor_
    score_lof = to_percentile_score(raw_lof)
    
    # --- Algorithm 3: PCA Reconstruction Error ---
    # Choose the smallest number of components explaining >= 90% of variance,
    # capped at d-1 so a residual always exists.
    pca_full = PCA().fit(X)
    cum = np.cumsum(pca_full.explained_variance_ratio_)
    n_comp = int(min(max(1, np.searchsorted(cum, 0.90) + 1), X.shape[1] - 1))
    pca = PCA(n_components=n_comp).fit(X)
    if verbose:
        print(f"[PCA] explained variance ratio per component: {np.round(pca_full.explained_variance_ratio_, 3)}; "
              f"using {n_comp} component(s) (cumulative {cum[n_comp-1]:.3f})")
    X_recon = pca.inverse_transform(pca.transform(X))
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
    results.attrs["pca_n_components"] = n_comp
    results.attrs["pca_explained_variance_ratio"] = pca_full.explained_variance_ratio_.tolist()
    results.attrs["X_scaled"] = X
    
    return results

def assign_outlier_label(df: pd.DataFrame, X_scaled: np.ndarray, fragility: np.ndarray) -> pd.DataFrame:
    """
    Descriptive labels only. 'region' says where the object sits on the BPT plane;
    'driving_feature' says which scaled feature deviates most from the sample
    median (|robust z|). Neither asserts a physical mechanism; that requires
    follow-up data ([S II]/[O I] ratios, kinematics, imaging, multiwavelength).
    """
    x = df["log_nii_ha"].values
    y = df["log_oiii_hb"].values
    cls = df["bpt_class"].values

    region = np.full(len(df), "Outlier at ordinary BPT position", dtype=object)
    region[(cls == "AGN") & (y >= 0.9)] = "AGN region, very high [O III]/Hbeta (y >= 0.9)"
    region[(cls == "AGN") & (x >= 0.35) & (y < 0.9)] = "AGN region, high [N II]/Halpha (x >= 0.35)"
    region[(cls == "AGN") & (y <= -0.45) & (x >= 0.0)] = "AGN region, low [O III]/Hbeta (LINER-like ratios)"
    region[(x <= -0.7) & (y >= 0.7)] = "Low [N II]/Halpha, high [O III]/Hbeta (Green-Pea-like ratios; unconfirmed)"
    region[(x <= -1.1)] = "Very low [N II]/Halpha (metal-poor-like ratios; unconfirmed)"
    region[(cls == "Composite") & (region == "Outlier at ordinary BPT position")] = "Composite region"
    region[(cls == "Star-Forming") & (region == "Outlier at ordinary BPT position")] = "Star-forming region"
    region[fragility >= 0.75] = "Suspected measurement artifact (high fragility)"

    z = np.abs(X_scaled)  # RobustScaler output: (value - median) / IQR
    drv_idx = np.argmax(z, axis=1)
    driving = np.array(FEATURE_NAMES)[drv_idx]
    driving_z = np.round(z[np.arange(len(df)), drv_idx], 2)

    return pd.DataFrame({
        "bpt_region": region,
        "driving_feature": driving,
        "driving_feature_robust_z": driving_z
    }, index=df.index)

def stability_analysis(df: pd.DataFrame, seeds=(0, 1, 2, 3, 4), ks=(20, 35, 50), top_n=50,
                       out_csv="data/anomaly_stability.csv") -> pd.DataFrame:
    """
    Re-runs the ensemble for every (seed, k) combination and reports, per object,
    how many of the len(seeds)*len(ks) runs place it in the top `top_n` by consensus.
    Also prints the mean pairwise Spearman correlation of consensus scores and the
    mean pairwise Jaccard overlap of the top-`top_n` sets.
    Column 'n_runs_in_top' is used by run_bpt to require stability.
    """
    from scipy.stats import spearmanr
    from itertools import combinations

    runs = {}
    for seed in seeds:
        for k in ks:
            r = detect_anomalies_tri_ensemble(df, random_state=seed, lof_k=k, verbose=False)
            runs[(seed, k)] = r["consensus_anomaly_score"].values
    keys = list(runs)
    rhos, jacs = [], []
    tops = {key: set(np.argsort(-runs[key])[:top_n]) for key in keys}
    for a, b in combinations(keys, 2):
        rhos.append(spearmanr(runs[a], runs[b]).correlation)
        jacs.append(len(tops[a] & tops[b]) / len(tops[a] | tops[b]))
    n_in_top = np.zeros(len(df), dtype=int)
    for key in keys:
        n_in_top[list(tops[key])] += 1
    print(f"[Stability] {len(keys)} runs; mean Spearman(consensus) = {np.mean(rhos):.3f}; "
          f"mean Jaccard(top-{top_n}) = {np.mean(jacs):.3f}")
    res = pd.DataFrame({"n_runs_in_top": n_in_top, "n_runs_total": len(keys)}, index=df.index)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    res.to_csv(out_csv)
    return res

