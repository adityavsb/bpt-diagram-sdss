"""
Master Pipeline Runner: BPT Emission-Line Diagnostic & Anomaly Discovery (SDSS DR17).
Executes:
1. Catalog query with signal-to-noise quality filtering
2. Deterministic BPT classification (Kauffmann 2003, Kewley 2001)
3. Anchor galaxy comparison (M51 vs. NGC 1068)
4. Tri-Algorithmic Anomaly Detection (Isolation Forest, LOF, PCA Reconstruction)
5. Physical Anomaly vs. Pipeline Glitch Disentanglement with 95% Confidence Intervals
6. Top 20 Candidate Catalog generation and SkyServer navigation links
7. Visual diagnostic figure generation (Classic BPT + Anomaly Discovery)
"""
import sys
import os
import argparse
import pandas as pd
import numpy as np

# Add src to module search path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from query_data import query_sdss_catalog, ANCHOR_GALAXIES, fetch_anchor_ratios
from bpt_math import classify_bpt
from plot_bpt import plot_bpt_diagram
from glitch_filter import (
    compute_line_snr,
    compute_error_propagation,
    compute_fragility_and_outlier_scores,
    compute_boundary_sigma_distance,
)
from anomaly_detector import (
    detect_anomalies_tri_ensemble,
    assign_outlier_label,
)
from plot_anomalies import plot_bpt_anomalies

# Candidate selection thresholds (reported verbatim in README)
MAX_FRAGILITY = 0.25      # measurement_fragility <= this
MIN_SNR = 10.0            # minimum rescaled S/N over the four lines
MIN_VOTES = 2             # algorithms placing the object in their top 5%
MIN_STABLE_RUNS_FRAC = 1.0  # fraction of stability runs with object in top-50

def run_pipeline(
    sample_size=10000,
    force_refresh=False,
    classic_plot="assets/bpt_diagnostic_sdss.png",
    anomaly_plot="assets/bpt_anomalies_sdss.png",
    candidates_csv="data/sdss_anomaly_candidates.csv",
    top_n_candidates=20
):
    fetch_anchor_ratios()
    print("=" * 85)
    print("SDSS DR17 BPT DIAGNOSTIC & TRI-ALGORITHMIC CELESTIAL ANOMALY DISCOVERY PIPELINE")
    print("=" * 85)

    # 1. Fetch data
    df = query_sdss_catalog(sample_size=sample_size, force_refresh=force_refresh)
    print(f"\n[1] Loaded {len(df):,} galaxies meeting S/N >= 3 (MPA-JHU rescaled errors) on all 4 lines.")

    # 2. Deterministic BPT Classification
    log_nii_ha = df["log_nii_ha"].values
    log_oiii_hb = df["log_oiii_hb"].values
    df["bpt_class"] = classify_bpt(log_nii_ha, log_oiii_hb)

    # 2b. Validate against the MPA-JHU reference classification
    from bpt_math import validate_against_mpa_jhu
    agreement, _ = validate_against_mpa_jhu(df)
    if agreement < 0.95:
        raise RuntimeError(f"Classifier disagrees with MPA-JHU bptclass on {1-agreement:.1%} of objects; investigate before publishing numbers.")

    # 3. Print Demographics Table
    print("\n" + "=" * 85)
    print("1. DETERMINISTIC SPECTROSCOPIC DEMOGRAPHICS (CLASSICAL BPT)")
    print("=" * 85)
    print(f"{'Classification':<20} | {'Count':<10} | {'Percentage':<12} | {'Dominant Physical Mechanism':<30}")
    print("-" * 85)
    counts = df["bpt_class"].value_counts()
    mechanisms = {
        "Star-Forming": "OB Stellar Photons (< 10 Myr)",
        "Composite": "Transition (Starburst + AGN)",
        "AGN": "Supermassive Black Hole Accretion / Shocks"
    }
    for c in ["Star-Forming", "Composite", "AGN"]:
        cnt = counts.get(c, 0)
        pct = (cnt / len(df)) * 100.0
        mech = mechanisms.get(c, "Unknown")
        print(f"{c:<20} | {cnt:<10,} | {pct:<11.1f}% | {mech:<30}")
    print("=" * 85)

    # 4. Anchor Galaxies
    print("\n" + "=" * 85)
    print("2. ANCHOR GALAXIES (PHYSICAL GROUNDING)")
    print("=" * 85)
    print(f"{'Galaxy':<25} | {'log([N II]/Hα)':<16} | {'log([O III]/Hβ)':<16} | {'Class':<15}")
    print("-" * 85)
    for name, props in ANCHOR_GALAXIES.items():
        print(f"{name:<25} | {props['log_nii_ha']:<16.2f} | {props['log_oiii_hb']:<16.2f} | {props['type']:<15}")
    print("=" * 85)

    # 5. Render Classical BPT Plot
    print(f"\n[3] Generating classical BPT diagnostic figure -> {classic_plot}")
    plot_bpt_diagram(df, output_path=classic_plot)

    # 6. Tri-Algorithmic Anomaly Detection
    print("\n[4] Running Tri-Algorithmic Anomaly Ensemble (Isolation Forest + LOF + PCA Rec)...")
    anomaly_df = detect_anomalies_tri_ensemble(df, random_state=42)
    
    # 7. Error Propagation & Fragility Disentanglement
    print("[5] Propagating line flux errors & calculating 95% Confidence Intervals...")
    snr_df = compute_line_snr(df)
    err_df = compute_error_propagation(df)
    bnd_df = compute_boundary_sigma_distance(df, err_df)
    fragility_df = compute_fragility_and_outlier_scores(
        df,
        anomaly_df["consensus_anomaly_score"].values,
        snr_df,
        err_df
    )

    # Merge results
    df_merged = pd.concat([df, anomaly_df, snr_df, err_df, bnd_df, fragility_df], axis=1)

    # Assign outlier labels
    label_df = assign_outlier_label(
        df_merged,
        anomaly_df.attrs["X_scaled"],
        df_merged["measurement_fragility"].values
    )
    df_merged = pd.concat([df_merged, label_df], axis=1)

    # Add interactive SkyServer and DESI Legacy Survey links
    df_merged["skyserver_url"] = [
        f"http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid={sid}"
        for sid in df_merged["specObjID"]
    ]
    df_merged["legacy_survey_url"] = [
        f"https://www.legacysurvey.org/viewer?ra={ra:.5f}&dec={dec:.5f}&layer=ls-dr10&zoom=15"
        for ra, dec in zip(df_merged["ra"], df_merged["dec"])
    ]

    # Rank Candidates
    stability_cond = (
        (df_merged["n_runs_in_top"] >= MIN_STABLE_RUNS_FRAC * df_merged["n_runs_total"])
        if "n_runs_in_top" in df_merged.columns
        else True
    )
    sel = (
        (df_merged["measurement_fragility"] <= MAX_FRAGILITY)
        & (df_merged["snr_min"] >= MIN_SNR)
        & (df_merged["algorithm_votes"] >= MIN_VOTES)
        & stability_cond
    )
    print(f"[Selection] {int(sel.sum()):,} objects pass fragility<= {MAX_FRAGILITY}, S/N>= {MIN_SNR}, "
          f"votes>= {MIN_VOTES}")
    physical_candidates = (
        df_merged[sel]
        .sort_values(by="quality_weighted_outlier_score", ascending=False)
        .head(top_n_candidates)
        .copy()
    )
    physical_candidates["rank"] = range(1, len(physical_candidates) + 1)

    # Control Glitches: high fragility (>= 0.65), high raw anomaly score (look anomalous purely due to noise)
    glitch_candidates = (
        df_merged[df_merged["measurement_fragility"] >= 0.65]
        .sort_values(by="consensus_anomaly_score", ascending=False)
        .head(10)
        .copy()
    )

    # Save Catalog
    output_cols = [
        "rank", "specObjID", "ra", "dec", "z", "bpt_class",
        "log_nii_ha", "log_oiii_hb",
        "ci_95_nii_ha_str", "ci_95_oiii_hb_str",
        "consensus_anomaly_score", "algorithm_votes",
        "quality_weighted_outlier_score", "measurement_fragility",
        "snr_min", "n_sigma_nearest_boundary", "class_robust_3sigma",
        "bpt_region", "driving_feature", "driving_feature_robust_z",
        "skyserver_url", "legacy_survey_url"
    ]
    os.makedirs(os.path.dirname(candidates_csv), exist_ok=True)
    physical_candidates[output_cols].to_csv(candidates_csv, index=False)
    print(f"\n[6] Saved Top {top_n_candidates} Physical Anomaly Candidates catalog -> {candidates_csv}")

    # 8. Print Top 20 Candidates Table
    print("\n" + "=" * 105)
    print(f"3. TOP {top_n_candidates} LINE-RATIO OUTLIER CANDIDATES (SDSS DR17)")
    print("=" * 105)
    header = f"{'Rank':<4} | {'SpecObjID':<19} | {'z':<6} | {'log([NII]/Hα)':<13} | {'log([OIII]/Hβ)':<14} | {'Votes':<5} | {'Score':<6} | {'QWScore':<8} | {'Fragil':<7} | {'Region':<26}"
    print(header)
    print("-" * 105)
    for _, row in physical_candidates.iterrows():
        print(
            f"#{row['rank']:<3} | {row['specObjID']:<19} | {row['z']:<6.3f} | "
            f"{row['log_nii_ha']:<13.2f} | {row['log_oiii_hb']:<14.2f} | "
            f"{row['algorithm_votes']}/3   | {row['consensus_anomaly_score']:<6.1f} | "
            f"{row['quality_weighted_outlier_score']:<8.1f} | {row['measurement_fragility']:<7.3f} | "
            f"{row['bpt_region'][:26]:<26}"
        )
    print("=" * 105)

    # 9. Render Anomaly Plot
    print(f"\n[7] Generating publication-grade Anomaly Diagnostic Figure -> {anomaly_plot}...")
    plot_bpt_anomalies(
        df=df_merged,
        top_physical=physical_candidates,
        top_glitches=glitch_candidates,
        output_path=anomaly_plot
    )
    print("\n>>> Pipeline execution completed successfully. All deliverables ready.")
    return df_merged, physical_candidates

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SDSS DR17 BPT & Anomaly Pipeline")
    parser.add_argument("--sample-size", type=int, default=10000, help="Number of galaxies to evaluate")
    parser.add_argument("--force-refresh", action="store_true", help="Ignore local cache and re-query SkyServer")
    args = parser.parse_args()
    run_pipeline(sample_size=args.sample_size, force_refresh=args.force_refresh)
