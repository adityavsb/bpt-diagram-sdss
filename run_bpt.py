"""
Master Pipeline Runner: BPT Emission-Line Diagnostic & Anomaly Discovery (SDSS DR17).
Executes:
1. Catalog query with signal-to-noise quality filtering
2. Deterministic BPT classification (Kauffmann 2003, Kewley 2001, Kewley 2006)
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

from query_data import query_sdss_catalog, ANCHOR_GALAXIES
from bpt_math import classify_bpt
from plot_bpt import plot_bpt_diagram
from glitch_filter import (
    compute_line_snr,
    compute_error_propagation,
    compute_glitch_and_confidence_scores,
)
from anomaly_detector import (
    detect_anomalies_tri_ensemble,
    assign_physical_archetype,
)
from plot_anomalies import plot_bpt_anomalies

def run_pipeline(
    sample_size=10000,
    classic_plot="assets/bpt_diagnostic_sdss.png",
    anomaly_plot="assets/bpt_anomalies_sdss.png",
    candidates_csv="data/sdss_anomaly_candidates.csv",
    top_n_candidates=20
):
    print("=" * 85)
    print("SDSS DR17 BPT DIAGNOSTIC & TRI-ALGORITHMIC CELESTIAL ANOMALY DISCOVERY PIPELINE")
    print("=" * 85)

    # 1. Fetch data
    df = query_sdss_catalog(sample_size=sample_size)
    print(f"\n[1] Loaded {len(df):,} galaxies meeting S/N >= 3 cut on all 4 lines.")

    # 2. Deterministic BPT Classification
    log_nii_ha = df["log_nii_ha"].values
    log_oiii_hb = df["log_oiii_hb"].values
    df["bpt_class"] = classify_bpt(log_nii_ha, log_oiii_hb)

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
        "Seyfert (AGN)": "Supermassive Black Hole Accretion",
        "LINER": "Low-Ionization Shocks / Old Stars"
    }
    for c in ["Star-Forming", "Composite", "Seyfert (AGN)", "LINER"]:
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
    
    # 7. Error Propagation & Glitch Disentanglement
    print("[5] Propagating line flux errors & calculating 95% Confidence Intervals...")
    snr_df = compute_line_snr(df)
    err_df = compute_error_propagation(df)
    glitch_df = compute_glitch_and_confidence_scores(
        df,
        anomaly_df["consensus_anomaly_score"].values,
        snr_df,
        err_df
    )

    # Merge results
    df_merged = pd.concat([df, anomaly_df, snr_df, err_df, glitch_df], axis=1)

    # Assign physical archetype
    df_merged["physical_archetype"] = assign_physical_archetype(
        df_merged["log_nii_ha"].values,
        df_merged["log_oiii_hb"].values,
        df_merged["pipeline_glitch_risk"].values,
        df_merged["physical_anomaly_confidence"].values
    )

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
    # Top Physical: high physical_anomaly_confidence, low glitch risk (<= 50)
    physical_candidates = (
        df_merged[df_merged["pipeline_glitch_risk"] <= 50.0]
        .sort_values(by="physical_anomaly_confidence", ascending=False)
        .head(top_n_candidates)
        .copy()
    )
    physical_candidates["rank"] = range(1, len(physical_candidates) + 1)

    # Control Glitches: high glitch risk (>= 60), high raw anomaly score (look anomalous purely due to noise)
    glitch_candidates = (
        df_merged[df_merged["pipeline_glitch_risk"] >= 65.0]
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
        "physical_anomaly_confidence", "pipeline_glitch_risk",
        "snr_min", "physical_archetype", "skyserver_url", "legacy_survey_url"
    ]
    os.makedirs(os.path.dirname(candidates_csv), exist_ok=True)
    physical_candidates[output_cols].to_csv(candidates_csv, index=False)
    print(f"\n[6] Saved Top {top_n_candidates} Physical Anomaly Candidates catalog -> {candidates_csv}")

    # 8. Print Top 20 Candidates Table
    print("\n" + "=" * 105)
    print(f"3. TOP {top_n_candidates} CONFIRMED PHYSICAL CELESTIAL ANOMALIES (SDSS DR17)")
    print("=" * 105)
    header = f"{'Rank':<4} | {'SpecObjID':<19} | {'z':<6} | {'log([NII]/Hα)':<13} | {'log([OIII]/Hβ)':<14} | {'Votes':<5} | {'Score':<6} | {'PhysConf':<8} | {'Glitch%':<7} | {'Archetype':<26}"
    print(header)
    print("-" * 105)
    for _, row in physical_candidates.iterrows():
        print(
            f"#{row['rank']:<3} | {row['specObjID']:<19} | {row['z']:<6.3f} | "
            f"{row['log_nii_ha']:<13.2f} | {row['log_oiii_hb']:<14.2f} | "
            f"{row['algorithm_votes']}/3   | {row['consensus_anomaly_score']:<6.1f} | "
            f"{row['physical_anomaly_confidence']:<8.1f} | {row['pipeline_glitch_risk']:<7.1f} | "
            f"{row['physical_archetype'][:26]:<26}"
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
    args = parser.parse_args()
    run_pipeline(sample_size=args.sample_size)
