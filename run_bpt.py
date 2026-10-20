"""
Master Pipeline Runner: BPT Emission-Line Diagnostic (SDSS DR17).
Executes:
1. Catalog query with signal-to-noise quality filtering
2. Physical classification (Kauffmann 2003, Kewley 2001, Kewley 2006)
3. Anchor galaxy comparison (M51 vs. NGC 1068)
4. Figure generation in assets/bpt_diagnostic_sdss.png
"""
import sys
import os

# Add src to module search path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from query_data import query_sdss_catalog, ANCHOR_GALAXIES
from bpt_math import classify_bpt
from plot_bpt import plot_bpt_diagram

def main(sample_size=10000, output_plot="assets/bpt_diagnostic_sdss.png"):
    print("=" * 75)
    print("BPT DIAGNOSTIC PIPELINE: STARBURSTS VS. ACCRETING BLACK HOLES (SDSS DR17)")
    print("=" * 75)

    # 1. Fetch data
    df = query_sdss_catalog(sample_size=sample_size)
    print(f"\n[1] Loaded {len(df):,} galaxies meeting S/N >= 3 cut on all 4 lines.")

    # 2. Classify
    log_nii_ha = df["log_nii_ha"].values
    log_oiii_hb = df["log_oiii_hb"].values
    df["bpt_class"] = classify_bpt(log_nii_ha, log_oiii_hb)

    # 3. Print Demographics Table
    print("\n" + "=" * 75)
    print("SPECTROSCOPIC CLASSIFICATION DEMOGRAPHICS")
    print("=" * 75)
    print(f"{'Classification':<20} | {'Count':<10} | {'Percentage':<12} | {'Dominant Physical Mechanism':<30}")
    print("-" * 75)
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
    print("=" * 75)

    # 4. Print Anchor Galaxies
    print("\n" + "=" * 75)
    print("ANCHOR GALAXIES (PHYSICAL COMPARISON)")
    print("=" * 75)
    print(f"{'Galaxy':<25} | {'log([N II]/Hα)':<16} | {'log([O III]/Hβ)':<16} | {'Class':<15}")
    print("-" * 75)
    for name, props in ANCHOR_GALAXIES.items():
        print(f"{name:<25} | {props['log_nii_ha']:<16.2f} | {props['log_oiii_hb']:<16.2f} | {props['type']:<15}")
    print("=" * 75)

    # 5. Generate Visual Figure
    print("\n[5] Generating publication-quality diagnostic figure...")
    plot_bpt_diagram(df, output_path=output_plot)

if __name__ == "__main__":
    main()
