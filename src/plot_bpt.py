"""
Publication-quality BPT Diagnostic Visualization Generator.
Generates:
1. 2D Hexbin density distribution of SDSS galaxies with demarcation curves
2. Anchor galaxy callouts: M51 (Star-Forming) and NGC 1068 (Active Black Hole)
3. Classification population breakdown
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from bpt_math import kauffmann_boundary, kewley_boundary, classify_bpt
from query_data import ANCHOR_GALAXIES

def plot_bpt_diagram(df, output_path="assets/bpt_diagnostic_sdss.png"):
    log_nii_ha = df["log_nii_ha"].values
    log_oiii_hb = df["log_oiii_hb"].values
    if "bpt_class" not in df.columns:
        df = df.assign(bpt_class=classify_bpt(log_nii_ha, log_oiii_hb))

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), gridspec_kw={"width_ratios": [1.4, 1.0]})

    # -------------------------------------------------------------
    # Panel 1: The BPT Diagnostic Diagram (Density Hexbin + Curves)
    # -------------------------------------------------------------
    ax1 = axes[0]
    
    # 2D Hexbin density
    hb = ax1.hexbin(
        log_nii_ha, log_oiii_hb,
        gridsize=65, cmap="magma",
        bins="log", mincnt=1,
        extent=[-2.3, 0.8, -1.5, 1.6]
    )
    cb = fig.colorbar(hb, ax=ax1, orientation="vertical", pad=0.02)
    cb.set_label(r"$\log_{10}(\mathrm{Galaxy\ Count\ per\ Bin})$", fontsize=11)

    # 1. Kauffmann (2003) Star-Formation Boundary
    x_kauff = np.linspace(-2.3, 0.0, 300)
    y_kauff = kauffmann_boundary(x_kauff)
    ax1.plot(x_kauff, y_kauff, "k--", lw=2.2, label="Kauffmann et al. (2003) Pure Starburst Limit")

    # 2. Kewley (2001) Theoretical Extreme Starburst Upper Boundary
    x_kewley = np.linspace(-2.3, 0.45, 300)
    y_kewley = kewley_boundary(x_kewley)
    ax1.plot(x_kewley, y_kewley, "k-", lw=2.5, label="Kewley et al. (2001) Extreme Starburst Limit")

    # Region Labels
    ax1.text(-1.15, -0.20, "Star-Forming\n(OB Stars)", fontsize=11, fontweight="bold", color="#A3E4D7", ha="center")
    ax1.text(-0.22, -0.15, "Composite\n(Star + AGN)", fontsize=9.5, fontweight="bold", color="#FAD7A0", ha="center")
    ax1.text(0.18, 0.55, "AGN\n(Seyfert & LINER)", fontsize=11, fontweight="bold", color="#F1948A", ha="center")

    # Anchor Galaxies: M51 and NGC 1068
    for name, props in ANCHOR_GALAXIES.items():
        x = props["log_nii_ha"]
        y = props["log_oiii_hb"]
        ax1.scatter(
            x, y,
            s=220, color=props["color"],
            marker=props["marker"], edgecolor="white",
            linewidth=1.8, zorder=10,
            label=f"{name} ({props['type']})"
        )
        # Custom clean annotation placement to prevent edge clipping or collisions
        if "M51" in name:
            dx, dy, ha = -0.06, -0.15, "right"
        else:  # NGC 1068
            dx, dy, ha = -0.08, 0.14, "right"

        ax1.annotate(
            name, (x, y),
            xytext=(x + dx, y + dy),
            arrowprops=dict(arrowstyle="->", color=props["color"], lw=1.5),
            fontsize=10, fontweight="bold", color="white",
            bbox=dict(boxstyle="round,pad=0.25", facecolor=props["color"], alpha=0.9),
            ha=ha
        )

    ax1.set_xlim(-2.3, 0.8)
    ax1.set_ylim(-1.5, 1.6)
    ax1.set_xlabel(r"$\log_{10}([\mathrm{N\ II}]\,\lambda 6584 \,/\, \mathrm{H}\alpha)$", fontsize=12)
    ax1.set_ylabel(r"$\log_{10}([\mathrm{O\ III}]\,\lambda 5007 \,/\, \mathrm{H}\beta)$", fontsize=12)
    ax1.set_title(r"$\mathbf{BPT\ Emission\text{-}Line\ Diagnostic\ (SDSS\ DR17)}$", fontsize=13)
    ax1.legend(loc="lower left", fontsize=8.5, framealpha=0.9)
    ax1.grid(True, alpha=0.2, linestyle=":")

    # -------------------------------------------------------------
    # Panel 2: Physical Classification Breakdown
    # -------------------------------------------------------------
    ax2 = axes[1]
    counts = df["bpt_class"].value_counts()
    percentages = (counts / len(df)) * 100.0

    classes = ["Star-Forming", "Composite", "AGN"]
    colors = ["#2ECC71", "#F39C12", "#E74C3C"]
    y_pos = np.arange(len(classes))
    vals = [percentages.get(c, 0.0) for c in classes]
    raw_counts = [counts.get(c, 0) for c in classes]

    bars = ax2.barh(y_pos, vals, color=colors, height=0.55, edgecolor="black", alpha=0.85)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(classes, fontsize=11, fontweight="bold")
    ax2.invert_yaxis()
    ax2.set_xlabel("Percentage of Galaxy Sample (%)", fontsize=11)
    ax2.set_xlim(0, max(vals) * 1.25)
    ax2.set_title(f"Sample Classification Demographics\n({len(df):,} SDSS DR17 Galaxies, S/N ≥ 3)", fontsize=12, fontweight="bold")
    ax2.grid(True, alpha=0.25, axis="x", linestyle=":")

    for bar, val, raw in zip(bars, vals, raw_counts):
        w = bar.get_width()
        ax2.text(w + 1.2, bar.get_y() + bar.get_height() / 2, f"{val:.1f}% ({raw:,} galaxies)", va="center", fontsize=9.5, fontweight="bold")

    # Explanatory text box
    info_text = (
        "Core Physical Takeaways:\n"
        "• Star-Forming (OB Stars): Photoionized by short-lived massive stars (< 10 Myr).\n"
        "• Composite: Transition systems with mixed starburst and accretion energetics.\n"
        "• AGN (Active Nuclei): Hard power-law ionization from supermassive black holes.\n"
        "• Standard Note: Seyfert vs. LINER separation requires [S II]/Hα or [O I]/Hα."
    )
    ax2.text(
        0.05, 0.06, info_text,
        transform=ax2.transAxes,
        fontsize=9, family="sans-serif",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#F8F9FA", edgecolor="#BDC3C7", alpha=0.95)
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    print(f"\n[✓] Publication-quality BPT figure saved to: {os.path.abspath(output_path)}")

if __name__ == "__main__":
    from query_data import query_sdss_catalog
    df = query_sdss_catalog(sample_size=5000)
    plot_bpt_diagram(df)
