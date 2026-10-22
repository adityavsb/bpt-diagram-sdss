"""
Publication-quality visualization for SDSS DR17 BPT Anomaly Detection.
Overlays:
- 2D hexbin density of the bulk galaxy population
- Kauffmann (2003), Kewley (2001), and Kewley (2006) demarcation curves
- Top 20 Confirmed Physical Anomalies with 95% Confidence Interval error ellipses
- Control Pipeline Glitches with large measurement jitter error bars
- Real anchor galaxies: M51 & NGC 1068
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from bpt_math import (
    kauffmann_boundary,
    kewley_boundary,
    seyfert_liner_boundary,
)
from query_data import ANCHOR_GALAXIES

def plot_bpt_anomalies(
    df: pd.DataFrame,
    top_physical: pd.DataFrame,
    top_glitches: pd.DataFrame,
    output_path: str = "assets/bpt_anomalies_sdss.png"
):
    """
    Renders publication-ready figure of SDSS DR17 galaxy anomalies on BPT plane.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Dark modern astronomical palette
    fig, ax = plt.subplots(figsize=(12, 10), facecolor="#0B0F19")
    ax.set_facecolor("#0F1424")

    # 1. Background Hexbin Density
    x = df["log_nii_ha"].values
    y = df["log_oiii_hb"].values
    
    hb = ax.hexbin(
        x, y,
        gridsize=85,
        cmap="viridis",
        bins="log",
        mincnt=1,
        alpha=0.80,
        edgecolors="none",
        zorder=2
    )
    cb = fig.colorbar(hb, ax=ax, orientation="vertical", pad=0.02, shrink=0.82)
    cb.set_label(r"Galaxy Density $\log_{10}(N)$", color="#E2E8F0", fontsize=11, labelpad=10)
    cb.ax.tick_params(colors="#A0AEC0", labelsize=9)
    cb.outline.set_edgecolor("#2D3748")

    # 2. Demarcation Boundary Lines
    x_k03 = np.linspace(-1.5, 0.0, 300)
    y_k03 = kauffmann_boundary(x_k03)
    ax.plot(
        x_k03, y_k03,
        color="#38BDF8", linestyle="--", linewidth=2.0,
        label=r"Kauffmann+03 (Pure Starburst Boundary)",
        zorder=3
    )

    x_k01 = np.linspace(-1.5, 0.42, 300)
    y_k01 = kewley_boundary(x_k01)
    ax.plot(
        x_k01, y_k01,
        color="#F43F5E", linestyle="-", linewidth=2.2,
        label=r"Kewley+01 (Theoretical Starburst Upper Limit)",
        zorder=3
    )

    x_liner = np.linspace(-0.25, 0.5, 200)
    y_liner = seyfert_liner_boundary(x_liner)
    ax.plot(
        x_liner, y_liner,
        color="#F59E0B", linestyle=":", linewidth=2.0,
        label=r"Kewley+06 (Seyfert / LINER Demarcation)",
        zorder=3
    )

    # 3. Control Pipeline Glitches (Demonstrate measurement error jitter)
    if not top_glitches.empty:
        gx = top_glitches["log_nii_ha"].values
        gy = top_glitches["log_oiii_hb"].values
        gdx = top_glitches["delta_x_95"].values
        gdy = top_glitches["delta_y_95"].values
        
        # Plot large error bars for pipeline glitches
        ax.errorbar(
            gx, gy, xerr=gdx, yerr=gdy,
            fmt="x", color="#EF4444", ecolor="#EF4444", elinewidth=1.2,
            capsize=3, capthick=1.0, markersize=8, markeredgewidth=2.0,
            alpha=0.75, label="Pipeline Glitches (Large 95% CI Error Jitter)",
            zorder=4
        )

    # 4. Top Physical Celestial Anomalies with 95% Confidence Ellipses
    px = top_physical["log_nii_ha"].values
    py = top_physical["log_oiii_hb"].values
    pdx = top_physical["delta_x_95"].values
    pdy = top_physical["delta_y_95"].values
    ranks = top_physical["rank"].values

    for i in range(len(top_physical)):
        # 95% Confidence Ellipse (width = 2 * delta_x, height = 2 * delta_y)
        ell = Ellipse(
            xy=(px[i], py[i]),
            width=2 * pdx[i],
            height=2 * pdy[i],
            angle=0,
            edgecolor="#FDE047",
            facecolor="#FDE047",
            alpha=0.18,
            linewidth=1.2,
            zorder=5
        )
        ax.add_patch(ell)
        
        # Number annotation for Top 10
        if ranks[i] <= 10:
            ax.annotate(
                f"#{ranks[i]}",
                xy=(px[i], py[i]),
                xytext=(px[i] + 0.04, py[i] + 0.04),
                fontsize=8,
                fontweight="bold",
                color="#FEF08A",
                zorder=7,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="#1E293B", edgecolor="#FDE047", alpha=0.85)
            )

    # Scatter points for physical anomalies
    ax.scatter(
        px, py,
        marker="*", s=140, facecolors="#FFD700", edgecolors="#FFFFFF",
        linewidths=0.9, label=r"Top 20 Physical Anomalies (Tight 95% CI)",
        zorder=6
    )

    # 5. Anchor Galaxies: M51 and NGC 1068
    for name, props in ANCHOR_GALAXIES.items():
        ax.scatter(
            props["log_nii_ha"], props["log_oiii_hb"],
            color=props["color"], marker=props["marker"],
            s=220, edgecolors="white", linewidths=1.8,
            label=f"{name} ({props['type']})",
            zorder=8
        )
        ax.annotate(
            name.split("(")[0].strip(),
            xy=(props["log_nii_ha"], props["log_oiii_hb"]),
            xytext=(props["log_nii_ha"] + 0.05, props["log_oiii_hb"] - 0.08),
            fontsize=10, fontweight="bold", color=props["color"],
            zorder=9,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#0B0F19", edgecolor=props["color"], alpha=0.9)
        )

    # 6. Physical Regime Labels
    ax.text(-1.10, -0.65, "STARBURST LOCUS\n(OB Stellar Nurseries)",
            color="#38BDF8", fontsize=10, fontweight="bold", alpha=0.85, ha="center")
    ax.text(0.20, 1.05, "SEYFERT AGN WING\n(Black Hole Accretion)",
            color="#F43F5E", fontsize=10, fontweight="bold", alpha=0.85, ha="center")
    ax.text(0.28, -0.15, "LINERs\n(Shocks / Old Stars)",
            color="#F59E0B", fontsize=9, fontweight="bold", alpha=0.85, ha="center")
    ax.text(-0.15, 0.05, "COMPOSITE\n(Star + AGN)",
            color="#CBD5E1", fontsize=9, style="italic", alpha=0.75, ha="center")
    
    # Anomaly discovery arrows/tags
    ax.text(-1.15, 0.95, "▲ GREEN PEA / EELG\nCANDIDATES",
            color="#FDE047", fontsize=9, fontweight="bold", ha="left")

    # 7. Axis formatting & bounds
    ax.set_xlim(-1.6, 0.6)
    ax.set_ylim(-1.2, 1.4)
    ax.set_xlabel(r"$\log_{10}([\mathrm{N\ II}]\,\lambda 6584 \,/\, \mathrm{H}\alpha)$",
                  color="#F1F5F9", fontsize=13, labelpad=10)
    ax.set_ylabel(r"$\log_{10}([\mathrm{O\ III}]\,\lambda 5007 \,/\, \mathrm{H}\beta)$",
                  color="#F1F5F9", fontsize=13, labelpad=10)
    
    ax.set_title(
        "SDSS DR17 Galaxy Diagnostics: Tri-Algorithmic Celestial Anomaly Detection\n"
        r"Disentangling Physical Discoveries vs. Pipeline Glitches via 95% Confidence Intervals",
        color="#F8FAFC", fontsize=14, fontweight="bold", pad=16
    )

    ax.tick_params(colors="#94A3B8", labelsize=10)
    for spine in ax.spines.values():
        spine.set_color("#334155")
    ax.grid(True, linestyle=":", color="#334155", alpha=0.45)

    # Legend
    legend = ax.legend(
        loc="lower left", facecolor="#0B0F19", edgecolor="#334155",
        fontsize=9, labelcolor="#F1F5F9", framealpha=0.92
    )
    legend.get_frame().set_linewidth(1.0)

    # 8. Methodology Annotation Box
    info_text = (
        r"$\bf{Tri-Algorithmic\ Ensemble}$" + "\n"
        "• Isolation Forest + LOF + PCA Rec\n"
        "• Equal-Weighted Consensus Voting\n"
        r"$\bf{Glitch\ Disentanglement}$" + "\n"
        "• Analytical Gaussian Error Propagation\n"
        "• 95% CI Ellipses (1.96σ)\n"
        "• S/N Thresholding & Locus Overlap"
    )
    ax.text(
        0.98, 0.03, info_text,
        transform=ax.transAxes,
        fontsize=8.5,
        color="#E2E8F0",
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#111827", edgecolor="#374151", alpha=0.92)
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close()
    print(f"[Plot] Saved anomaly diagnostic diagram to: {output_path}")
