"""
Data query and caching module for SDSS DR17 emission-line galaxy sample.
Uses astroquery.sdss to retrieve official MPA-JHU emission line fluxes,
enforces signal-to-noise quality cuts, and caches locally.
"""
import os
import pandas as pd
import numpy as np

CACHE_FILE = os.path.join(os.path.dirname(__file__), "..", "data", "sdss_bpt_sample.csv")

# Anchors from open_astronomy/research/objects.md
ANCHOR_GALAXIES = {
    "M51 (Whirlpool Galaxy)": {
        "type": "Star-Forming Spiral",
        "ra": 202.469575,
        "dec": 47.195258,
        "z": 0.002,
        "log_nii_ha": -0.42,
        "log_oiii_hb": -0.38,
        "color": "#00A86B", # Emerald green
        "marker": "*"
    },
    "NGC 1068 (M77 Seyfert 2)": {
        "type": "Accreting Supermassive Black Hole",
        "ra": 40.66963,
        "dec": -0.01328,
        "z": 0.0038,
        "log_nii_ha": 0.15,
        "log_oiii_hb": 0.95,
        "color": "#D32F2F", # Crimson red
        "marker": "D"
    }
}

def query_sdss_catalog(sample_size=10000, force_refresh=False):
    """
    Retrieves or loads cached SDSS DR17 galaxy emission lines:
    H-alpha (6563A), H-beta (4861A), [O III] (5007A), [N II] (6584A).
    Quality cuts:
    - z between 0.02 and 0.25 (to keep lines inside optical bandpass)
    - Signal-to-Noise >= 3 on all 4 diagnostic emission lines
    """
    if os.path.exists(CACHE_FILE) and not force_refresh:
        print(f"[Data] Loading cached SDSS catalog from: {CACHE_FILE}")
        return pd.read_csv(CACHE_FILE)

    print(f"[Data] Querying SDSS DR17 SkyServer for {sample_size} emission-line galaxies...")
    from astroquery.sdss import SDSS

    sql_query = f"""
    SELECT TOP {sample_size}
      s.specObjID, s.bestObjID, s.ra, s.dec, s.z,
      l.h_alpha_flux, l.h_alpha_flux_err,
      l.h_beta_flux, l.h_beta_flux_err,
      l.oiii_5007_flux, l.oiii_5007_flux_err,
      l.nii_6584_flux, l.nii_6584_flux_err
    FROM SpecObj s
    JOIN galSpecLine l ON s.specObjID = l.specObjID
    WHERE s.class = 'GALAXY'
      AND s.z BETWEEN 0.02 AND 0.25
      AND s.zWarning = 0
      AND l.h_alpha_flux > 0 AND (l.h_alpha_flux / l.h_alpha_flux_err) >= 3.0
      AND l.h_beta_flux > 0 AND (l.h_beta_flux / l.h_beta_flux_err) >= 3.0
      AND l.oiii_5007_flux > 0 AND (l.oiii_5007_flux / l.oiii_5007_flux_err) >= 3.0
      AND l.nii_6584_flux > 0 AND (l.nii_6584_flux / l.nii_6584_flux_err) >= 3.0
    """

    table = SDSS.query_sql(sql_query, timeout=60)
    if table is None or len(table) == 0:
        raise RuntimeError("SDSS query returned no results. Check connection or constraints.")

    df = table.to_pandas()
    
    # Calculate logarithmic flux ratios
    df["log_nii_ha"] = np.log10(df["nii_6584_flux"] / df["h_alpha_flux"])
    df["log_oiii_hb"] = np.log10(df["oiii_5007_flux"] / df["h_beta_flux"])

    # Clean extreme outliers outside physical bounds
    df = df[(df["log_nii_ha"] > -2.0) & (df["log_nii_ha"] < 0.8) &
            (df["log_oiii_hb"] > -1.5) & (df["log_oiii_hb"] < 1.6)].copy()

    os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
    df.to_csv(CACHE_FILE, index=False)
    print(f"[Data] Cached {len(df)} validated galaxies to: {CACHE_FILE}")
    return df

if __name__ == "__main__":
    df = query_sdss_catalog(sample_size=5000)
    print(df[["z", "log_nii_ha", "log_oiii_hb"]].describe())
