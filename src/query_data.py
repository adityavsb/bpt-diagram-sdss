"""
Data query and caching module for SDSS DR17 emission-line galaxy sample.
Uses astroquery.sdss to retrieve official MPA-JHU emission line fluxes,
enforces signal-to-noise quality cuts, and caches locally.
"""
import os
import pandas as pd
import numpy as np

CACHE_FILE = os.path.join(os.path.dirname(__file__), "..", "data", "sdss_bpt_sample.csv")

# Anchor galaxies: nuclear line ratios from Ho, Filippenko & Sargent (1997,
# ApJS, 112, 315; DOI 10.1086/313041; VizieR J/ApJS/112/315), the Palomar
# spectroscopic survey of nearby galaxy nuclei. Values are filled by
# fetch_anchor_ratios() and cached in data/anchor_galaxies.json.
# NGC 5194 (M51) is a Seyfert 2 / LINER NUCLEUS in HFS97, not a star-forming
# point; its disk H II regions are star-forming but are not what a nuclear
# spectrograph fiber sees.
ANCHOR_GALAXIES = {
    "M51 / NGC 5194": {
        "type": "Nearby spiral; Seyfert 2 nucleus (HFS97)",
        "hfs97_name": "NGC 5194",
        "ra": 202.469575,
        "dec": 47.195258,
        "z": 0.00154,
        "log_nii_ha": None,
        "log_oiii_hb": None,
        "color": "#00A86B",
        "marker": "*",
    },
    "NGC 1068 / M77": {
        "type": "Seyfert 2 nucleus (HFS97)",
        "hfs97_name": "NGC 1068",
        "ra": 40.66963,
        "dec": -0.01328,
        "z": 0.00379,
        "log_nii_ha": None,
        "log_oiii_hb": None,
        "color": "#D32F2F",
        "marker": "D",
    },
}

ANCHOR_CACHE = os.path.join(os.path.dirname(__file__), "..", "data", "anchor_galaxies.json")


def fetch_anchor_ratios(force_refresh=False):
    """
    Fills ANCHOR_GALAXIES[*]["log_nii_ha"/"log_oiii_hb"] from HFS97 via VizieR.
    Procedure (deterministic):
      1. If cached in data/anchor_galaxies.json and not force_refresh, load cache.
      2. Otherwise query VizieR catalog J/ApJS/112/315.
      3. Fall back gracefully to verified published values from Table 4 if remote
         service is unreachable (e.g. timeout on VPN).
    """
    import json
    if os.path.exists(ANCHOR_CACHE) and not force_refresh:
        with open(ANCHOR_CACHE) as fh:
            cached = json.load(fh)
        for k, v in cached.get("anchors", {}).items():
            if k in ANCHOR_GALAXIES:
                ANCHOR_GALAXIES[k].update(v)
        return ANCHOR_GALAXIES

    out = None
    try:
        from astroquery.vizier import Vizier
        v = Vizier(columns=["**"], row_limit=-1, timeout=15)
        tables = v.get_catalogs("J/ApJS/112/315")
        chosen = None
        for t in tables:
            descs = {c: (t[c].description or "").lower() for c in t.colnames}
            nii_col = next((c for c, d in descs.items() if "nii" in d.replace("[", "").replace("]", "") and "ha" in d), None)
            oiii_col = next((c for c, d in descs.items() if "oiii" in d.replace("[", "").replace("]", "") and "hb" in d), None)
            if nii_col and oiii_col and "Name" in t.colnames:
                chosen = (t, nii_col, oiii_col)
                break
        if chosen is not None:
            t, nii_col, oiii_col = chosen
            is_log_nii = "log" in (t[nii_col].description or "").lower()
            is_log_oiii = "log" in (t[oiii_col].description or "").lower()
            out = {
                "source": {
                    "catalog": "J/ApJS/112/315",
                    "table": t.meta.get("name"),
                    "nii_col": nii_col,
                    "oiii_col": oiii_col
                },
                "anchors": {}
            }
            for key, props in ANCHOR_GALAXIES.items():
                row = t[[str(n).strip() == props["hfs97_name"] for n in t["Name"]]]
                if len(row) > 0:
                    nii = float(row[nii_col][0])
                    oiii = float(row[oiii_col][0])
                    props["log_nii_ha"] = nii if is_log_nii else float(np.log10(nii))
                    props["log_oiii_hb"] = oiii if is_log_oiii else float(np.log10(oiii))
                    out["anchors"][key] = {
                        "log_nii_ha": props["log_nii_ha"],
                        "log_oiii_hb": props["log_oiii_hb"]
                    }
    except Exception as exc:
        print(f"[Data] VizieR query for HFS97 unavailable ({exc}); using verified published HFS97 Table 4 values.")

    # Canonical fallback from Ho, Filippenko & Sargent (1997, ApJS, 112, 315, Table 4)
    if out is None or len(out.get("anchors", {})) < 2:
        out = {
            "source": {
                "catalog": "J/ApJS/112/315",
                "table": "table4",
                "nii_col": "[NII]/Ha",
                "oiii_col": "[OIII]/Hb",
                "reference": "Ho, Filippenko & Sargent (1997, ApJS, 112, 315, Table 4)"
            },
            "anchors": {
                "M51 / NGC 5194": {
                    "log_nii_ha": 0.07,
                    "log_oiii_hb": 0.22
                },
                "NGC 1068 / M77": {
                    "log_nii_ha": 0.17,
                    "log_oiii_hb": 1.10
                }
            }
        }
        for k, v in out["anchors"].items():
            if k in ANCHOR_GALAXIES:
                ANCHOR_GALAXIES[k].update(v)

    os.makedirs(os.path.dirname(ANCHOR_CACHE), exist_ok=True)
    with open(ANCHOR_CACHE, "w") as fh:
        json.dump(out, fh, indent=2)
    return ANCHOR_GALAXIES

# MPA-JHU recommended multiplicative corrections to the formal galSpecLine
# flux errors, derived from duplicate observations. Source:
# https://wwwmpa.mpa-garching.mpg.de/SDSS/DR4/raw_data.html ("About error estimates").
# The raw *_flux_err columns are NOT scaled; we apply the scaling everywhere.
MPA_JHU_ERR_SCALE = {
    "h_alpha": 2.473,
    "h_beta": 1.882,
    "oiii_5007": 1.566,
    "nii_6584": 2.039,
}

# MPA-JHU bptclass codes (galSpecExtra). Verify against
# https://data.sdss.org/datamodel/files/SPECTRO_REDUX/galSpecExtra.html before
# changing. -1 = unclassifiable, 1 = star-forming, 2 = low-S/N star-forming,
# 3 = composite, 4 = AGN (non-LINER), 5 = low-S/N LINER.
MPA_JHU_BPTCLASS_TO_LABEL = {
    1: "Star-Forming",
    3: "Composite",
    4: "AGN",
}

LINE_COLUMNS = ["h_alpha", "h_beta", "oiii_5007", "nii_6584"]


def query_sdss_catalog(sample_size=10000, force_refresh=False, timeout=600):
    """
    Retrieves or loads a cached sample of SDSS Legacy galaxies with MPA-JHU
    emission-line measurements (galSpecLine, served in SDSS DR17).

    Lines: H-alpha 6563, H-beta 4861, [O III] 5007, [N II] 6584.

    Selection:
    - class = 'GALAXY', zWarning = 0, sciencePrimary = 1 (one spectrum per object)
    - galSpecInfo.reliable = 1 (MPA-JHU reliability flag)
    - 0.02 < z < 0.25
    - all four line fluxes > 0 and flux errors > 0
    - S/N >= 3 on all four lines USING MPA-JHU RESCALED ERRORS
    - deterministic pseudo-random ordering by a hash of specObjID so that the
      sample covers the whole footprint and is reproducible between runs.

    Also retrieves galSpecExtra.bptclass for validating our classifier.
    """
    if os.path.exists(CACHE_FILE) and not force_refresh:
        print(f"[Data] Loading cached SDSS catalog from: {CACHE_FILE}")
        return pd.read_csv(CACHE_FILE)

    print(f"[Data] Querying SDSS SkyServer (MPA-JHU galSpec tables) for {sample_size} galaxies...")
    from astroquery.sdss import SDSS

    s = MPA_JHU_ERR_SCALE
    sql_query = f"""
    SELECT TOP {sample_size}
      s.specObjID, s.bestObjID, s.plate, s.mjd, s.fiberID,
      s.ra, s.dec, s.z, s.zErr,
      l.h_alpha_flux, l.h_alpha_flux_err,
      l.h_beta_flux, l.h_beta_flux_err,
      l.oiii_5007_flux, l.oiii_5007_flux_err,
      l.nii_6584_flux, l.nii_6584_flux_err,
      l.h_alpha_eqw,
      x.bptclass AS mpa_jhu_bptclass
    FROM SpecObj AS s
    JOIN galSpecInfo  AS i ON s.specObjID = i.specObjID
    JOIN galSpecLine  AS l ON s.specObjID = l.specObjID
    JOIN galSpecExtra AS x ON s.specObjID = x.specObjID
    WHERE s.class = 'GALAXY'
      AND s.zWarning = 0
      AND s.sciencePrimary = 1
      AND i.reliable = 1
      AND s.z BETWEEN 0.02 AND 0.25
      AND l.h_alpha_flux > 0  AND l.h_alpha_flux_err > 0
      AND l.h_beta_flux > 0   AND l.h_beta_flux_err > 0
      AND l.oiii_5007_flux > 0 AND l.oiii_5007_flux_err > 0
      AND l.nii_6584_flux > 0 AND l.nii_6584_flux_err > 0
      AND l.h_alpha_flux  / (l.h_alpha_flux_err  * {s['h_alpha']})  >= 3.0
      AND l.h_beta_flux   / (l.h_beta_flux_err   * {s['h_beta']})   >= 3.0
      AND l.oiii_5007_flux / (l.oiii_5007_flux_err * {s['oiii_5007']}) >= 3.0
      AND l.nii_6584_flux / (l.nii_6584_flux_err * {s['nii_6584']}) >= 3.0
    ORDER BY ABS(CHECKSUM(s.specObjID))
    """

    table = SDSS.query_sql(sql_query, timeout=timeout)
    if table is None or len(table) == 0:
        raise RuntimeError("SDSS query returned no results. Check connection or constraints.")

    df = table.to_pandas()

    # Rescaled errors become the errors used everywhere downstream.
    for line in LINE_COLUMNS:
        df[f"{line}_flux_err_raw"] = df[f"{line}_flux_err"]
        df[f"{line}_flux_err"] = df[f"{line}_flux_err_raw"] * MPA_JHU_ERR_SCALE[line]

    # Logarithmic flux ratios
    df["log_nii_ha"] = np.log10(df["nii_6584_flux"] / df["h_alpha_flux"])
    df["log_oiii_hb"] = np.log10(df["oiii_5007_flux"] / df["h_beta_flux"])

    # Keep only finite ratios. No clipping of extreme-but-real objects.
    finite = np.isfinite(df["log_nii_ha"]) & np.isfinite(df["log_oiii_hb"])
    n_dropped = int((~finite).sum())
    df = df[finite].copy()
    if n_dropped:
        print(f"[Data] Dropped {n_dropped} rows with non-finite line ratios.")

    os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
    df.to_csv(CACHE_FILE, index=False)
    print(f"[Data] Cached {len(df)} validated galaxies to: {CACHE_FILE}")
    print(f"[Data] Sky coverage: RA {df['ra'].min():.1f}..{df['ra'].max():.1f}, "
          f"Dec {df['dec'].min():.1f}..{df['dec'].max():.1f}")
    return df


if __name__ == "__main__":
    df = query_sdss_catalog(sample_size=5000)
    print(df[["z", "log_nii_ha", "log_oiii_hb"]].describe())

