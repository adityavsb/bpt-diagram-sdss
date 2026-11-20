"""
Data access for the single-object deep dive on candidate #2 (LEDA 4572013).

Every remote call is cached under data/deep_dive/, so the notebook re-runs
offline once the cache exists. Delete a cache file (or pass refresh=True) to
re-fetch it. Nothing here interprets the data; that happens in the notebook.

Sources
  - SDSS SkyServer DR17 SQL (galSpecLine / galSpecExtra / galSpecInfo / SpecObj / PhotoObj)
  - SDSS SAS: the spectrum FITS file itself
  - VizieR (Harvard mirror, ASU interface): ALFALFA, AllWISE, GALEX GUVcat_AIS
  - Legacy Survey viewer: grz JPEG cutouts
  - SIMBAD TAP (Harvard mirror): titles of the papers SIMBAD associates with the object
  - Crossref: verifies the class-literature references (title + DOI) instead of trusting memory
"""
import html
import io
import os
import re
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "deep_dive"

# The target is fixed by the candidate catalogue (data/sdss_anomaly_candidates.csv, rank 2).
TARGET = {
    "name": "LEDA 4572013",
    "alt_names": ["HS 1248+3523", "AGC 227547", "SDSS J125116.57+350648.4"],
    "specObjID": 2237320683103217664,
    "bestObjID": 1237665025444741310,
    "plate": 1987,
    "mjd": 53765,
    "fiberID": 573,
    "run2d": 26,
    "ra": 192.81909,
    "dec": 35.113453,
    "z": 0.03337908,
    "agc": 227547,
}

SKYSERVER_SQL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
VIZIER_ASU = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
SIMBAD_TAP = "https://simbad.harvard.edu/simbad/sim-tap"
HEADERS = {"User-Agent": "bpt-diagram-sdss-deep-dive/0.1"}


def _path(name: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    return CACHE / name


def sql(query: str, name: str, refresh: bool = False, timeout: int = 180) -> pd.DataFrame:
    """Run a SkyServer SQL query, cache the CSV, return a DataFrame."""
    p = _path(f"{name}.csv")
    if p.exists() and not refresh:
        return pd.read_csv(p)
    r = requests.get(SKYSERVER_SQL, params={"format": "csv", "cmd": query},
                     headers=HEADERS, timeout=timeout)
    r.raise_for_status()
    text = r.text
    if text.lstrip().lower().startswith(("error", "<")) or "ERROR" in text[:200]:
        raise RuntimeError(f"SkyServer rejected query '{name}': {text[:300]}")
    df = pd.read_csv(io.StringIO(text), comment="#")
    df.to_csv(p, index=False)
    return df


# ---------------------------------------------------------------- the target
def target_tables(refresh: bool = False) -> dict:
    """All SDSS/MPA-JHU tables for the target, one row each."""
    sid = TARGET["specObjID"]
    out = {
        "lines": sql(f"SELECT * FROM galSpecLine WHERE specObjID={sid}", "target_galSpecLine", refresh),
        "extra": sql(f"SELECT * FROM galSpecExtra WHERE specObjID={sid}", "target_galSpecExtra", refresh),
        "info": sql(f"SELECT * FROM galSpecInfo WHERE specObjID={sid}", "target_galSpecInfo", refresh),
        "spec": sql(f"SELECT * FROM SpecObj WHERE specObjID={sid}", "target_SpecObj", refresh),
        "photo": sql(
            "SELECT objID, ra, dec, type, flags, "
            "modelMag_u, modelMag_g, modelMag_r, modelMag_i, modelMag_z, "
            "modelMagErr_u, modelMagErr_g, modelMagErr_r, modelMagErr_i, modelMagErr_z, "
            "cModelMag_r, psfMag_r, petroR50_r, petroR90_r, fracDeV_r, "
            "extinction_u, extinction_g, extinction_r, extinction_i, extinction_z "
            f"FROM PhotoObjAll WHERE objID={TARGET['bestObjID']}", "target_photo", refresh),
    }
    return out


def spectrum_path(refresh: bool = False) -> Path:
    """Download the SDSS 'lite' spectrum FITS (about 170 kB) if it is not cached."""
    t = TARGET
    p = _path(f"spec-{t['plate']}-{t['mjd']}-{t['fiberID']:04d}.fits")
    if p.exists() and not refresh:
        return p
    url = (f"https://data.sdss.org/sas/dr17/sdss/spectro/redux/{t['run2d']}/spectra/lite/"
           f"{t['plate']}/spec-{t['plate']}-{t['mjd']}-{t['fiberID']:04d}.fits")
    r = requests.get(url, headers=HEADERS, timeout=120)
    r.raise_for_status()
    if not r.content.startswith(b"SIMPLE"):
        raise RuntimeError("SAS did not return a FITS file")
    p.write_bytes(r.content)
    return p


def spec_neighbours(radius_deg: float = 0.6, dv_kms: float = 2000.0, refresh: bool = False) -> pd.DataFrame:
    """SDSS spectroscopic neighbours near the target and within +-dv_kms of its redshift."""
    import numpy as np
    t = TARGET
    dz = dv_kms / 299792.458 * (1 + t["z"])
    dra = radius_deg / np.cos(np.radians(t["dec"]))
    q = (
        "SELECT s.specObjID, s.ra, s.dec, s.z, s.class, s.subClass, "
        "p.modelMag_r, p.modelMag_g, p.petroR50_r "
        "FROM SpecObj AS s JOIN PhotoObjAll AS p ON p.objID = s.bestObjID "
        f"WHERE s.ra BETWEEN {t['ra']-dra} AND {t['ra']+dra} "
        f"AND s.dec BETWEEN {t['dec']-radius_deg} AND {t['dec']+radius_deg} "
        f"AND s.z BETWEEN {t['z']-dz} AND {t['z']+dz} AND s.sciencePrimary = 1"
    )
    return sql(q, "target_spec_neighbours", refresh)


# ---------------------------------------------------------------- controls
_CONTROL_FILTER = (
    "s.class = 'GALAXY' AND s.zWarning = 0 AND s.sciencePrimary = 1 AND i.reliable = 1 "
)
_CONTROL_FROM = (
    "FROM SpecObj AS s JOIN galSpecInfo AS i ON s.specObjID = i.specObjID "
    "JOIN galSpecExtra AS e ON s.specObjID = e.specObjID "
    "JOIN galSpecLine AS l ON s.specObjID = l.specObjID "
)


def mzr_control(refresh: bool = False) -> pd.DataFrame:
    """
    Star-forming galaxies with MPA-JHU stellar mass and gas-phase O/H.
    Two pieces: a general random draw and a dedicated low-mass slice (otherwise
    the dwarf regime is nearly empty in a random draw).
    """
    cols = ("s.specObjID, s.z, e.lgm_tot_p50, e.oh_p16, e.oh_p50, e.oh_p84, "
            "e.sfr_tot_p50, e.specsfr_tot_p50, l.h_alpha_eqw, l.h_alpha_flux, l.h_beta_flux, "
            "l.nii_6584_flux, l.oiii_5007_flux, l.oii_3726_flux, l.oii_3729_flux")
    common = (f"WHERE {_CONTROL_FILTER} AND e.bptclass = 1 AND e.oh_p50 > 0 "
              "AND e.lgm_tot_p50 > 0 AND s.z BETWEEN 0.01 AND 0.12 ")
    general = sql(f"SELECT TOP 30000 {cols} {_CONTROL_FROM} {common} ORDER BY ABS(CHECKSUM(s.specObjID))",
                  "control_mzr_general", refresh)
    low = sql(f"SELECT TOP 6000 {cols} {_CONTROL_FROM} {common} AND e.lgm_tot_p50 BETWEEN 7.3 AND 8.9 "
              "ORDER BY ABS(CHECKSUM(s.specObjID))", "control_mzr_lowmass", refresh)
    df = pd.concat([general.assign(slice="general"), low.assign(slice="lowmass")])
    return df.drop_duplicates("specObjID").reset_index(drop=True)


def te_control(refresh: bool = False) -> pd.DataFrame:
    """
    Star-forming galaxies with a detected [O III] 4363 line (formal S/N >= 5):
    the control set on which the direct-Te code is validated and the target is compared.
    Note: MPA-JHU gives no recommended error rescaling for 4363, so S/N here is formal.
    """
    cols = ("s.specObjID, s.z, i.e_bv_sfd, e.lgm_tot_p50, e.oh_p50, e.sfr_tot_p50, "
            "l.oii_3726_flux, l.oii_3726_flux_err, l.oii_3729_flux, l.oii_3729_flux_err, "
            "l.h_beta_flux, l.h_beta_flux_err, l.h_alpha_flux, l.h_alpha_flux_err, "
            "l.h_gamma_flux, l.h_gamma_flux_err, "
            "l.oiii_4363_flux, l.oiii_4363_flux_err, l.oiii_4959_flux, l.oiii_5007_flux, "
            "l.oiii_5007_flux_err, l.nii_6584_flux, l.nii_6584_flux_err, "
            "l.sii_6717_flux, l.sii_6731_flux, l.h_alpha_eqw")
    where = (f"WHERE {_CONTROL_FILTER} AND e.bptclass = 1 AND s.z BETWEEN 0.03 AND 0.30 "
             "AND l.oiii_4363_flux > 0 AND l.oiii_4363_flux_err > 0 "
             "AND l.oiii_4363_flux / l.oiii_4363_flux_err >= 5 "
             "AND l.h_beta_flux_err > 0 AND l.h_beta_flux / l.h_beta_flux_err >= 10 "
             "AND l.h_alpha_flux_err > 0 AND l.oii_3726_flux > 0 AND l.oii_3729_flux > 0 "
             "AND l.oiii_5007_flux > 0 AND l.oiii_4959_flux > 0 AND l.nii_6584_flux > 0 "
             "AND l.nii_6584_flux_err > 0 AND l.h_gamma_flux > 0")
    return sql(f"SELECT TOP 8000 {cols} {_CONTROL_FROM} {where} ORDER BY ABS(CHECKSUM(s.specObjID))",
               "control_te", refresh)


# ---------------------------------------------------------------- VizieR
def vizier(source: str, name: str, ra=None, dec=None, radius_arcsec=None, where: str = "",
           cols: str = "*", max_rows: int = 20, refresh: bool = False) -> pd.DataFrame:
    """Query one VizieR table through the ASU TSV interface (Harvard mirror). Empty frame if no rows."""
    p = _path(f"{name}.csv")
    if p.exists() and not refresh:
        try:
            return pd.read_csv(p)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    params = {"-source": source, "-out.max": max_rows, "-out": cols, "-oc.form": "d"}
    if ra is not None:
        params["-c"] = f"{ra:.6f} {dec:+.6f}"
        if radius_arcsec is not None:
            params["-c.rs"] = radius_arcsec
    q = "&".join(f"{k}={requests.utils.quote(str(v))}" for k, v in params.items())
    r = requests.get(f"{VIZIER_ASU}?{q}" + (f"&{where}" if where else ""), headers=HEADERS, timeout=90)
    r.raise_for_status()
    lines = [ln for ln in r.text.splitlines() if ln.strip() and not ln.startswith("#")]
    # layout: header, units, dashes, then data rows
    if len(lines) < 4:
        pd.DataFrame().to_csv(p, index=False)
        return pd.DataFrame()
    header = [h.strip() for h in lines[0].split("\t")]
    rows = [[c.strip() for c in ln.split("\t")] for ln in lines[3:]]
    df = pd.DataFrame(rows, columns=header).replace("", pd.NA)
    for c in df.columns:
        conv = pd.to_numeric(df[c], errors="coerce")
        if conv.notna().sum() == df[c].notna().sum():
            df[c] = conv
    df.to_csv(p, index=False)
    return df


def alfalfa(refresh: bool = False):
    """ALFALFA (Haynes+18) row for AGC 227547, and every ALFALFA source within 10 arcmin."""
    t = TARGET
    row = vizier("J/ApJ/861/49/table2", "alfalfa_target", where=f"AGC={t['agc']}", refresh=refresh)
    near = vizier("J/ApJ/861/49/table2", "alfalfa_within_10arcmin", ra=t["ra"], dec=t["dec"],
                  radius_arcsec=600, cols="AGC,Name,RAJ2000,DEJ2000,RAO,DEO,Vhel,W50,HIflux,SNR,Dist,logMHI,HI",
                  max_rows=50, refresh=refresh)
    return row, near


HI_BOX = {"ra_min": 175.0, "ra_max": 210.0, "dec_min": 20.0, "dec_max": 40.0}


def hi_control(refresh: bool = False):
    """
    H I control sample for the gas-fraction test: every ALFALFA source in a sky box around
    the target, plus MPA-JHU stellar masses of SDSS galaxies in the same box. The notebook
    matches them on the ALFALFA optical-counterpart position, so the control has the same
    ALFALFA depth and the same SDSS mass estimator as the target.
    """
    b = HI_BOX
    ra_c, dec_c = (b["ra_min"] + b["ra_max"]) / 2, (b["dec_min"] + b["dec_max"]) / 2
    w, h = b["ra_max"] - b["ra_min"], b["dec_max"] - b["dec_min"]
    # ASU box search is in the tangent plane; widen a little and clip in the notebook
    hi = vizier("J/ApJ/861/49/table2", "alfalfa_box", ra=ra_c, dec=dec_c, radius_arcsec=None,
                where=f"-c.bd={w:.1f}x{h:.1f}",
                cols="AGC,RAO,DEO,Vhel,W50,HIflux,SNR,Dist,logMHI,e_logMHI,HI",
                max_rows=20000, refresh=refresh)
    sdss = sql(
        "SELECT s.specObjID, s.ra, s.dec, s.z, e.lgm_tot_p50, e.bptclass, e.sfr_tot_p50 "
        "FROM SpecObj AS s JOIN galSpecInfo AS i ON s.specObjID = i.specObjID "
        "JOIN galSpecExtra AS e ON s.specObjID = e.specObjID "
        f"WHERE s.class = 'GALAXY' AND s.zWarning = 0 AND s.sciencePrimary = 1 AND i.reliable = 1 "
        f"AND s.z BETWEEN 0.001 AND 0.06 AND e.lgm_tot_p50 > 0 "
        f"AND s.ra BETWEEN {b['ra_min']} AND {b['ra_max']} AND s.dec BETWEEN {b['dec_min']} AND {b['dec_max']}",
        "sdss_masses_hi_box", refresh)
    return hi, sdss


def allwise(refresh: bool = False) -> pd.DataFrame:
    t = TARGET
    return vizier("II/328/allwise", "allwise", ra=t["ra"], dec=t["dec"], radius_arcsec=6,
                  cols="AllWISE,RAJ2000,DEJ2000,W1mag,e_W1mag,W2mag,e_W2mag,W3mag,e_W3mag,W4mag,e_W4mag,ccf,ext",
                  refresh=refresh)


def galex(refresh: bool = False) -> pd.DataFrame:
    t = TARGET
    return vizier("II/335/galex_ais", "galex_ais", ra=t["ra"], dec=t["dec"], radius_arcsec=15,
                  cols="Name,RAJ2000,DEJ2000,FUVmag,e_FUVmag,NUVmag,e_NUVmag,Nafl,Nexf", refresh=refresh)


# ---------------------------------------------------------------- images
def legacy_cutout(name: str, pixscale: float, size: int, refresh: bool = False) -> Path:
    t = TARGET
    p = _path(f"{name}.jpg")
    if p.exists() and not refresh:
        return p
    r = requests.get("https://www.legacysurvey.org/viewer/jpeg-cutout",
                     params={"ra": t["ra"], "dec": t["dec"], "layer": "ls-dr10",
                             "pixscale": pixscale, "size": size},
                     headers=HEADERS, timeout=90)
    r.raise_for_status()
    p.write_bytes(r.content)
    return p


# ---------------------------------------------------------------- literature
def simbad_references(bibcodes, refresh: bool = False) -> pd.DataFrame:
    """Titles / years / DOIs for the bibcodes SIMBAD lists for the object."""
    p = _path("simbad_references.csv")
    if p.exists() and not refresh:
        return pd.read_csv(p)
    import pyvo
    quoted = ",".join("'" + b.replace("'", "''") + "'" for b in bibcodes)
    svc = pyvo.dal.TAPService(SIMBAD_TAP)
    res = svc.search(f'SELECT bibcode, title, "year", journal, doi FROM ref WHERE bibcode IN ({quoted})')
    df = res.to_table().to_pandas()
    df.to_csv(p, index=False)
    return df


def crossref_lookup(queries: dict, refresh: bool = False) -> pd.DataFrame:
    """
    For each {key: free-text citation} return Crossref's best match (title, first author,
    year, DOI). The caller must read the result and confirm it is the intended paper;
    this is a lookup, not a verification by itself.
    """
    p = _path("crossref_lookup.csv")
    if p.exists() and not refresh:
        return pd.read_csv(p)
    rows = []
    for key, text in queries.items():
        for attempt in range(5):
            r = requests.get("https://api.crossref.org/works",
                             params={"query.bibliographic": text, "rows": 1,
                                     "select": "title,author,issued,DOI,container-title"},
                             headers=HEADERS, timeout=60)
            if r.status_code != 429:
                break
            time.sleep(float(r.headers.get("Retry-After", 5)) + 2 * attempt)
        r.raise_for_status()
        items = r.json()["message"]["items"]
        if not items:
            rows.append({"key": key, "query": text})
            continue
        it = items[0]

        def clean(s):  # Crossref titles can carry HTML tags and entities
            return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", s))).strip()

        rows.append({
            "key": key, "query": text,
            "title": clean((it.get("title") or [""])[0]),
            "first_author": (it.get("author") or [{}])[0].get("family", ""),
            "year": (it.get("issued", {}).get("date-parts") or [[None]])[0][0],
            "journal": clean((it.get("container-title") or [""])[0]),
            "doi": it.get("DOI", ""),
        })
        time.sleep(1.5)
    df = pd.DataFrame(rows)
    df.to_csv(p, index=False)
    return df


def fetch_all(refresh: bool = False) -> None:
    """Populate the whole cache (used once; the notebook only reads)."""
    target_tables(refresh)
    spectrum_path(refresh)
    spec_neighbours(refresh=refresh)
    mzr_control(refresh)
    te_control(refresh)
    alfalfa(refresh)
    hi_control(refresh)
    allwise(refresh)
    galex(refresh)
    legacy_cutout("legacy_120arcsec", 0.4, 300, refresh)
    legacy_cutout("legacy_8arcmin", 1.2, 400, refresh)


if __name__ == "__main__":
    fetch_all()
    print("cache:", sorted(p.name for p in CACHE.iterdir()))
