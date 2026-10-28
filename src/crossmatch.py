"""
SIMBAD cross-match for outlier candidates. Produces, per candidate, the main
identifier, SIMBAD object type, redshift, all identifiers, and the list of
bibcodes SIMBAD associates with the object. Only facts that appear in this
output may be used in README candidate write-ups.
"""
import time
import os
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad


def crossmatch_candidates(candidates: pd.DataFrame, radius_arcsec: float = 5.0,
                          out_csv: str = "data/sdss_anomaly_candidates_simbad.csv") -> pd.DataFrame:
    """
    Cross-matches candidates against the SIMBAD astronomical database.
    Queries the official mirror at Harvard for high reliability.
    """
    s = Simbad()
    s.server = "simbad.harvard.edu"
    s.TIMEOUT = 60
    s.add_votable_fields("otype", "rvz_redshift", "morphtype", "ids", "nbref")

    rows = []
    for _, c in candidates.iterrows():
        coord = SkyCoord(ra=c["ra"] * u.deg, dec=c["dec"] * u.deg, frame="icrs")
        try:
            res = s.query_region(coord, radius=radius_arcsec * u.arcsec)
        except Exception as e:
            print(f"[Crossmatch] Query warning for candidate {c.get('rank')}: {e}")
            res = None

        time.sleep(0.3)  # be polite to CDS / Harvard
        if res is None or len(res) == 0:
            rows.append({
                "rank": c["rank"],
                "specObjID": c["specObjID"],
                "simbad_main_id": None,
                "simbad_otype": None,
                "simbad_z": None,
                "simbad_morphtype": None,
                "simbad_ids": None,
                "simbad_n_bibcodes": 0,
                "simbad_bibcodes": "",
                "simbad_url": None,
            })
            continue

        r = res[0]
        main_id = str(r["main_id"]).strip()
        ids = str(r["ids"]).replace("|", "; ")
        nbref = int(r["nbref"]) if "nbref" in res.colnames and r["nbref"] is not None else 0
        morph = str(r["morph_type"]).strip() if "morph_type" in res.colnames and r["morph_type"] is not None else None
        if morph == "--" or morph == "":
            morph = None
        rvz = float(r["rvz_redshift"]) if "rvz_redshift" in res.colnames and r["rvz_redshift"] is not None else None

        bib_list = []
        if nbref > 0:
            try:
                q = (
                    "SELECT ref.bibcode FROM ref "
                    "JOIN has_ref ON ref.oidbib = has_ref.oidbibref "
                    "JOIN basic ON basic.oid = has_ref.oidref "
                    f"WHERE basic.main_id = '{main_id}'"
                )
                tap_res = s.tap.search(q).to_table()
                if tap_res is not None and "bibcode" in tap_res.colnames:
                    bib_list = [str(b).strip() for b in tap_res["bibcode"]]
            except Exception as e:
                print(f"[Crossmatch] TAP bibcode fetch warning for {main_id}: {e}")

        rows.append({
            "rank": c["rank"],
            "specObjID": c["specObjID"],
            "simbad_main_id": main_id,
            "simbad_otype": str(r["otype"]).strip() if "otype" in res.colnames else None,
            "simbad_z": rvz,
            "simbad_morphtype": morph,
            "simbad_ids": ids,
            "simbad_n_bibcodes": len(bib_list) if bib_list else nbref,
            "simbad_bibcodes": ";".join(bib_list),
            "simbad_url": "https://simbad.cds.unistra.fr/simbad/sim-id?Ident=" + main_id.replace(" ", "+"),
        })

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    out.to_csv(out_csv, index=False)
    print(f"[Crossmatch] SIMBAD results for {len(out)} candidates -> {out_csv}")
    return out


if __name__ == "__main__":
    candidates_csv = "data/sdss_anomaly_candidates.csv"
    if os.path.exists(candidates_csv):
        df_c = pd.read_csv(candidates_csv)
        crossmatch_candidates(df_c)
    else:
        print(f"Candidates file {candidates_csv} not found.")
