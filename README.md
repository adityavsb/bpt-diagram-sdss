# Mapping the Cosmic Boundary: Starburst Galaxies vs. Accreting Black Holes (SDSS DR17)

[![Python](https://img.shields.io/badge/Python-3.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://www.python.org/)
[![Astropy](https://img.shields.io/badge/Astropy-6.x%2B-orange.svg)](https://www.astropy.org/)
[![Astroquery](https://img.shields.io/badge/Astroquery-SDSS%20DR17-purple.svg)](https://astroquery.readthedocs.io/)
[![Survey](https://img.shields.io/badge/Survey-SDSS%20DR17%20%2F%20MPA--JHU-teal.svg)](https://www.sdss.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An observational astrophysics project implementing the **Baldwin, Phillips & Terlevich (BPT, 1981)** emission-line diagnostic on ~10,000 galaxies from the **Sloan Digital Sky Survey**, using the MPA-JHU emission-line measurements (Brinchmann et al. 2004) of SDSS Legacy spectra as distributed in the `galSpecLine` table of **SDSS DR17**.

The project investigates a fundamental question in extragalactic astronomy:  
*How do we know whether a galaxy’s glow is driven by newborn stars or by a feeding supermassive black hole?*

---

## The Astrophysical Principle: Reading the Gas

When high-energy radiation strikes interstellar gas, atoms ionize and recombine, emitting photons at discrete spectral wavelengths. By taking the logarithmic ratio of specific forbidden collision lines to hydrogen recombination lines, we construct an excitation diagnostic that is independent of distance or dust reddening:

* **y-axis:** **log₁₀([O III] λ5007 / Hβ)** — tracks the ionization parameter and radiation hardness of the ionizing field.
* **x-axis:** **log₁₀([N II] λ6584 / Hα)** — tracks gas-phase metallicity and nitrogen abundance.

```text
               High [O III] / H-beta  (Hard Ionizing Photons)
                         ^
                         |       ^
                         |      /   AGN Wing
                         |     /    [Hard UV/X-ray accretion]
                         |    /
                         |   /
                         |  /
       Star-Forming      | /       Transition (Composite)
       (OB Star Nurseries)
                         |
                         v
        Low [N II] / H-alpha <-------------------> High [N II] / H-alpha
                                           (High Metallicity)
```

---

## Diagnostic Boundaries

We implement the standard theoretical and empirical classification curves from the literature on the optical $[\text{N II}]\lambda 6584/\text{H}\alpha$ vs. $[\text{O III}]\lambda 5007/\text{H}\beta$ diagnostic plane:

1. **Kauffmann et al. (2003) Pure Star-Formation Limit:**

$$
\log_{10}\left(\frac{[\text{O III}]}{\text{H}\beta}\right) = \frac{0.61}{\log_{10}([\text{N II}]/\text{H}\alpha) - 0.05} + 1.30 \quad (\text{for } x < 0.05)
$$

   Galaxies below this dashed curve are powered solely by young massive OB stars (< 10 Myr).

2. **Kewley et al. (2001) Extreme Starburst Theoretical Limit:**

$$
\log_{10}\left(\frac{[\text{O III}]}{\text{H}\beta}\right) = \frac{0.61}{\log_{10}([\text{N II}]/\text{H}\alpha) - 0.47} + 1.19 \quad (\text{for } x < 0.47)
$$

   Galaxies above this solid curve cannot be ionized by starbursts alone (it is the upper envelope of the Kewley et al. 2001 grid of starburst photoionization models spanning sub-solar to super-solar metallicity and a wide range of ionization parameters, their "maximum starburst line"), requiring hard ionizing radiation from an Active Galactic Nucleus (AGN) or radiative shocks.

3. **Composite Transition Sequence:**
   Galaxies situated between the Kauffmann (2003) and Kewley (2001) curves exhibit intermediate line ratios resulting from mixed excitation—typically stellar nurseries co-existing with weak nuclear activity or low-velocity shocks.

> [!NOTE]
> **On Seyfert vs. LINER separation.** Kewley et al. (2006) draw the Seyfert/LINER line on the [S II]/Hα and [O I]/Hα diagrams (their Eqs. 1–2), where the two populations separate most cleanly; on the [N II]/Hα diagram they overlap more. [N II]-only separations do exist — Schawinski et al. (2007): $\log([\text{O III}]/\text{H}\beta) = 1.05\,\log([\text{N II}]/\text{H}\alpha) + 0.45$; Cid Fernandes et al. (2010): $1.01\,x + 0.48$ — and the WHAN diagram (Cid Fernandes et al. 2010) separates true LINERs from "retired" galaxies ionised by old stars. This project follows Kewley et al. (2006) and reports three classes; adding [S II] and [O I] is the planned Iteration 3.

---

## Observational Results (10,000 galaxies, deterministic hash-ordered sample)

Selection: `class = GALAXY`, `zWarning = 0`, `sciencePrimary = 1`, MPA-JHU `reliable = 1`, $0.02 < z < 0.25$, and $S/N \ge 3$ on all four lines after applying the MPA-JHU recommended error rescaling (Hα ×2.473, Hβ ×1.882, [O III] ×1.566, [N II] ×2.039). The sample is the first 10,000 rows of the selection ordered by a hash of `specObjID`, i.e. a reproducible pseudo-random draw across the full Legacy footprint, not a complete census. Percentages below describe this sample.

| Classification | Galaxy Count | Population % | Dominant Physical Mechanism |
| :--- | :---: | :---: | :--- |
| **Star-Forming** | `7,344` | **`73.4%`** | Photoionization by short-lived OB stars in stellar nurseries |
| **Composite** | `1,783` | **`17.8%`** | Transition systems undergoing simultaneous starbursts and AGN accretion |
| **AGN** | `873` | **`8.7%`** | Hard power-law ionizing radiation from supermassive black hole accretion / shocks |

### Validation against MPA-JHU

Cross-tabulating our classification against the official MPA-JHU `bptclass` column (`galSpecExtra`) across the 10,000 galaxies yields 100.0% agreement (7,344 Star-Forming, 1,783 Composite, 873 AGN; see `data/classifier_validation.csv`).

To verify that my implementation of the classical Kewley et al. (2001) and Kauffmann et al. (2003) demarcation equations was numerically sound, I cross-validated the entire 10,000-galaxy sample against the official MPA-JHU `bptclass` catalog flags. Achieving 100.0% classification agreement across all 7,344 Star-Forming, 1,783 Composite, and 873 AGN galaxies verified that applying standard Balmer/forbidden width groupings and recommended S/N thresholds exactly reproduces the standard astronomical catalog benchmark.

![BPT Diagnostic Diagram](assets/bpt_diagnostic_sdss.png)

---

## Two Anchor Galaxies: M51 vs. NGC 1068

Nuclear line ratios from the Palomar survey of nearby galaxy nuclei (Ho, Filippenko & Sargent 1997, ApJS, 112, 315), which observed the central few arcseconds of each galaxy, the same kind of aperture an SDSS fiber provides at the sample's redshifts.

| Anchor Galaxy | HFS97 nuclear class | log([N II]/Hα) | log([O III]/Hβ) | Where it falls |
| :--- | :--- | :---: | :---: | :--- |
| **M51 / NGC 5194** | Seyfert 2 (earlier works: LINER) | `+0.07` | `+0.22` | AGN. The famous star-forming spiral arms are outside a nuclear aperture; the *nucleus* hosts a low-luminosity AGN. |
| **NGC 1068 / M77** | Seyfert 2 | `+0.17` | `+1.10` | AGN. The archetypal obscured AGN. |

---

## Project Structure

```text
bpt-diagram-sdss/
+-- README.md               # Scientific documentation & observational results
+-- PLAN.md                 # Technical architecture & astrophysics plan
+-- CITATION.cff            # Citation metadata for academic and open-source citing
+-- REFERENCES.bib          # NASA/ADS BibTeX database for LaTeX / paper references
+-- NOTES.md                # Working notes, decisions, and corrections
+-- LICENSE                 # MIT License
+-- .gitignore              # Defense-in-depth rules (excludes FITS caches & secrets)
+-- requirements.txt        # astropy, astroquery, pandas, matplotlib, numpy, scikit-learn
+-- src/
|   +-- bpt_math.py         # Boundary equations and MPA-JHU validation
|   +-- query_data.py       # SDSS DR17 query & local caching with quality cuts
|   +-- plot_bpt.py         # 2D hexbin density visualization & demographics generator
|   +-- glitch_filter.py    # Analytical Gaussian error propagation & 95% CI modeling
|   +-- anomaly_detector.py # Ensemble outlier detector (iForest + LOF + PCA) & stability
|   +-- crossmatch.py       # Automated SIMBAD cross-match resolver via astroquery
|   \-- plot_anomalies.py   # Anomaly diagnostic map with 2-D error ellipses
+-- assets/
|   +-- bpt_diagnostic_sdss.png # Classical BPT diagram
|   \-- bpt_anomalies_sdss.png  # Outlier detection diagnostic map
+-- data/
|   +-- sdss_bpt_sample.csv                  # Tracked sample of 10,000 galaxies (reproducibility)
|   +-- sdss_anomaly_candidates.csv          # Outlier candidates for follow-up (high ML consensus, low fragility)
|   +-- sdss_anomaly_candidates_simbad.csv   # SIMBAD cross-match results and reference counts
|   +-- classifier_validation.csv            # Cross-tabulation vs. official MPA-JHU bptclass (100% agreement)
|   +-- anomaly_stability.csv                # Stability across 15 (seed, k) runs
|   \-- anchor_galaxies.json                 # Palomar HFS97 nuclear emission line ratios
# Derived from public SDSS DR17 / MPA-JHU catalog data; redistributed under the SDSS data policy with the acknowledgements below.
\-- run_bpt.py              # Master pipeline runner (Classical BPT + Outlier Detection)
```

---

## Iteration 2: Ensemble line-ratio outlier detection

Classical demarcation curves partition the diagram; objects far from the bulk population in line-ratio space are worth a second look, whether because of unusual physics or because of measurement problems. Iteration 2 ranks objects by how unusual their *line ratios* are, using three unsupervised detectors on four distance-independent features: $\log([\text{N II}]/\text{H}\alpha)$, $\log([\text{O III}]/\text{H}\beta)$, $\log(\text{H}\alpha/\text{H}\beta)$ and $\log \text{EW}(\text{H}\alpha)$. Absolute fluxes are excluded on purpose. The output is a **candidate list for follow-up**, not a list of new astrophysical claims.

![SDSS DR17 BPT Anomaly Detection](assets/bpt_anomalies_sdss.png)

### 1. The Three-Algorithm Voting Ensemble

To eliminate single-algorithm bias, we employ an equal-weighted consensus ensemble across three complementary unsupervised paradigms:
1. **Isolation Forest (`iForest`)**: Measures path-partitioning depth in feature space to identify global spectral outliers.
2. **Local Outlier Factor (`LOF`, k = 35; k = 20 and 50 used for the stability check)**: Density ratio relative to the k nearest neighbours in feature space.
3. **PCA reconstruction residual** (components chosen to explain ≥ 90% of variance; see console output): Distance from the low-dimensional manifold traced by the bulk population.

Each raw score is calibrated to a percentile scale ($0\text{--}100$). The **Consensus Anomaly Score** is their equal-weighted average:

$$
S_{\text{consensus}} = \frac{1}{3} S_{\text{iForest}} + \frac{1}{3} S_{\text{LOF}} + \frac{1}{3} S_{\text{PCA}}
$$

### 2. Measurement fragility and uncertainty

In survey astronomy, the raw top outliers of any unsupervised algorithm are often **data artifacts** (e.g., cosmic ray hits, sky subtraction residuals, or faint-line flux denominators approaching zero).

To evaluate how sensitive each outlier is to measurement uncertainty, we implement:
* **Analytical Gaussian Error Propagation:**

$$
\sigma_x = \frac{1}{\ln 10} \sqrt{\left(\frac{\sigma_{[\text{N II}]}}{[\text{N II}]}\right)^2 + \left(\frac{\sigma_{\text{H}\alpha}}{\text{H}\alpha}\right)^2}, \quad \sigma_y = \frac{1}{\ln 10} \sqrt{\left(\frac{\sigma_{[\text{O III}]}}{[\text{O III}]}\right)^2 + \left(\frac{\sigma_{\text{H}\beta}}{\text{H}\beta}\right)^2}
$$

* **Uncertainty geometry:** 1-D marginal 95% intervals are $x \pm 1.96\,\sigma_x$, $y \pm 1.96\,\sigma_y$ (the `[lo, hi]` strings in the table). The plotted 2-D 95% error ellipses use semi-axes $2.448\,\sigma_x$, $2.448\,\sigma_y$ ($\sqrt{\chi^2_{2}(0.95)}$); a 1.96σ box would enclose only ~85% of a 2-D Gaussian. Errors in the two ratios are treated as uncorrelated (the MPA-JHU fits share a Balmer and a forbidden-line width, so this is an approximation).
* **`measurement_fragility` (0–1, heuristic):** rises toward 1 as the minimum rescaled line S/N approaches the selection floor, as the 2-D error ellipse grows relative to the sample's 95th-percentile ellipse, and as Hβ/[O III] S/N falls. It is not a probability.
* **`quality_weighted_outlier_score` (0–100):** consensus percentile × (1 − fragility).
* **`n_sigma_nearest_boundary`:** distance, in units of the object's own σ, to the nearest demarcation curve. Objects with ≥ 3 keep their class under a 3σ fluctuation.
* **Stability:** the ensemble is re-run for 5 Isolation Forest seeds × 3 LOF neighbourhood sizes; a candidate must be in the top 50 in every run.

#### Stability across seeds and neighbourhoods

Evaluating the ensemble across 15 runs (5 Isolation Forest seeds × 3 LOF neighbourhood sizes $k \in \{20, 35, 50\}$) yields a mean pairwise Spearman correlation of 0.990 and mean top-50 Jaccard overlap of 0.893 (see `data/anomaly_stability.csv`). 42 galaxies appear in the top 50 in all 15 runs.

Unsupervised outlier detectors can often be brittle, easily perturbed by random tree partitions or neighborhood parameter choices. By testing the ensemble across 15 independent runs (varying Isolation Forest random seeds and LOF $k$-neighborhoods), I verified that the consensus ranking is highly consistent (mean Spearman rank correlation of 0.990, top-50 Jaccard overlap of 0.893). Exactly 42 galaxies remained in the top 50 across every single run, allowing me to filter out fragile borderline cases and retain only the 8 candidates that meet strict stability and physical signal-to-noise criteria.

### 3. Top outlier candidates for follow-up

Selection: `measurement_fragility <= 0.25`, minimum rescaled S/N >= 10, >= 2 of 3 algorithm votes, stable in 15/15 runs; ranked by `quality_weighted_outlier_score`. These are statistically unusual line-ratio combinations; none has been verified as physically unusual by this project.

| Rank | SpecObjID | z | log([N II]/Hα) | log([O III]/Hβ) | 1-D 95% CI | Consensus | QW score | Fragility | n_σ boundary | BPT region | Driving feature | Links |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| **#1** | `1522331895189235712` | `0.235` | `-0.17` | `+0.93` | `[-0.20, -0.14]`, `[0.89, 0.96]` | `99.6` | **`99.0`** | `0.007` | `30.8σ` | AGN region, very high [O III]/Hbeta (y >= 0.9) | `log_oiii_hb` (z=+2.35) | [SDSS](https://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=1522331895189235712) • [Legacy](https://www.legacysurvey.org/viewer?ra=225.57194&dec=38.05195&layer=ls-dr10&zoom=15) |
| **#2** | `2237320683103217664` | `0.033` | `-1.58` | `+0.61` | `[-1.64, -1.52]`, `[0.59, 0.63]` | `98.5` | **`97.6`** | `0.009` | `22.0σ` | Very low [N II]/Halpha (low [N II]/Hα; unverified) | `log_nii_ha` (z=+4.75) | [SDSS](https://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=2237320683103217664) • [Legacy](https://www.legacysurvey.org/viewer?ra=192.81909&dec=35.11345&layer=ls-dr10&zoom=15) |
| **#3** | `2293593682860861440` | `0.230` | `-0.97` | `+0.52` | `[-1.00, -0.93]`, `[0.51, 0.54]` | `97.5` | **`97.1`** | `0.004` | `12.2σ` | Star-forming region | `log_ew_ha` (z=+3.12) | [SDSS](https://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=2293593682860861440) • [Legacy](https://www.legacysurvey.org/viewer?ra=170.60098&dec=37.76703&layer=ls-dr10&zoom=15) |
| **#4** | `1432306634738132992` | `0.176` | `-0.35` | `+0.17` | `[-0.37, -0.33]`, `[0.13, 0.20]` | `97.5` | **`97.0`** | `0.005` | `11.1σ` | Composite region | `log_ew_ha` (z=+2.31) | [SDSS](https://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=1432306634738132992) • [Legacy](https://www.legacysurvey.org/viewer?ra=137.25270&dec=33.43255&layer=ls-dr10&zoom=15) |
| **#5** | `1679928197229078528` | `0.039` | `-0.65` | `+0.16` | `[-0.72, -0.58]`, `[0.09, 0.23]` | `97.6` | **`94.1`** | `0.036` | `4.1σ` | Star-forming region | `log_ew_ha` (z=+1.99) | [SDSS](https://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=1679928197229078528) • [Legacy](https://www.legacysurvey.org/viewer?ra=4.00171&dec=-0.30269&layer=ls-dr10&zoom=15) |
| **#6** | `2247426838710216704` | `0.035` | `-1.71` | `+0.66` | `[-1.78, -1.63]`, `[0.64, 0.68]` | `99.5` | **`91.1`** | `0.085` | `18.6σ` | Very low [N II]/Halpha (low [N II]/Hα; unverified) | `log_nii_ha` (z=+5.29) | [SDSS](https://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=2247426838710216704) • [Legacy](https://www.legacysurvey.org/viewer?ra=171.18009&dec=39.76517&layer=ls-dr10&zoom=15) |
| **#7** | `512391142269020160` | `0.031` | `-0.62` | `+0.86` | `[-0.70, -0.53]`, `[0.78, 0.94]` | `97.7` | **`85.0`** | `0.130` | `5.2σ` | Outlier at ordinary BPT position | `log_oiii_hb` (z=+2.21) | [SDSS](https://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=512391142269020160) • [Legacy](https://www.legacysurvey.org/viewer?ra=38.29496&dec=-7.80372&layer=ls-dr10&zoom=15) |
| **#8** | `2518725292126660608` | `0.051` | `-0.35` | `-0.82` | `[-0.38, -0.33]`, `[-0.91, -0.74]` | `96.9` | **`80.8`** | `0.166` | `8.2σ` | Star-forming region | `log_oiii_hb` (z=+1.29) | [SDSS](https://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=2518725292126660608) • [Legacy](https://www.legacysurvey.org/viewer?ra=188.59066&dec=29.65590&layer=ls-dr10&zoom=15) |

*Full candidate catalog exported with complete coordinates and multi-wavelength links in [`data/sdss_anomaly_candidates.csv`](data/sdss_anomaly_candidates.csv).*

---

## What SIMBAD says about the top candidates

The table below is generated from `data/sdss_anomaly_candidates_simbad.csv` (5″ cone match, queried 2026-10-28). "Driving feature" is the scaled feature that deviates most from the sample median for that object (see §Iteration 2). Where a published paper discusses the object individually, it is cited; catalog-level inclusions are not described as studies of the object.

| Rank | SIMBAD main ID | SIMBAD type | z | BPT region | Driving feature (standardized z) | SIMBAD refs | Notes |
| :---: | :--- | :--- | :---: | :--- | :--- | :---: | :--- |
| **#1** | `2MASS J15021725+3803069` | Sy2 | `0.235` | AGN region, very high [O III]/Hbeta (y >= 0.9) | `log_oiii_hb` (+2.35) | 7 | Seyfert 2 (SIMBAD). Optically selected Type 2 quasar candidate studied in Reyes et al. (2008) and Kong & Ho (2018). Radio-detected in FIRST. |
| **#2** | `LEDA 4572013` | EmG | `0.033` | Very low [N II]/Halpha (low [N II]/Hα; unverified) | `log_nii_ha` (+4.75) | 6 | Emission-line galaxy (SIMBAD; HS 1248+3523 in Hopp et al. 2000). Cataloged in the ALFALFA H I survey (Haynes et al. 2018). Extremely low [N II]/Hα. |
| **#3** | *No cross-match* | -- | `0.230` | Star-forming region | `log_ew_ha` (+3.12) | 0 | No cross-match within 5″ in SIMBAD. High equivalent width Hα emission in an otherwise standard star-forming locus. |
| **#4** | *No cross-match* | -- | `0.176` | Composite region | `log_ew_ha` (+2.31) | 0 | No cross-match within 5″ in SIMBAD. Sits in composite transition track; flagged for elevated equivalent width Hα. |
| **#5** | *No cross-match* | -- | `0.039` | Star-forming region | `log_ew_ha` (+1.99) | 0 | No cross-match within 5″ in SIMBAD. Low-redshift star-forming galaxy flagged for elevated equivalent width Hα. |
| **#6** | *No cross-match* | -- | `0.035` | Very low [N II]/Halpha (low [N II]/Hα; unverified) | `log_nii_ha` (+5.29) | 0 | No cross-match within 5″ in SIMBAD. Extreme low-metallicity line ratio signature (log([N II]/Hα) = -1.71). |
| **#7** | `2MASX J02331076-0748131` | Sy2 | `0.031` | Outlier at ordinary BPT position | `log_oiii_hb` (+2.21) | 16 | Seyfert 2 (SIMBAD). Identified as harboring a low-mass black hole in Greene & Ho (2007) and Jiang et al. (2011). X-ray source in Chandra survey (Bi et al. 2020). |
| **#8** | `2MASX J12342176+2939215` | SBG | `0.051` | Star-forming region | `log_oiii_hb` (+1.29) | 11 | Starburst galaxy with morphological type Sc (SIMBAD). Cataloged in KPNO International Spectroscopic Survey (KISSR 78; Salzer et al. 2001, Wegner et al. 2003). |

**Rules used to write the Notes column** (kept here so readers can judge them): (1) a paper is cited only if its bibcode is in SIMBAD's reference list for the object *and* the paper's text discusses the object individually; (2) survey membership is stated with the paper that defines the sample; (3) no physical mechanism (outflow, torus, IMBH, rejuvenation) is asserted unless a cited paper asserts it for this object; (4) morphology comes from SIMBAD's `morphtype` field with its reference.

## Scope and limitations

* **What this is:** a reproducible implementation of the [N II]/Hα BPT classification on a hash-ordered sample of SDSS Legacy galaxies with MPA-JHU line measurements, plus an unsupervised ranking of unusual line-ratio combinations with propagated uncertainties.
* **What this is not:** a paper claiming new astrophysical sources. No candidate has been followed up; "outlier" means statistically unusual relative to this sample in four line-derived features, not verified through dedicated observations.
* **Known approximations:** errors in the two ratios are treated as uncorrelated; the MPA-JHU error rescaling factors are sample-wide averages; the fragility score is a hand-tuned heuristic; the sample is 10,000 objects, not the full survey.
* **Validation performed:** agreement with the MPA-JHU `bptclass` column (`data/classifier_validation.csv`); seed/neighbourhood stability of the outlier ranking (`data/anomaly_stability.csv`); SIMBAD cross-match of candidates (`data/sdss_anomaly_candidates_simbad.csv`); unit tests (`pytest`).
* **Literature statements** about individual candidates are limited to papers that SIMBAD associates with the object and that discuss it individually; see §What SIMBAD says about the top candidates.
* **Working notes:** See [`NOTES.md`](NOTES.md) for dated notes on design decisions, corrections, and dropped candidates.

---

## What I learned

Developing this pipeline and conducting a thorough code and methodology audit challenged several of my initial assumptions about machine learning in observational astronomy:

1. **Unsupervised ML without domain constraints optimizes for catalog artifacts:** In my first iteration, I tried to leverage machine learning by feeding the models as much data as possible—including raw emission-line fluxes alongside line ratios. Even though "garbage in, garbage out" is a standard rule, I assumed tree partitions and distance metrics would naturally handle the multiple dimensions. Instead, the models took the path of least resistance: they voted for the brightest objects in the catalog with extreme, uncalibrated raw fluxes ($> 10^9$) and local coordinates, flagging bright lightbulbs rather than exotic astrophysics. Permanently adopting a 4D scale-invariant feature space (ratios, Balmer decrement, and equivalent width) was essential to ground the detector in intrinsic gas thermodynamics.
2. **The danger of confirmation bias driven by Malmquist bias:** When my initial candidates were cross-referenced against NED and SIMBAD, finding dozens of papers discussing radio cores and starbursts gave me false confidence that the algorithm was pursuing the right direction. In reality, those objects had extensive literature simply because they were historically prominent, bright targets (astronomical Malmquist bias); on the BPT diagram, their line ratios sat squarely inside the ordinary composite cloud. True validation requires inspecting feature attributions and stability metrics, not celebrating search hits on famous galaxies.
3. **Statistical geometry in multi-dimensional astronomy:** When plotting 2-D error ellipses, I initially used a $1.96\sigma$ radius—a common mistake carried over from 1-D normal distributions. In a 2-D joint Gaussian distribution, the squared distance follows a $\chi^2_2$ distribution, where $1.96\sigma$ encloses only ~85% of the probability. Adjusting the ellipse multiplier to $\sqrt{\chi^2_2(0.95)} \approx 2.448\sigma$ is critical to avoid classifying 15% of ordinary galaxies as false outliers due to measurement noise.

---

## Quickstart & Reproduction

### 1. Set Up Environment

```bash
git clone https://github.com/adityavsb/bpt-diagram-sdss.git
cd bpt-diagram-sdss

# Create and activate virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

### 2. Run the Diagnostic & Outlier Detection Pipeline

```bash
python run_bpt.py
```
*Executes the end-to-end pipeline in ~4 seconds: loads the tracked sample from `data/` (pass `--force-refresh` to re-query SkyServer), performs classical BPT classification, runs the ensemble outlier detector, computes 95% confidence intervals, and outputs figures to `assets/`.*

### 3. Run Automated Tests

```bash
pip install -r requirements-dev.txt
pytest -q
```
*Executes regression tests covering theoretical boundaries, BPT region classification, error propagation, confidence ellipses, and Mahalanobis boundary detachment.*

---

## Scientific References & Literature Attribution

This project implements analytical boundaries, survey datasets, and observational findings from the following literature:

1. **The BPT Diagnostic Foundation:**  
   Baldwin, J. A., Phillips, M. M., & Terlevich, R. (1981). *Classification parameters for the emission-line spectra of extragalactic objects.* **PASP**, 93, 5–19. [DOI: 10.1086/130766](https://doi.org/10.1086/130766)
2. **Empirical Star-Formation Demarcation:**  
   Kauffmann, G., Heckman, T. M., Tremonti, C., et al. (2003). *The host galaxies of active galactic nuclei.* **MNRAS**, 346(4), 1055–1077. [DOI: 10.1111/j.1365-2966.2003.07154.x](https://doi.org/10.1111/j.1365-2966.2003.07154.x)
3. **Extreme Starburst Theoretical Limit:**  
   Kewley, L. J., Dopita, M. A., Sutherland, R. S., et al. (2001). *Theoretical Modeling of Starburst Galaxies.* **ApJ**, 556(1), 121–140. [DOI: 10.1086/321545](https://doi.org/10.1086/321545)
4. **AGN Host Galaxies & Multi-Line Diagnostics:**  
   Kewley, L. J., Groves, B., Kauffmann, G., & Heckman, T. (2006). *The host galaxies and classification of active galactic nuclei.* **MNRAS**, 372(3), 961–976. [DOI: 10.1111/j.1365-2966.2006.10859.x](https://doi.org/10.1111/j.1365-2966.2006.10859.x)
5. **MPA-JHU Emission-Line Measurements:**  
   Brinchmann, J., Charlot, S., White, S. D. M., et al. (2004). *The physical properties of star-forming galaxies in the low-redshift Universe.* **MNRAS**, 351(4), 1151–1179. [DOI: 10.1111/j.1365-2966.2004.07881.x](https://doi.org/10.1111/j.1365-2966.2004.07881.x)
6. **SDSS-IV Data Release 17:**  
   Abdurro'uf, Accetta, K., Aerts, C., et al. (2022). *The Seventeenth Data Release of the Sloan Digital Sky Surveys: Complete Release of MaNGA, MaStar, and APOGEE-2 Data.* **ApJS**, 259(2), 35. [DOI: 10.3847/1538-4365/ac4414](https://doi.org/10.3847/1538-4365/ac4414)
7. **Palomar Survey Nuclear Emission Line Ratios (Anchor Galaxies):**  
   Ho, L. C., Filippenko, A. V., & Sargent, W. L. W. (1997). *A Search for "Dwarf" Seyfert Nuclei. III. Spectroscopic Parameters and Properties of the Host Galaxies.* **ApJS**, 112(2), 315–390. [DOI: 10.1086/313041](https://doi.org/10.1086/313041)
8. **AGN Feedback & [N II] Seyfert/LINER Demarcation:**  
   Schawinski, K., Thomas, D., Sarzi, M., et al. (2007). *Observational evidence for AGN feedback in early-type galaxies.* **MNRAS**, 382(4), 1415–1431. [DOI: 10.1111/j.1365-2966.2007.12487.x](https://doi.org/10.1111/j.1365-2966.2007.12487.x)
9. **Alternative Diagnostic Diagrams & WHAN Classification:**  
   Cid Fernandes, R., Stasińska, G., Schlickmann, M. S., et al. (2010). *Alternative diagnostic diagrams and the 'forgotten' population of weak line galaxies in the SDSS.* **MNRAS**, 403(2), 1036–1053. [DOI: 10.1111/j.1365-2966.2009.16185.x](https://doi.org/10.1111/j.1365-2966.2009.16185.x)
10. **Type 2 Quasar Selection (Candidate #1):**  
    Reyes, R., Zakamska, N. L., Strauss, M. A., et al. (2008). *Space Density of Optically Selected Type 2 Quasars.* **AJ**, 136(6), 2373–2390. [DOI: 10.1088/0004-6256/136/6/2373](https://doi.org/10.1088/0004-6256/136/6/2373)
11. **Low-Mass Black Holes in Active Galaxies (Candidate #7):**  
    Greene, J. E., & Ho, L. C. (2007). *A New Sample of Low-Mass Black Holes in Active Galaxies.* **ApJ**, 670(1), 92–104. [DOI: 10.1086/522082](https://doi.org/10.1086/522082)
12. **KPNO International Spectroscopic Survey (Candidate #8):**  
    Salzer, J. J., Gronwall, C., Lipovetsky, V. A., et al. (2001). *The KPNO International Spectroscopic Survey. II. Hα-selected Survey List 1.* **AJ**, 121(1), 66–85. [DOI: 10.1086/318040](https://doi.org/10.1086/318040)

Full BibTeX records are available in [`REFERENCES.bib`](REFERENCES.bib).

---

## Data Rights & Database Acknowledgements

1. **Sloan Digital Sky Survey IV (SDSS-IV):**  
   Funding for the Sloan Digital Sky Survey IV has been provided by the Alfred P. Sloan Foundation, the U.S. Department of Energy Office of Science, and the Participating Institutions. SDSS-IV is managed by the Astrophysical Research Consortium for the Participating Institutions of the SDSS Collaboration ([www.sdss.org](https://www.sdss.org)).

2. **NASA/IPAC Extragalactic Database (NED):**  
   This research has made use of the NASA/IPAC Extragalactic Database (NED), funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology.

3. **CDS SIMBAD Astronomical Database:**  
   This research has made use of the SIMBAD database, operated at CDS, Strasbourg, France.

---

## License & Academic Attribution Policy

* **Code License:** [MIT License](LICENSE) © 2026 Aditya Bhadoria.
* **Primary Scientific Attribution (Requested):** This repository is an educational reproduction and open-source demonstration. Any scientific publication, paper, or research work utilizing the BPT diagnostic, classification boundaries, or SDSS catalog **should cite the original literature** listed in [Scientific References](#scientific-references--literature-attribution) and [`REFERENCES.bib`](REFERENCES.bib) (Baldwin et al. 1981; Kauffmann et al. 2003; Kewley et al. 2001, 2006; Brinchmann et al. 2004; Abdurro'uf et al. 2022).
* **Software Implementation Citation (Optional):** If you specifically adapt or reuse this Python pipeline, code structure, or plotting workflow in an educational tutorial, university lab, or coding demo, you may additionally cite the software repository via [`CITATION.cff`](CITATION.cff) or GitHub's **"Cite this repository"** button.
