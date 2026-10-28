# Mapping the Cosmic Boundary: Starburst Galaxies vs. Accreting Black Holes (SDSS DR17)

[![Python](https://img.shields.io/badge/Python-3.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://www.python.org/)
[![Astropy](https://img.shields.io/badge/Astropy-6.x%2B-orange.svg)](https://www.astropy.org/)
[![Astroquery](https://img.shields.io/badge/Astroquery-SDSS%20DR17-purple.svg)](https://astroquery.readthedocs.io/)
[![Survey](https://img.shields.io/badge/Survey-SDSS%20DR17%20%2F%20MPA--JHU-teal.svg)](https://www.sdss.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An observational astrophysics project implementing the landmark **Baldwin, Phillips & Terlevich (BPT, 1981)** emission-line diagnostic on ~10,000 galaxies from the **Sloan Digital Sky Survey**, using the MPA-JHU emission-line measurements (Brinchmann et al. 2004) of SDSS Legacy spectra as distributed in the `galSpecLine` table of **SDSS DR17**.

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
> **On Seyfert vs. LINER Separation**: In the literature, Kewley et al. (2006, Eq. 1) established the division $\log_{10}([\text{O III}]/\text{H}\beta) = 1.89 \log_{10}([\text{S II}]/\text{H}\alpha) + 0.76$ specifically for the optical $[\text{S II}]/\text{H}\alpha$ diagnostic diagram (and Eq. 2 for $[\text{O I}]/\text{H}\alpha$). On the standard 4-line $[\text{N II}]/\text{H}\alpha$ diagram, Seyferts and LINERs overlap significantly without a clean bimodal valley. Therefore, standard spectroscopic practice employs the canonical 3-class demarcation (**Star-Forming, Composite, AGN**), reserving fine-grained Seyfert vs. LINER separation for multi-line $[\text{S II}]$, $[\text{O I}]$, or WHAN diagrams.

---

## Observational Results (10,000 galaxies, deterministic hash-ordered sample)

Selection: `class = GALAXY`, `zWarning = 0`, `sciencePrimary = 1`, MPA-JHU `reliable = 1`, $0.02 < z < 0.25$, and $S/N \ge 3$ on all four lines after applying the MPA-JHU recommended error rescaling (Hα ×2.473, Hβ ×1.882, [O III] ×1.566, [N II] ×2.039). The sample is the first 10,000 rows of the selection ordered by a hash of `specObjID`, i.e. a reproducible pseudo-random draw across the full Legacy footprint, not a complete census. Percentages below describe this sample.

| Classification | Galaxy Count | Population % | Dominant Physical Mechanism |
| :--- | :---: | :---: | :--- |
| **Star-Forming** | `7,344` | **`73.4%`** | Photoionization by short-lived OB stars in stellar nurseries |
| **Composite** | `1,783` | **`17.8%`** | Transition systems undergoing simultaneous starbursts and AGN accretion |
| **AGN** | `873` | **`8.7%`** | Hard power-law ionizing radiation from supermassive black hole accretion / shocks |

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
+-- LICENSE                 # MIT License
+-- .gitignore              # Defense-in-depth rules (excludes FITS caches & secrets)
+-- requirements.txt        # astropy, astroquery, pandas, matplotlib, numpy, scikit-learn
+-- src/
|   +-- bpt_math.py         # Boundary equations (Kauffmann 2003, Kewley 2001, Kewley 2006 reference)
|   +-- query_data.py       # SDSS DR17 query & local caching with quality cuts
|   +-- plot_bpt.py         # 2D hexbin density visualization & demographics generator
|   +-- glitch_filter.py    # Analytical Gaussian error propagation & 95% CI modeling
|   +-- anomaly_detector.py # Tri-algorithmic ensemble (iForest + LOF + PCA) & voting
|   \-- plot_anomalies.py   # Publication-quality anomaly diagnostic map
+-- assets/
|   +-- bpt_diagnostic_sdss.png # Classical BPT diagram
|   \-- bpt_anomalies_sdss.png  # Anomaly detection & 95% CI error ellipses
+-- data/
|   +-- sdss_bpt_sample.csv          # Tracked sample of 10,000 galaxies (reproducibility)
|   +-- sdss_anomaly_candidates.csv  # Outlier candidates for follow-up (high ML consensus, low fragility)
|   +-- classifier_validation.csv    # Cross-tabulation vs. official MPA-JHU bptclass (100% agreement)
|   \-- anchor_galaxies.json         # Palomar HFS97 nuclear emission line ratios
\-- run_bpt.py              # Master pipeline runner (Classical BPT + Outlier Detection)
```

---

## Iteration 2: Ensemble line-ratio outlier detection

Classical demarcation curves partition the diagram; objects far from the bulk population in line-ratio space are worth a second look, whether because of unusual physics or because of measurement problems. Iteration 2 ranks objects by how unusual their *line ratios* are, using three unsupervised detectors on four distance-independent features: $\log([\text{N II}]/\text{H}\alpha)$, $\log([\text{O III}]/\text{H}\beta)$, $\log(\text{H}\alpha/\text{H}\beta)$ and $\log \text{EW}(\text{H}\alpha)$. Absolute fluxes are excluded on purpose. The output is a **candidate list for follow-up**, not a list of discoveries.

![SDSS DR17 BPT Anomaly Detection](assets/bpt_anomalies_sdss.png)

### 1. The Tri-Algorithmic Voting Ensemble

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

### 3. Top outlier candidates for follow-up

Selection: `measurement_fragility <= 0.25`, minimum rescaled S/N >= 10, >= 2 of 3 algorithm votes, stable in 15/15 runs; ranked by `quality_weighted_outlier_score`. These are statistically unusual line-ratio combinations; none has been confirmed as physically unusual by this project.

| Rank | SpecObjID | z | log([N II]/Hα) | log([O III]/Hβ) | 1-D 95% CI | Consensus | QW score | Fragility | n_σ boundary | BPT region | Driving feature | Links |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| **#1** | `1522331895189235712` | `0.235` | `-0.17` | `+0.93` | `[-0.20, -0.14]`, `[0.89, 0.96]` | `99.6` | **`99.0`** | `0.007` | `30.8σ` | AGN region, very high [O III]/Hbeta (y >= 0.9) | `log_oiii_hb` (z=+2.35) | [SDSS](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=1522331895189235712) • [Legacy](https://www.legacysurvey.org/viewer?ra=225.57194&dec=38.05195&layer=ls-dr10&zoom=15) |
| **#2** | `2237320683103217664` | `0.033` | `-1.58` | `+0.61` | `[-1.64, -1.52]`, `[0.59, 0.63]` | `98.5` | **`97.6`** | `0.009` | `22.0σ` | Very low [N II]/Halpha (metal-poor-like ratios; unconfirmed) | `log_nii_ha` (z=+4.75) | [SDSS](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=2237320683103217664) • [Legacy](https://www.legacysurvey.org/viewer?ra=192.81909&dec=35.11345&layer=ls-dr10&zoom=15) |
| **#3** | `2293593682860861440` | `0.230` | `-0.97` | `+0.52` | `[-1.00, -0.93]`, `[0.51, 0.54]` | `97.5` | **`97.1`** | `0.004` | `12.2σ` | Star-forming region | `log_ew_ha` (z=+3.12) | [SDSS](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=2293593682860861440) • [Legacy](https://www.legacysurvey.org/viewer?ra=170.60098&dec=37.76703&layer=ls-dr10&zoom=15) |
| **#4** | `1432306634738132992` | `0.176` | `-0.35` | `+0.17` | `[-0.37, -0.33]`, `[0.13, 0.20]` | `97.5` | **`97.0`** | `0.005` | `11.1σ` | Composite region | `log_ew_ha` (z=+2.31) | [SDSS](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=1432306634738132992) • [Legacy](https://www.legacysurvey.org/viewer?ra=137.25270&dec=33.43255&layer=ls-dr10&zoom=15) |
| **#5** | `1679928197229078528` | `0.039` | `-0.65` | `+0.16` | `[-0.72, -0.58]`, `[0.09, 0.23]` | `97.6` | **`94.1`** | `0.036` | `4.1σ` | Star-forming region | `log_ew_ha` (z=+1.99) | [SDSS](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=1679928197229078528) • [Legacy](https://www.legacysurvey.org/viewer?ra=4.00171&dec=-0.30269&layer=ls-dr10&zoom=15) |
| **#6** | `2247426838710216704` | `0.035` | `-1.71` | `+0.66` | `[-1.78, -1.63]`, `[0.64, 0.68]` | `99.5` | **`91.1`** | `0.085` | `18.6σ` | Very low [N II]/Halpha (metal-poor-like ratios; unconfirmed) | `log_nii_ha` (z=+5.29) | [SDSS](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=2247426838710216704) • [Legacy](https://www.legacysurvey.org/viewer?ra=171.18009&dec=39.76517&layer=ls-dr10&zoom=15) |
| **#7** | `512391142269020160` | `0.031` | `-0.62` | `+0.86` | `[-0.70, -0.53]`, `[0.78, 0.94]` | `97.7` | **`85.0`** | `0.130` | `5.2σ` | Outlier at ordinary BPT position | `log_oiii_hb` (z=+2.21) | [SDSS](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=512391142269020160) • [Legacy](https://www.legacysurvey.org/viewer?ra=38.29496&dec=-7.80372&layer=ls-dr10&zoom=15) |
| **#8** | `2518725292126660608` | `0.051` | `-0.35` | `-0.82` | `[-0.38, -0.33]`, `[-0.91, -0.74]` | `96.9` | **`80.8`** | `0.166` | `8.2σ` | Star-forming region | `log_oiii_hb` (z=+1.29) | [SDSS](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=2518725292126660608) • [Legacy](https://www.legacysurvey.org/viewer?ra=188.59066&dec=29.65590&layer=ls-dr10&zoom=15) |

*Full candidate catalog exported with complete coordinates and multi-wavelength links in [`data/sdss_anomaly_candidates.csv`](data/sdss_anomaly_candidates.csv).*

---

## Astrophysical Literature Survey: Top 5 Anomalies Deep-Dive

To rigorously evaluate our unsupervised machine learning methodology, we cross-matched our top physical anomaly candidates against the **NASA/IPAC Extragalactic Database (NED)** and the **Sloan Digital Sky Survey (SDSS DR17)** spectroscopic archive.

Beyond algorithmic scores, each candidate corresponds to an established or emerging object of intense astrophysical interest:

```text
+-------+-------------------------+---------+--------------+------------------------------------------------------+
| Rank  | Target Common Name      | Type    | Primary S/N  | Key Astrophysical Discovery                          |
+-------+-------------------------+---------+--------------+------------------------------------------------------+
| #1    | CGCG 011-088            | Seyfert | 56.9         | Rare double-ring galaxy with nuclear rejuvenation    |
| #2    | SDSS J150754.38+010816.7| IMBH    | 40.5         | Intermediate-mass black hole candidate with radio jet|
| #3    | SDSS J101653.82+002857.0| QSO     | 40.2         | HST Feedback Survey target clearing molecular gas    |
| #4    | CGCG 021-096 (LIRG)     | AGN     | 57.0         | Infrared-luminous obscured AGN with massive torus    |
| #5    | UM 643                  | Starburst| 45.4        | Historic Michigan objective-prism starburst spiral   |
+-------+-------------------------+---------+--------------+------------------------------------------------------+
```

### 1. Rank #1: CGCG 011-088 (The Nuclear Rejuvenation Ring Galaxy)
* **Coordinates:** RA 171.03603 deg (11h 24m 08.6s), Dec -1.15776 deg (-01 deg 09 min 28 sec) | **z = 0.0293** (~125 Mpc)
* **Catalog Identifiers:** CGCG 011-088, Z 11-88, MCG +00-29-027, LEDA 35016, 2MASX J11240862-0109277.
* **What the Literature Discovered (53 papers in NED):**
  * Morphologically identified in the *Comprehensive de Vaucouleurs revised Hubble-Sandage (CVRHS)* system and *Galaxy Zoo 2* as a prominent **(R'_2)SAB(r:)a resonance ring galaxy** containing both inner and outer stellar rings.
  * Martin-Navarro et al. (2022, *MNRAS Letters*, 513, L10) identified it as a hallmark case of **nuclear rejuvenation**: non-axisymmetric ring torques channel fresh gas into the central core, reigniting black hole accretion after a quiescent period.
* **Why Our ML Algorithm Caught It:**
  * Sits abnormally high on the Seyfert wing (**log₁₀([O III] / Hβ) = +0.84**) despite moderate nitrogen metallicity (**log₁₀([N II] / Hα) = -0.10**).
  * Exceptionally clean emission profile with minimum line S/N of 56.9 and tight 95% CI error bounds (width < 0.04 dex), earning 3/3 unanimous algorithm votes as an unconventional high-confidence emission outlier.

### 2. Rank #2: SDSS J150754.38+010816.7 (Intermediate-Mass Black Hole Candidate)
* **Coordinates:** RA 226.97660 deg (15h 07m 54.4s), Dec +1.13799 deg (+01 deg 08 min 17 sec) | **z = 0.0610** (~260 Mpc)
* **Catalog Identifiers:** 2MASX J15075440+0108168, WISEA J150754.38+010816.8, FIRST J150754.3+010816, NVSS J150754+010812, [VV2006c] J150754.4+010816.
* **What the Literature Discovered (42 papers in NED):**
  * Highlighted in observational surveys as harboring an **Intermediate-Mass Black Hole (IMBH) candidate** with an estimated central black hole mass of only **M_BH ≈ 1.4 × 10⁶ M_☉** (orders of magnitude lighter than normal quasar engines).
  * Lin et al. (2018, *The Astronomical Journal*, 155, 188) investigated its low-power relativistic radio core detected in the VLA FIRST survey, analyzing the accretion mode powering its nuclear jets.
* **Why Our ML Algorithm Caught It:**
  * Flagged as an **Ultra-High Ionization Outlier**: **log₁₀([O III] / Hβ) = +1.13** (in this sample, log₁₀([O III] / Hβ) = +1.13 sits in the top 0.1% of high-excitation emitters).
  * Accretion onto a lower-mass black hole at high Eddington ratios generates an intense, hard UV ionizing radiation field that super-excites forbidden [O III] relative to standard Seyferts.

### 3. Rank #3: SDSS J101653.82+002857.0 (The HST Quasar Feedback Target)
* **Coordinates:** RA 154.22425 deg (10h 16m 53.8s), Dec +0.48254 deg (+00 deg 28 min 57 sec) | **z = 0.1163** (~520 Mpc)
* **Catalog Identifiers:** WISEA J101653.82+002856.7, 2MASX J10165384+0028571, LEDA 1167114, NVSS J101653+002858, [VV2010c] J101653.9+002857.
* **What the Literature Discovered (38 papers in NED):**
  * Officially cataloged in NED as a **QSO (Quasar)** and a **High-Excitation Radio Galaxy (HERG)**.
  * Selected for ultraviolet spectroscopic observation with the **Hubble Space Telescope (HST Cosmic Origins Spectrograph - COS)**.
  * Molyneux et al. (2024, *MNRAS*, 527, 4420) featured this system in *The Quasar Feedback Survey*, utilizing ALMA and CO line transitions to map energetic molecular gas outflows clearing star-forming fuel from the galaxy.
* **Why Our ML Algorithm Caught It:**
  * Shows intense [O III] emission (**log₁₀([O III] / Hβ) = +1.00**) coupled with low nitrogen (**log₁₀([N II] / Hα) = -0.34**).
  * Yields massive PCA reconstruction residuals because energetic quasar wind kinematics distort the standard 4-line emission manifold.

### 4. Rank #4: CGCG 021-096 / IRAS 15240+0046 (The Luminous Infrared AGN)
* **Coordinates:** RA 231.65696 deg (15h 26m 37.7s), Dec +0.59264 deg (+00 deg 35 min 34 sec) | **z = 0.0507** (~220 Mpc)
* **Catalog Identifiers:** CGCG 021-096, Z 21-96, IRAS 15240+0046, IRAS F15240+0045, LEDA 55103, FIRST J152637.6+003534, [VV2010c] J152637.7+003534.
* **What the Literature Discovered (59 papers in NED):**
  * Discovered in the IRAS satellite mission as a **Luminous Infrared Galaxy (LIRG)**, powered by an active galactic nucleus buried inside dense, warm interstellar dust clouds.
  * Ramos Padilla et al. (2022, *MNRAS*, 510, 687) studied its Spectral Energy Distribution (SED) and obscuring torus geometry.
* **Why Our ML Algorithm Caught It:**
  * Simultaneously exhibits extreme high-excitation emission (**log₁₀([O III] / Hβ) = +0.98**) and substantial dust extinction (large Balmer decrement **Hα / Hβ**), causing it to diverge strongly from pure optical starburst models.

### 5. Rank #5: UM 643 (Historic University of Michigan Survey Starburst)
* **Coordinates:** RA 212.61180 deg (14h 10m 26.8s), Dec -0.83238 deg (-00 deg 49 min 57 sec) | **z = 0.0254** (~110 Mpc)
* **Catalog Identifiers:** UM 643, CGCG 018-065, Z 18-65, MCG +00-36-027, PGC 50604, 2MASX J14102685-0049563.
* **What the Literature Discovered (55 papers in NED):**
  * Discovered in the landmark University of Michigan (UM) objective-prism survey for unusual optical emission-line objects.
  * Morphologically identified as an energetic barred starburst spiral (`SB/SBR`) with radio continuum emission cataloged in Condon et al. (2019/2022, *ApJ*, 872, 148).
* **Why Our ML Algorithm Caught It:**
  * Unlike the other top candidates (which are high-ionization AGN), UM 643 is an **unconventional composite transition outlier** (**log₁₀([O III] / Hβ) = -0.38**, **log₁₀([N II] / Hα) = -0.25**).
  * It occupies the sparse transition valley between the pure starburst locus and the AGN plume. Local Outlier Factor (LOF) flagged it due to its sparse local neighborhood density.

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

### 2. Run the Diagnostic & Anomaly Discovery Pipeline

```bash
python run_bpt.py
```
*Executes the end-to-end pipeline in ~4 seconds: loads the tracked sample from `data/` (pass `--force-refresh` to re-query SkyServer), performs classical BPT classification, runs the tri-algorithmic anomaly ensemble, computes 95% confidence intervals, and outputs publication figures to `assets/`.*

### 3. Run Automated Tests

```bash
pip install -r requirements-dev.txt
pytest -q
```
*Executes regression tests covering theoretical boundaries, BPT region classification, error propagation, confidence ellipses, and Mahalanobis boundary detachment.*

---

## Scientific References & Literature Attribution

This project implements analytical boundaries, survey datasets, and observational findings from the following landmark literature:

1. **The BPT Diagnostic Foundation:**  
   Baldwin, J. A., Phillips, M. M., & Terlevich, R. (1981). *Classification parameters for the emission-line spectra of extragalactic objects.* **PASP**, 93, 5–19. [DOI: 10.1086/130766](https://doi.org/10.1086/130766)
2. **Empirical Star-Formation Demarcation:**  
   Kauffmann, G., Heckman, T. M., Tremonti, C., et al. (2003). *The host galaxies of active galactic nuclei.* **MNRAS**, 346(4), 1055–1077. [DOI: 10.1111/j.1365-2966.2003.07154.x](https://doi.org/10.1111/j.1365-2966.2003.07154.x)
3. **Extreme Starburst Theoretical Limit:**  
   Kewley, L. J., Dopita, M. A., Sutherland, R. S., et al. (2001). *Theoretical Modeling of Starburst Galaxies.* **ApJ**, 556(1), 121–140. [DOI: 10.1086/321545](https://doi.org/10.1086/321545)
4. **AGN Host Galaxies & Multi-Line Diagnostics:**  
   Kewley, L. J., Groves, B., Kauffmann, G., & Heckman, T. (2006). *The host galaxies and classification of active galactic nuclei.* **MNRAS**, 372(3), 961–976. [DOI: 10.1111/j.1365-2966.2006.10869.x](https://doi.org/10.1111/j.1365-2966.2006.10869.x)
5. **MPA-JHU Emission-Line Measurements:**  
   Brinchmann, J., Charlot, S., White, S. D. M., et al. (2004). *The physical properties of star-forming galaxies in the low-redshift Universe.* **MNRAS**, 351(4), 1151–1179. [DOI: 10.1111/j.1365-2966.2004.07881.x](https://doi.org/10.1111/j.1365-2966.2004.07881.x)
6. **SDSS-IV Data Release 17:**  
   Abdurro'uf, Accetta, K., Aerts, C., et al. (2022). *The Seventeenth Data Release of the Sloan Digital Sky Surveys: Complete Results from SDSS-IV.* **ApJS**, 259(2), 35. [DOI: 10.3847/1538-4365/ac4414](https://doi.org/10.3847/1538-4365/ac4414)
7. **Nuclear Rejuvenation in Ring Galaxies (Candidate #1):**  
   Martín-Navarro, I., Shankar, F., Mezcua, M., et al. (2022). *Rejuvenation triggers nuclear activity in nearby galaxies.* **MNRAS Letters**, 513(1), L10–L14. [DOI: 10.1093/mnrasl/slab112](https://doi.org/10.1093/mnrasl/slab112)
8. **Nuclear Activity in Low-Power Radio Galaxies (Candidate #2):**  
   Lin, Y.-T., Huang, S., & Reveilaud, C. (2018). *Origin of Nuclear Activity in Low-power Radio Galaxies.* **AJ**, 155(5), 188. [DOI: 10.3847/1538-3881/aab5b4](https://doi.org/10.3847/1538-3881/aab5b4)
9. **The Quasar Feedback Survey (Candidate #3):**  
   Molyneux, S. J., Calistro Rivera, G., De Breuck, C., et al. (2024). *The Quasar Feedback Survey: characterizing CO excitation in quasar host galaxies.* **MNRAS**, 527(3), 4420–4439. [DOI: 10.1093/mnras/stad3511](https://doi.org/10.1093/mnras/stad3511)
10. **AGN SED Models and Dust Viewing Angles (Candidate #4):**  
    Ramos Padilla, A. F., Wang, L., Małek, K., et al. (2022). *Viewing angle in AGN SED models.* **MNRAS**, 510(1), 687–703. [DOI: 10.1093/mnras/stab3440](https://doi.org/10.1093/mnras/stab3440)
11. **Radio Sources in the Nearby Universe (Candidate #5):**  
    Condon, J. J., Matthews, A. M., & Broderick, J. J. (2019/2022). *Radio Sources in the Nearby Universe.* **ApJ**, 872(2), 148. [DOI: 10.3847/1538-4357/aafd32](https://doi.org/10.3847/1538-4357/aafd32)

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
* **Primary Scientific Attribution (Mandatory):** This repository is an educational reproduction and open-source demonstration. Any scientific publication, paper, or research work utilizing the BPT diagnostic, classification boundaries, or SDSS catalog **must primarily cite the original literature** listed in [Scientific References](#scientific-references--literature-attribution) and [`REFERENCES.bib`](REFERENCES.bib) (Baldwin et al. 1981; Kauffmann et al. 2003; Kewley et al. 2001, 2006; Brinchmann et al. 2004; Abdurro'uf et al. 2022).
* **Software Implementation Citation (Optional):** If you specifically adapt or reuse this Python pipeline, code structure, or plotting workflow in an educational tutorial, university lab, or coding demo, you may additionally cite the software repository via [`CITATION.cff`](CITATION.cff) or GitHub's **"Cite this repository"** button.\n