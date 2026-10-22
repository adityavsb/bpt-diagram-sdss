# Mapping the Cosmic Boundary: Starburst Galaxies vs. Accreting Black Holes (SDSS DR17)

[![Python](https://img.shields.io/badge/Python-3.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://www.python.org/)
[![Astropy](https://img.shields.io/badge/Astropy-6.x%2B-orange.svg)](https://www.astropy.org/)
[![Astroquery](https://img.shields.io/badge/Astroquery-SDSS%20DR17-purple.svg)](https://astroquery.readthedocs.io/)
[![Survey](https://img.shields.io/badge/Survey-SDSS%20DR17%20%2F%20MPA--JHU-teal.svg)](https://www.sdss.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An observational astrophysics project implementing the landmark **Baldwin, Phillips & Terlevich (BPT, 1981)** emission-line diagnostic on ~10,000 galaxies from the **Sloan Digital Sky Survey (SDSS DR17)**.

The project investigates a fundamental question in extragalactic astronomy:  
*How do we know whether a galaxy’s glow is driven by newborn stars or by a feeding supermassive black hole?*

---

## The Astrophysical Principle: Reading the Gas

When high-energy radiation strikes interstellar gas, atoms ionize and recombine, emitting photons at discrete spectral wavelengths. By taking the logarithmic ratio of specific forbidden collision lines to hydrogen recombination lines, we construct an excitation diagnostic that is independent of distance or dust reddening:

* **$y$-axis:** $\log_{10}([\mathrm{O\ III}]\,\lambda 5007 \,/\, \mathrm{H}\beta)$ — tracks the ionization parameter and hardness of the ionizing radiation field.
* **$x$-axis:** $\log_{10}([\mathrm{N\ II}]\,\lambda 6584 \,/\, \mathrm{H}\alpha)$ — tracks gas-phase metallicity and nitrogen abundance.

```text
               High [O III] / Hβ  (Hard Ionizing Photons)
                         ▲
                         │       ▲
                         │      /   AGN Wing (Seyferts)
                         │     /    [Hard power-law accretion]
                         │    /
                         │   /      ▲
                         │  /      /   LINERs (Shocks / Old Stars)
       Star-Forming      │ /      /
       (OB Star Nurseries)
                         │
                         ▼
        Low [N II] / Hα ◄─────────────────► High [N II] / Hα
                                          (High Metallicity)
```

---

## Diagnostic Boundaries

We implement the standard theoretical and empirical classification curves from the literature:

1. **Kauffmann et al. (2003) Pure Star-Formation Limit:**
   $$\log_{10}\left(\frac{[\mathrm{O\ III}]}{\mathrm{H}\beta}\right) = \frac{0.61}{\log_{10}([\mathrm{N\ II}]/\mathrm{H}\alpha) - 0.05} + 1.30 \quad (\text{for } x < 0.05)$$
   Galaxies below this dashed curve are powered solely by young massive OB stars ($< 10\text{ Myr}$).

2. **Kewley et al. (2001) Extreme Starburst Theoretical Limit:**
   $$\log_{10}\left(\frac{[\mathrm{O\ III}]}{\mathrm{H}\beta}\right) = \frac{0.61}{\log_{10}([\mathrm{N\ II}]/\mathrm{H}\alpha) - 0.47} + 1.19 \quad (\text{for } x < 0.47)$$
   Galaxies above this solid curve cannot be ionized by starbursts alone, requiring an Active Galactic Nucleus (AGN).

3. **Kewley et al. (2006) Seyfert vs. LINER Demarcation:**
   $$\log_{10}\left(\frac{[\mathrm{O\ III}]}{\mathrm{H}\beta}\right) = 1.89 \log_{10}\left(\frac{[\mathrm{N\ II}]}{\mathrm{H}\alpha}\right) + 0.76$$
   Separates high-ionization Seyferts (active black holes) from low-ionization LINERs.

---

## Observational Results (9,999 SDSS DR17 Galaxies)

Filtered for validated spectroscopic galaxies ($0.02 < z < 0.25$) with Signal-to-Noise $S/N \ge 3$ across all 4 diagnostic emission lines:

| Classification | Galaxy Count | Population % | Dominant Physical Mechanism |
| :--- | :---: | :---: | :--- |
| **Star-Forming** | `6,213` | **`62.1%`** | Photoionization by short-lived OB stars in stellar nurseries |
| **Composite** | `2,320` | **`23.2%`** | Transition systems undergoing simultaneous starbursts and AGN accretion |
| **Seyfert (AGN)** | `226` | **`2.3%`** | Hard UV/X-ray ionizing radiation from supermassive black hole accretion disks |
| **LINER** | `1,240` | **`12.4%`** | Low-ionization shock excitation or evolved post-AGB stellar populations |

![BPT Diagnostic Diagram](assets/bpt_diagnostic_sdss.png)

---

## Two Anchor Galaxies: M51 vs. NGC 1068

To anchor the catalog statistics to real objects you can see in the sky:

| Anchor Galaxy | Type | $\log([\mathrm{N\ II}]/\mathrm{H}\alpha)$ | $\log([\mathrm{O\ III}]/\mathrm{H}\beta)$ | Diagnostic Position |
| :--- | :--- | :---: | :---: | :--- |
| **M51 (Whirlpool Galaxy)** | Star-Forming Spiral | `-0.42` | `-0.38` | Deep in the lower-left starburst locus. Its glow is driven by massive young stars. |
| **NGC 1068 (M77)** | Seyfert 2 Active Nucleus | `+0.15` | `+0.95` | Prominently in the upper-right AGN wing. Its emission lines require an accreting supermassive black hole. |

---

## Project Structure

```text
bpt-diagram-sdss/
├── README.md               # Scientific documentation & observational results
├── PLAN.md                 # Technical architecture & astrophysics plan
├── CITATION.cff            # Citation metadata for academic and open-source citing
├── REFERENCES.bib          # NASA/ADS BibTeX database for LaTeX / paper references
├── LICENSE                 # MIT License
├── .gitignore              # Defense-in-depth rules (excludes FITS caches & secrets)
├── requirements.txt        # astropy, astroquery, pandas, matplotlib, numpy
├── src/
│   ├── bpt_math.py         # Boundary equations (Kauffmann 2003, Kewley 2001/2006)
│   ├── query_data.py       # SDSS DR17 query & local caching with quality cuts
│   ├── plot_bpt.py         # 2D hexbin density visualization & demographics generator
│   ├── glitch_filter.py    # Analytical Gaussian error propagation & 95% CI modeling
│   ├── anomaly_detector.py # Tri-algorithmic ensemble (iForest + LOF + PCA) & voting
│   └── plot_anomalies.py   # Publication-quality anomaly diagnostic map
├── assets/
│   ├── bpt_diagnostic_sdss.png # Classical BPT diagram
│   └── bpt_anomalies_sdss.png  # Anomaly detection & 95% CI error ellipses
├── data/
│   ├── sdss_bpt_sample.csv          # Cached 9,999 SDSS DR17 galaxies
│   └── sdss_anomaly_candidates.csv  # Ranked Top 20 physical celestial candidates
└── run_bpt.py              # Master pipeline runner (Classical BPT + Anomaly Discovery)
```

---

## Iteration 2: Tri-Algorithmic Celestial Anomaly Detection & Glitch Disentanglement

While classical demarcation curves partition galaxies into deterministic boxes (*Star-Forming, Composite, Seyfert, LINER*), real-world survey astronomy makes breakthrough discoveries (e.g. **Green Pea galaxies, dual supermassive black holes, intermediate-mass black holes, tidal disruption events**) precisely where objects break these standard rules.

Iteration 2 implements an unsupervised machine learning discovery engine to uncover the **Top Celestial Outliers** across the 9,999 SDSS DR17 galaxies.

![SDSS DR17 BPT Anomaly Detection](assets/bpt_anomalies_sdss.png)

### 1. The Tri-Algorithmic Voting Ensemble

To eliminate single-algorithm bias, we employ an equal-weighted consensus ensemble across three complementary unsupervised paradigms:
1. **Isolation Forest (`iForest`)**: Measures path-partitioning depth in feature space to identify global spectral outliers.
2. **Local Outlier Factor (`LOF`, $k=35$)**: Density-ratio metric that uncovers objects lingering in the low-density "composite desert" between the starburst locus and the AGN plume.
3. **PCA Reconstruction Residuals (`PCA-Rec`)**: Projects the multi-line emission space onto principal components; galaxies with high $L_2$ reconstruction error exhibit unusual ionization physics unexplained by standard stellar/AGN models.

Each raw score is calibrated to a percentile scale ($0\text{--}100$). The **Consensus Anomaly Score** is their equal-weighted average:
$$S_{\text{consensus}} = \frac{1}{3} S_{\text{iForest}} + \frac{1}{3} S_{\text{LOF}} + \frac{1}{3} S_{\text{PCA}}$$

### 2. Disentangling Physical Anomalies vs. Pipeline Glitches

In survey astronomy, the raw top outliers of any unsupervised algorithm are often **data artifacts** (e.g., cosmic ray hits, sky subtraction residuals, or faint-line flux denominators approaching zero).

To ensure candidates are scientifically genuine, we implement:
* **Analytical Gaussian Error Propagation:**
  $$\sigma_x = \frac{1}{\ln 10} \sqrt{\left(\frac{\sigma_{[\mathrm{N\ II}]}}{[\mathrm{N\ II}]}\right)^2 + \left(\frac{\sigma_{\mathrm{H}\alpha}}{\mathrm{H}\alpha}\right)^2}, \quad \sigma_y = \frac{1}{\ln 10} \sqrt{\left(\frac{\sigma_{[\mathrm{O\ III}]}}{[\mathrm{O\ III}]}\right)^2 + \left(\frac{\sigma_{\mathrm{H}\beta}}{\mathrm{H}\beta}\right)^2}$$
* **$95\%$ Confidence Interval Ellipses:** $\text{CI}_{95\%} = [x \pm 1.96\,\sigma_x, \; y \pm 1.96\,\sigma_y]$.
* **`pipeline_glitch_risk` ($0\text{--}100\%$):** Flags objects whose outlier status is fragile to measurement noise ($S/N \to 3.0$ or huge error ellipse area).
* **`physical_anomaly_confidence` ($0\text{--}100\%$):** Confirms high ML consensus backed by high signal-to-noise ($S/N \gg 10$), compact error ellipses, and statistical detachment ($>3\sigma$) from standard tracks.

### 3. Top Confirmed Physical Celestial Candidates

Filtered for validated low glitch risk (`pipeline_glitch_risk` $\le 0.5\%$, $S/N_{\min} \ge 40$), the top candidates exhibit undeniable physical oddities:

| Rank | SpecObjID | $z$ | $\log([\mathrm{N\ II}]/\mathrm{H}\alpha)$ | $\log([\mathrm{O\ III}]/\mathrm{H}\beta)$ | $95\%$ CI Ellipse | Consensus Score | Phys. Conf. | Physical Archetype | Optical Sky Links |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| **#1** | `315265744856180736` | `0.029` | `-0.10` | `+0.84` | `[-0.11, -0.09]`, `[+0.82, +0.86]` | `99.2` | **`99.0%`** | Unconventional Emission Outlier | [SDSS Explore](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=315265744856180736) • [DESI Legacy](https://www.legacysurvey.org/viewer?ra=171.03603&dec=-1.15776&layer=ls-dr10&zoom=15) |
| **#2** | `350255779153995776` | `0.061` | `-0.17` | `+1.13` | `[-0.18, -0.16]`, `[+1.11, +1.14]` | `99.4` | **`99.0%`** | Ultra-High Ionization Outlier | [SDSS Explore](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=350255779153995776) • [DESI Legacy](https://www.legacysurvey.org/viewer?ra=226.97660&dec=1.13799&layer=ls-dr10&zoom=15) |
| **#3** | `305239577747023872` | `0.116` | `-0.34` | `+1.00` | `[-0.35, -0.33]`, `[+0.99, +1.02]` | `99.2` | **`98.9%`** | Ultra-High Ionization Outlier | [SDSS Explore](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=305239577747023872) • [DESI Legacy](https://www.legacysurvey.org/viewer?ra=154.22425&dec=0.48254&layer=ls-dr10&zoom=15) |
| **#4** | `352546336886777856` | `0.051` | `-0.10` | `+0.98` | `[-0.11, -0.09]`, `[+0.96, +1.00]` | `99.1` | **`98.8%`** | Ultra-High Ionization Outlier | [SDSS Explore](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=352546336886777856) • [DESI Legacy](https://www.legacysurvey.org/viewer?ra=231.65696&dec=0.59264&layer=ls-dr10&zoom=15) |
| **#9** | `301908877954279424` | `0.032` | `-0.23` | `+0.20` | `[-0.24, -0.22]`, `[+0.18, +0.22]` | `98.1` | **`98.0%`** | Transition Composite / Void Outlier | [SDSS Explore](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=301908877954279424) • [DESI Legacy](https://www.legacysurvey.org/viewer?ra=146.90798&dec=-0.62646&layer=ls-dr10&zoom=15) |
| **#10** | `341263147589462016` | `0.053` | `+0.01` | `+0.93` | `[0.00, +0.02]`, `[+0.91, +0.95]` | `98.2` | **`97.8%`** | Extreme Seyfert / Fast Accretion Outflow | [SDSS Explore](http://skyserver.sdss.org/dr17/en/tools/explore/summary.aspx?sid=341263147589462016) • [DESI Legacy](https://www.legacysurvey.org/viewer?ra=216.59478&dec=-0.56942&layer=ls-dr10&zoom=15) |

*Full 20-candidate catalog exported with complete coordinates and multi-wavelength links in [`data/sdss_anomaly_candidates.csv`](data/sdss_anomaly_candidates.csv).*


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

### 2. Run the Diagnostic Pipeline

```bash
python run_bpt.py
```
*Queries SDSS DR17 (or loads the cached validated catalog), calculates emission-line ratios, applies demarcation curves, and outputs the publication-ready diagnostic plot to `assets/bpt_diagnostic_sdss.png`.*

---

## Scientific References & Literature Attribution

This project directly implements analytical boundaries and datasets established in the following landmark papers:

1. **The BPT Diagnostic Foundation:**  
   Baldwin, J. A., Phillips, M. M., & Terlevich, R. (1981). *Classification parameters for the emission-line spectra of extragalactic objects.* **PASP**, 93, 5–19. [DOI: 10.1086/130766](https://doi.org/10.1086/130766)
2. **Empirical Star-Formation Demarcation:**  
   Kauffmann, G., Heckman, T. M., Tremonti, C., et al. (2003). *The host galaxies of active galactic nuclei.* **MNRAS**, 346(4), 1055–1077. [DOI: 10.1111/j.1365-2966.2003.07154.x](https://doi.org/10.1111/j.1365-2966.2003.07154.x)
3. **Extreme Starburst Theoretical Limit:**  
   Kewley, L. J., Dopita, M. A., Sutherland, R. S., et al. (2001). *Theoretical Modeling of Starburst Galaxies.* **ApJ**, 556(1), 121–140. [DOI: 10.1086/321545](https://doi.org/10.1086/321545)
4. **Seyfert vs. LINER Optical Classification:**  
   Kewley, L. J., Groves, B., Kauffmann, G., & Heckman, T. (2006). *The host galaxies and classification of active galactic nuclei.* **MNRAS**, 372(3), 961–976. [DOI: 10.1111/j.1365-2966.2006.10869.x](https://doi.org/10.1111/j.1365-2966.2006.10869.x)
5. **MPA-JHU Emission-Line Measurements:**  
   Brinchmann, J., Charlot, S., White, S. D. M., et al. (2004). *The physical properties of star-forming galaxies in the low-redshift Universe.* **MNRAS**, 351(4), 1151–1179. [DOI: 10.1111/j.1365-2966.2004.07881.x](https://doi.org/10.1111/j.1365-2966.2004.07881.x)
6. **SDSS-IV Data Release 17:**  
   Abdurro'uf, Accetta, K., Aerts, C., et al. (2022). *The Seventeenth Data Release of the Sloan Digital Sky Surveys: Complete Results from SDSS-IV.* **ApJS**, 259(2), 35. [DOI: 10.3847/1538-4365/ac4414](https://doi.org/10.3847/1538-4365/ac4414)

Full BibTeX records are available in [`REFERENCES.bib`](REFERENCES.bib).

---

## Data Rights & SDSS Acknowledgement

This project uses public data from the **Sloan Digital Sky Survey IV (SDSS-IV)**. 

> *Funding for the Sloan Digital Sky Survey IV has been provided by the Alfred P. Sloan Foundation, the U.S. Department of Energy Office of Science, and the Participating Institutions. SDSS-IV acknowledges support and resources from the Center for High-Performance Computing at the University of Utah. The SDSS web site is [www.sdss.org](https://www.sdss.org).*
>
> *SDSS-IV is managed by the Astrophysical Research Consortium for the Participating Institutions of the SDSS Collaboration.*

---

## License & Academic Attribution Policy

* **Code License:** [MIT License](LICENSE) © 2026 Aditya Bhadoria.
* **Primary Scientific Attribution (Mandatory):** This repository is an educational reproduction and open-source demonstration. Any scientific publication, paper, or research work utilizing the BPT diagnostic, classification boundaries, or SDSS catalog **must primarily cite the original literature** listed in [Scientific References](#scientific-references--literature-attribution) and [`REFERENCES.bib`](REFERENCES.bib) (Baldwin et al. 1981; Kauffmann et al. 2003; Kewley et al. 2001, 2006; Brinchmann et al. 2004; Abdurro'uf et al. 2022).
* **Software Implementation Citation (Optional):** If you specifically adapt or reuse this Python pipeline, code structure, or plotting workflow in an educational tutorial, university lab, or coding demo, you may additionally cite the software repository via [`CITATION.cff`](CITATION.cff) or GitHub's **"Cite this repository"** button.
