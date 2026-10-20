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
├── LICENSE                 # MIT License
├── .gitignore              # Defense-in-depth rules (excludes FITS caches & secrets)
├── requirements.txt        # astropy, astroquery, pandas, matplotlib, numpy
├── src/
│   ├── bpt_math.py         # Boundary equations (Kauffmann 2003, Kewley 2001/2006)
│   ├── query_data.py       # SDSS DR17 query & local caching with quality cuts
│   └── plot_bpt.py         # 2D hexbin density visualization & demographics generator
├── assets/
│   └── bpt_diagnostic_sdss.png
└── run_bpt.py              # Master pipeline runner
```

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

## License

This project is open-source and available under the [MIT License](LICENSE).
