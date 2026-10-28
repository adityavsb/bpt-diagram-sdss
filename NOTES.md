# Working notes

Dated, first-person record of decisions and corrections. Newest first.
These are notes, not documentation; see README for the current method.

## 2026-10 — Review and rework

### Decisions
* 2026-10-28 — Adopted 4D scale-invariant feature space: dropped raw Hα and Hβ line fluxes and moved to 4 pure ratios/widths (log([N II]/Hα), log([O III]/Hβ), Balmer decrement Hα/Hβ, log EW(Hα)). Eliminates Malmquist brightness bias and catalog scaling glitches (metrics_log.csv: row 8).
* 2026-10-28 — Dynamic PCA variance threshold (≥ 90%): automatically select minimum components (3 components, 92.8% variance) rather than fixing components arbitrarily.
* 2026-10-28 — 2-D 95% error ellipses with k = 2.448: corrected error ellipse semi-axes from 1.96σ (which covers only 85.3% in 2D) to 2.448σ (chi-square 2 d.o.f. for 95% joint confidence).
* 2026-10-28 — Deterministic hash sampling: ordered 10,000 SDSS galaxies by SHA-256 hash of specObjID to guarantee identical, reproducible cross-platform sampling without relying on arbitrary `TOP 10000` catalog ordering.
* 2026-10-28 — Automated SIMBAD cross-match with Harvard TAP mirror: replaced manual NED notes with reproducible cone queries via astroquery.simbad, adding automated fallback to the Harvard mirror when the Strasbourg CDS mirror experiences latency.

### Things I got wrong
* Malmquist bias masquerading as anomaly detection: In Iteration 1/2, I fed raw line fluxes into the unsupervised models. The algorithms voted for extreme brightness outliers (median Hα flux ~7.67 billion) rather than unusual line ratios.
* Confirmation bias from literature searches: Seeing dozens of papers in NED/SIMBAD for candidates like UM 643 gave me false confidence that the detector was discovering rare physics. In truth, those galaxies had papers because they were historically bright targets; their BPT line ratios were completely unremarkable and sat right in the center of the normal composite cloud.
* 1-D vs. 2-D error ellipse geometry: I initially used 1.96σ for 2-D error ellipses, forgetting that a bivariate Gaussian requires 2.448σ for 95% joint coverage. The 1.96σ box left ~14.7% of genuine normal galaxies outside the ellipse.
* M51 nucleus misidentification: In early documentation, I conflated the companion galaxy NGC 5195 with the primary nucleus of NGC 5194 (M51), misattributing the nuclear line ratios.
* Unimplemented S/N ≥ 40 filter: The README previously claimed a strict S/N ≥ 40 cut, but the code actually applied S/N ≥ 3 on unscaled flux errors. Now cleanly unified to S/N ≥ 3 with recommended MPA-JHU error rescaling factors.

### Candidates I stopped believing in
* UM 643 (Z 18-65): Flagged in the baseline 6D detector as an exotic transition valley starburst. In the 4D ratio space, it sits squarely inside the ordinary composite locus; it was selected in the baseline solely because of its extreme uncalibrated brightness.
* CGCG 011-088 (Z 11-88): Dropped from the top candidate tier because its baseline outlier score was driven by extreme flux scale rather than anomalous gas excitation.
* Distant high-redshift candidates without morphological resolution: Several candidates at z > 0.15 in early runs were artifacts of 3-arcsecond fiber aperture mixing host galaxy disk light with nuclear emission.

### Open questions
* Fixed-fiber aperture bias across redshift: In SDSS, the 3-arcsecond fiber covers ~1.8 kpc at z = 0.03, but over ~6 kpc at z = 0.12, inevitably mixing disk star formation into nuclear spectra at higher redshifts.
* Composite/LINER shock degeneracy: Whether low-luminosity outliers at the boundary of the composite and LINER regions are powered by low-velocity shocks from galaxy interactions or by radiatively inefficient accretion flows (RIAFs).

## 2026-10 — Iteration 2 (anomaly detection)
* Original idea: use unsupervised machine learning (Isolation Forest, LOF, PCA) to discover galaxies whose gas excitation departs from standard photoionization equilibrium.
* The first candidate list was dominated by extreme brightness outliers due to including raw fluxes; top candidates clustered at high flux scales and local distances.

## 2026-10 — Iteration 1 (classical BPT)
* Project origin: implement classical BPT diagnostic curves on modern SDSS DR17 data.
* Standardized to the 3-class standard (Star-Forming, Composite, AGN) on the [N II]/Hα diagram per Kewley et al. (2006), reserving Seyfert/LINER separation for 4-line diagnostic sets.
