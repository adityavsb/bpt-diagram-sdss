"""
Mathematical boundaries and classification logic for the BPT (Baldwin, Phillips & Terlevich 1981) Diagram.
Implements:
1. Kauffmann et al. (2003) empirical star-formation demarcation
2. Kewley et al. (2001) theoretical extreme starburst upper limit
3. Kewley et al. (2006) Seyfert vs. LINER division
"""
import numpy as np

def kauffmann_boundary(log_nii_ha):
    """
    Kauffmann et al. (2003, MNRAS, 346, 1055, Eq. 1) empirical star-formation limit:
        log([O III] / Hβ) = 0.61 / (log([N II] / Hα) - 0.05) + 1.30
    
    Valid for log([N II] / Hα) < 0.05. Galaxies below this line are purely
    star-forming (photoionized by massive OB stars).
    Returns NaN for x >= 0.05.
    """
    x = np.asarray(log_nii_ha)
    y = np.full_like(x, np.nan, dtype=np.float64)
    valid = x < 0.05
    y[valid] = (0.61 / (x[valid] - 0.05)) + 1.30
    return y

def kewley_boundary(log_nii_ha):
    """
    Kewley et al. (2001, ApJ, 556, 121, Eq. 5) extreme starburst theoretical limit:
        log([O III] / Hβ) = 0.61 / (log([N II] / Hα) - 0.47) + 1.19
    
    Valid for log([N II] / Hα) < 0.47. Represents the maximum energy contribution
    achievable by stellar photoionization alone without black hole accretion or shocks.
    Returns NaN for x >= 0.47.
    """
    x = np.asarray(log_nii_ha)
    y = np.full_like(x, np.nan, dtype=np.float64)
    valid = x < 0.47
    y[valid] = (0.61 / (x[valid] - 0.47)) + 1.19
    return y

def seyfert_liner_boundary(log_nii_ha):
    """
    Kewley et al. (2006, MNRAS, 372, 961, Eq. 1) demarcation separating
    high-ionization Seyfert nuclei from low-ionization LINERs:
        log([O III] / Hβ) = 1.89 * log([N II] / Hα) + 0.76
    """
    x = np.asarray(log_nii_ha)
    return 1.89 * x + 0.76

def classify_bpt(log_nii_ha, log_oiii_hb):
    """
    Classifies a galaxy into one of 4 physical categories:
    - 'Star-Forming': Pure photoionization by young OB stars
    - 'Composite': Mixture of starburst and black hole accretion
    - 'Seyfert': Hard power-law ionization from supermassive black hole
    - 'LINER': Low-ionization nuclear emission-line region (shocks / old stars)
    """
    kauff_y = kauffmann_boundary(log_nii_ha)
    kewley_y = kewley_boundary(log_nii_ha)
    liner_y = seyfert_liner_boundary(log_nii_ha)

    is_star_forming = (log_nii_ha < 0.05) & (log_oiii_hb < kauff_y)
    is_composite = (log_nii_ha < 0.47) & (log_oiii_hb < kewley_y) & ~is_star_forming
    is_agn = ~is_star_forming & ~is_composite

    is_seyfert = is_agn & (log_oiii_hb >= liner_y)
    is_liner = is_agn & (log_oiii_hb < liner_y)

    categories = np.empty(len(log_nii_ha), dtype=object)
    categories[is_star_forming] = "Star-Forming"
    categories[is_composite] = "Composite"
    categories[is_seyfert] = "Seyfert (AGN)"
    categories[is_liner] = "LINER"
    return categories
