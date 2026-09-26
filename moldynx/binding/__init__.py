"""
Binding free energy by MM-GBSA / MM-PBSA with gmx_MMPBSA.

MolDynX prepares, gates, runs and analyses calculations performed with
**gmx_MMPBSA** (https://github.com/Valdes-Tresanca-MS/gmx_MMPBSA; Valdés-Tresanca
et al., J. Chem. Theory Comput. 2021, 17, 6281–6291), which drives AmberTools'
MMPBSA.py (Miller et al., J. Chem. Theory Comput. 2012, 8, 3314–3321). Cite both
when you use these results.

* :mod:`moldynx.binding.prepare` -- protein-only complex system, input files,
  derived ionic strength, explicit decomposition residue set, run script with a
  probe gate.
* :mod:`moldynx.binding.analyse` -- parse the outputs, autocorrelation-corrected
  statistics per window, drift, GB vs PB, hotspots, closure and entropy gates.
"""

GMX_MMPBSA_URL = "https://github.com/Valdes-Tresanca-MS/gmx_MMPBSA"
CITATIONS = [
    "Valdés-Tresanca, M. S.; Soler, M. A.; Moreno, E. et al. gmx_MMPBSA: A New Tool to "
    "Perform End-State Free Energy Calculations with GROMACS. J. Chem. Theory Comput. 2021, "
    "17, 6281–6291. https://doi.org/10.1021/acs.jctc.1c00645",
    "Miller, B. R.; McGee, T. D.; Swails, J. M. et al. MMPBSA.py: An Efficient Program for "
    "End-State Free Energy Calculations. J. Chem. Theory Comput. 2012, 8, 3314–3321. "
    "https://doi.org/10.1021/ct300418h",
]
