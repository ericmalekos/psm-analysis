"""Physical constants and modification mass deltas used throughout psmplot.

These are chemistry, not design — they should essentially never change. Visual
choices live in ``theme.py``; project-specific choices (which MS2PIP model goes
with which dataset) live in ``datasets.py``.
"""
from pyteomics import mass

# monoisotopic residue masses (Da), from pyteomics' standard table
AA = dict(mass.std_aa_mass)

PROTON = 1.00727646688      # mass of a proton (Da)
H2O    = 18.0105646863      # monoisotopic water (Da)
ISO    = 1.00286            # ~C13-C12 spacing, used to tolerate off-isotope precursor picks

# modification mass deltas (Da)
CAM = 57.021464            # Carbamidomethyl on C   (fixed, iodoacetamide alkylation)
OXM = 15.994915            # Oxidation on M         (variable)
TMT = 229.162932           # TMT 6-plex on N-term and K (only PXD029155 in this project)
