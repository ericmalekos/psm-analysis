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
C13    = 1.003355           # exact 13C-12C spacing (first isotope of the TMT reporter, tag and complement ions)

CO     = 27.994915          # carbon monoxide (Da), lost with the reporter ion from a TMT complement ion

# modification mass deltas (Da)
CAM = 57.021464            # Carbamidomethyl on C   (fixed, iodoacetamide alkylation)
OXM = 15.994915            # Oxidation on M         (variable)
TMT = 229.162932           # TMT 6/10/11-plex on N-term and K (only PXD029155 in this project)
TMTPRO = 304.207146        # TMTpro 16/18-plex on N-term and K

# TMT reporter ions (m/z, 1+). TMT 6/10/11-plex use the first 11 (TMT6: 126, 127N, 128C, 129N, 130C, 131N),
# TMTpro 16/18-plex all 18. They carry no sequence information and are removed before annotation and grading
# (annotate.drop_tmt), with the intact tag ion and the complement ions.
TMT_REPORTERS = (126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471,
                 129.137790, 130.134825, 130.141145, 131.138180, 131.144499)
TMTPRO_REPORTERS = TMT_REPORTERS + (132.141535, 132.147855, 133.144890, 133.151210,
                                    134.148245, 134.154565, 135.151600)
