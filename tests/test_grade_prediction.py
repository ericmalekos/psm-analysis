"""Tests of the rubric's immonium ions. Run with ``python tests/test_grade_prediction.py`` (or pytest)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import psmplot  # noqa: E402
from psmplot.grade import explained  # noqa: E402

# immonium m/z (residue - CO + H+)
IMM = {"L": 86.09643, "H": 110.07127, "F": 120.08078, "Y": 136.07569, "K": 101.10732, "P": 70.06513}


def test_immonium():
    pep = "LHFYPK"
    emz, elab = explained(pep, {}, 2, 400.0)
    imm = {l[3:]: m for m, l in zip(emz, elab) if l.startswith("imm")}
    for a, m in IMM.items():
        assert abs(imm[a] - m) < 1e-3, (a, imm[a], m)
    # residues carry their modifications: oxidised Met 120.0478, carbamidomethyl Cys 133.0431
    emz, elab = explained("MC", {0: 15.994915, 1: 57.021464}, 1, 300.0)
    imm = sorted(m for m, l in zip(emz, elab) if l.startswith("imm"))
    assert abs(imm[0] - 120.0478) < 1e-3 and abs(imm[1] - 133.0431) < 1e-3, imm


def test_immonium_explains_base_peak():
    # a His immonium base peak at the C-terminal residue is explained (it is no internal fragment)
    pep, mods = "AGLSVAKH", {}
    ions = [m for t, n, c, m in psmplot.by_ions(pep, mods, 1)]
    mz = np.array(sorted(ions + [IMM["H"]]))
    it = np.where(np.abs(mz - IMM["H"]) < 1e-6, 10.0, 1.0)
    pmz = (psmplot.neutral_mass(pep, mods) + 2 * psmplot.PROTON) / 2
    g = psmplot.grade(dict(mz=mz, it=it, pmz=pmz, pz=2), pep, mods, 2, 0, 0.0, tol=("Da", 0.025))
    assert g["base_exp"] and g["base_assign"] == "immH", g["base_assign"]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
