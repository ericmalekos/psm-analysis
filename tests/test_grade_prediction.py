"""Tests of the rubric's immonium ions and of the MS2PIP prediction. Run with ``python tests/test_grade_prediction.py``
(or pytest). The prediction tests need MS2PIP 4 and are skipped without it; they use HCD2019, which needs no model files."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import psmplot  # noqa: E402
import psmplot.prediction as pred  # noqa: E402
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


def test_proforma():
    assert psmplot.proforma("PEPCK", psmplot.base_mods("PEPCK", tmt=True), 2) == \
        "[+229.162932]-PEPC[+57.021464]K[+229.162932]/2"
    assert psmplot.proforma("PEPMK", {3: 15.994915}, 3) == "PEPM[+15.994915]K/3"
    assert psmplot.proforma("PEPK", {}, 1) == "PEPK/1"


def _ms2pip4():
    try:
        import ms2pip  # noqa: F401
        return not pred._has_v3()
    except ImportError:
        return False


def test_predict_ms2pip4():
    if not _ms2pip4():
        print("skip: MS2PIP 4 not installed")
        return
    pep, mods = "PEPTIDEK", {}
    mz, it, lab = psmplot.predict_ms2pip(pep, mods, 2, "HCD2019")
    assert mz is not None and len(mz) == 14 and sorted(lab) == sorted([f"b{i}" for i in range(1, 8)] + [f"y{i}" for i in range(1, 8)])
    frag = {f"{t}{n}": m for t, n, c, m in psmplot.by_ions(pep, mods, 1)}
    assert all(abs(m - frag[l]) < 1e-9 for m, l in zip(mz, lab))      # our own b/y m/z
    assert (it >= 0).all() and it.max() > 0
    # modified peptide: masses carry the mods, the prediction runs
    m2 = psmplot.base_mods("PEPCK", tmt="tmtpro")
    mz2, it2, lab2 = psmplot.predict_ms2pip("PEPCK", m2, 2, "HCD2019")
    assert mz2 is not None and abs(dict(zip(lab2, mz2))["y1"] - (147.11280 + 304.207146)) < 1e-3


def test_predict_ms2pip4_batch():
    if not _ms2pip4():
        print("skip: MS2PIP 4 not installed")
        return
    items = [("ACDEFGHIK", {1: 57.021464}, 2, "HCD2019"), ("LLMNPQR", {2: 15.994915}, 3, "HCD2019")]
    assert psmplot.predict_ms2pip_batch(items) == 2
    import ms2pip
    real = ms2pip.predict_single
    ms2pip.predict_single = lambda *a, **k: (_ for _ in ()).throw(AssertionError("cache not used"))
    try:
        for pep, mods, z, model in items:                               # served from the cache
            assert psmplot.predict_ms2pip(pep, mods, z, model)[0] is not None
    finally:
        ms2pip.predict_single = real
    single = psmplot.predict_ms2pip("ACDEFGHIK", {1: 57.021464}, 2, "HCD2019")
    pred._CACHE.clear()
    again = psmplot.predict_ms2pip("ACDEFGHIK", {1: 57.021464}, 2, "HCD2019")
    assert np.allclose(single[1], again[1])                             # batch and single agree


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
