"""Tests of the TMT/TMTpro handling. Run with ``python tests/test_tmt.py`` (or pytest)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import psmplot  # noqa: E402
from psmplot.constants import CO, PROTON, TMT, TMTPRO, TMT_REPORTERS, TMTPRO_REPORTERS  # noqa: E402


def test_tmt_kind():
    assert psmplot.tmt_kind(False) is None and psmplot.tmt_kind(None) is None and psmplot.tmt_kind("") is None
    assert psmplot.tmt_kind(True) == "tmt" and psmplot.tmt_kind("TMT") == "tmt"
    assert psmplot.tmt_kind("tmtpro") == "tmtpro"
    assert psmplot.tmt_kind(True, {"N-term": TMTPRO}) == "tmtpro"     # True resolves from the label mass
    try:
        psmplot.tmt_kind("itraq")
    except ValueError:
        pass
    else:
        raise AssertionError("an unknown label must raise")


def test_base_mods():
    assert psmplot.base_mods("PEPCK") == {3: 57.021464}                # unlabelled: as before
    m = psmplot.base_mods("PEPCK", tmt=True)                           # TMT 6/10/11-plex: as before
    assert m["N-term"] == TMT and m[4] == TMT and m[3] == 57.021464
    m = psmplot.base_mods("PEPCK", tmt="tmtpro")
    assert m["N-term"] == TMTPRO and m[4] == TMTPRO


def test_tmt_ions():
    assert psmplot.tmt_ions(None) == []
    ions = psmplot.tmt_ions("tmt")
    rep = [m for m, l in ions if l == "TMTrep"]
    assert rep[:11] == list(TMT_REPORTERS) and len(rep) == 22                    # and their 13C isotopes
    assert abs(rep[11] - 127.131081) < 1e-6 and abs(rep[-1] - 132.147854) < 1e-5
    tag = [m for m, l in ions if l == "TMTtag"]
    assert abs(tag[0] - 230.1702) < 1e-3 and abs(tag[1] - 248.1808) < 1e-3       # tag ion, tag + H2O
    assert abs(tag[2] - 231.1736) < 1e-3                                          # 13C isotope of the tag
    assert len([1 for _, l in psmplot.tmt_ions("tmtpro") if l == "TMTrep"]) == 36
    assert min(abs(m - 305.2144) for m, l in psmplot.tmt_ions("tmtpro") if l == "TMTtag") < 1e-3
    assert min(abs(m - 323.2250) for m, l in psmplot.tmt_ions("tmtpro") if l == "TMTtag") < 1e-3
    # complement ions of a TMTpro SLGAGPK 2+ (N-term and K labelled): precursor less reporter and CO, 1+
    mods = psmplot.base_mods("SLGAGPK", tmt="tmtpro")
    pmz = (psmplot.neutral_mass("SLGAGPK", mods) + 2 * PROTON) / 2
    comp = [m for m, l in psmplot.tmt_ions("tmtpro", pmz, 2) if l == "TMTc"]
    assert len(comp) == 2 * 18                                         # each reporter, with its 13C isotope
    assert min(abs(c - 1082.660) for c in comp) < 0.01                 # 128N channel; observed 1082.659
    assert not [1 for _, l in psmplot.tmt_ions("tmtpro", pmz, 1) if l == "TMTc"]   # none for 1+ precursors


def test_drop_tmt():
    pep = "SLGAGPK"
    mods = psmplot.base_mods(pep, tmt="tmtpro")
    pmz = (psmplot.neutral_mass(pep, mods) + 2 * PROTON) / 2
    seq = [m for t, n, c, m in psmplot.by_ions(pep, mods, 1)] + [132.10191]   # b/y ions and a y1 of Leu
    tmt = list(TMTPRO_REPORTERS) + [TMTPRO + PROTON, TMTPRO + PROTON + 1.003355, TMTPRO + PROTON + 18.010565,
                                    pmz * 2 - 128.128116 - CO]
    spec = dict(mz=np.array(sorted(seq + tmt)), it=np.ones(len(seq) + len(tmt)), pmz=pmz, pz=2)
    out = psmplot.drop_tmt(spec, "tmtpro", 2, 0.025)
    assert sorted(out["mz"].tolist()) == sorted(seq)                   # every TMT ion gone, every other peak kept
    assert psmplot.drop_tmt(spec, None) is spec                        # unlabelled: untouched
    # a complement ion that coincides with a sequence ion is kept (TDKLSFK 3+: y4 798.504, 128N complement 798.489)
    pep = "TDKLSFK"
    mods = psmplot.base_mods(pep, tmt="tmtpro")
    pmz = (psmplot.neutral_mass(pep, mods) + 3 * PROTON) / 3
    by = [m for t, n, c, m in psmplot.by_ions(pep, mods, 2)]
    y4 = [m for t, n, c, m in psmplot.by_ions(pep, mods, 1) if (t, n) == ("y", 4)][0]
    comp = (pmz * 3 - 128.128116 - CO) / 2
    assert abs(y4 - comp) < 0.025
    s2 = dict(mz=np.array([y4]), it=np.ones(1), pmz=pmz, pz=3)
    assert psmplot.drop_tmt(s2, "tmtpro", 3, 0.025)["mz"].tolist() == []                  # unprotected: removed
    assert psmplot.drop_tmt(s2, "tmtpro", 3, 0.025, keep_mz=by)["mz"].tolist() == [y4]    # protected: kept
    assert psmplot.drop_tmt(s2, "tmtpro", 3, 0.4)["mz"].tolist() == [y4]                  # low resolution: no complements


def _ladder_spectrum(pep, mods, z, keep_b, keep_y, extra):
    """synthetic spectrum: the chosen b/y ions at intensity 1, plus ``extra`` [(mz, intensity)]"""
    ions = [(m, 1.0) for t, n, c, m in psmplot.by_ions(pep, mods, 1)
            if (t == "b" and n in keep_b) or (t == "y" and n in keep_y)]
    pk = sorted(ions + extra)
    pmz = (psmplot.neutral_mass(pep, mods) + z * PROTON) / z
    return dict(mz=np.array([m for m, _ in pk]), it=np.array([i for _, i in pk]), pmz=pmz, pz=z)


def test_grade_tmt():
    pep = "AEATESAMEK"
    mods = psmplot.base_mods(pep, tmt="tmtpro")
    spec = _ladder_spectrum(pep, mods, 2, keep_b={2, 3, 4}, keep_y={1, 2, 3},        # 6 of 9 sites: cov 0.67
                            extra=[(r, 20.0) for r in TMTPRO_REPORTERS] + [(TMTPRO + PROTON, 8.0)])
    plain = psmplot.grade(spec, pep, mods, 2, 0, 0.0, tmt=False, tol=("Da", 0.025))
    labelled = psmplot.grade(spec, pep, mods, 2, 0, 0.0, tmt="tmtpro", tol=("Da", 0.025))
    assert plain["grade"] == "caution" and plain["base_assign"] == "UNEXPLAINED"   # reporter base peak
    assert labelled["tmt"] == "tmtpro" and labelled["base_exp"] and labelled["max_unexp"] == 0.0
    assert plain["ppm_ok"] is True and labelled["ppm_ok"] is True
    far = psmplot.grade(spec, pep, mods, 2, 0, 12.5, tol=("Da", 0.025))            # 12.5 ppm: only TMT allows it
    assert far["ppm_ok"] is False and psmplot.grade(spec, pep, mods, 2, 0, 12.5, tmt="tmtpro", tol=("Da", 0.025))["ppm_ok"] is True
    assert psmplot.grade(spec, pep, mods, 2, 0, 12.5, tmt=None, tol=("Da", 0.025))["ppm_ok"] is False   # a bool, not None
    assert labelled["grade"] == "good" and labelled["cov"] == plain["cov"]          # same ladder, no TMT ions
    # an unlabelled spectrum grades the same with and without the TMT branch untouched
    m0 = psmplot.base_mods(pep)
    s0 = _ladder_spectrum(pep, m0, 2, keep_b={2, 3, 4, 5, 6, 7, 8, 9}, keep_y=set(range(1, 10)), extra=[])
    assert psmplot.grade(s0, pep, m0, 2, 0, 1.0)["grade"] == "excellent"


def test_ms2pip_modstr():
    assert psmplot.ms2pip_modstr("PEPCK", psmplot.base_mods("PEPCK", tmt=True)) == "0|TMT6plexNterm|4|Carbamidomethyl|5|TMT6plex"
    assert psmplot.ms2pip_modstr("PEPCK", psmplot.base_mods("PEPCK", tmt="tmtpro")) == "0|TMTproNterm|4|Carbamidomethyl|5|TMTpro"
    assert psmplot.ms2pip_modstr("PEPCK", psmplot.base_mods("PEPCK")) == "4|Carbamidomethyl"


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
