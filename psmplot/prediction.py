"""MS2PIP fragment-intensity prediction and observed-vs-predicted correlation.

Key gotcha preserved from the original: MS2PIP's returned fragment m/z ignore our
modification masses, so we remap its per-ion *intensities* onto our own correctly
computed b/y m/z (via the ion labels). This is what fixed negative TMT correlations.
"""
import re
import numpy as np
from .constants import CAM, OXM, TMT, TMTPRO
from .annotate import by_ions

_SP = None


def _sp():
    """Lazily construct a single MS2PIP SinglePrediction with our mod definitions."""
    global _SP
    if _SP is None:
        from ms2pip.single_prediction import SinglePrediction
        _SP = SinglePrediction()
        _SP.mod_info.add_from_ms2pip_modstrings(
            [f"Carbamidomethyl,{CAM},opt,C", f"Oxidation,{OXM},opt,M",
             f"TMT6plex,{TMT},opt,K", f"TMT6plexNterm,{TMT},opt,N-term",
             f"TMTpro,{TMTPRO},opt,K", f"TMTproNterm,{TMTPRO},opt,N-term"])
    return _SP


def ms2pip_modstr(pep, mods):
    """Translate our mods dict into an MS2PIP modification string ('0|Name|i|Name|...')."""
    parts = []
    nt = mods.get("N-term")
    if nt and abs(nt - TMTPRO) < 0.5:
        parts.append("0|TMTproNterm")
    elif nt and abs(nt - TMT) < 0.5:
        parts.append("0|TMT6plexNterm")
    for i, a in enumerate(pep):
        d = mods.get(i, 0.0)
        if a == "C" and abs(d - CAM) < 0.5:
            parts.append(f"{i+1}|Carbamidomethyl")
        elif a == "M" and abs((d % TMT) - OXM) < 0.5 or (a == "M" and abs(d - OXM) < 0.5):
            parts.append(f"{i+1}|Oxidation")
        elif a == "K" and abs(d - TMTPRO) < 0.5:
            parts.append(f"{i+1}|TMTpro")
        elif a == "K" and abs(d - TMT) < 1:
            parts.append(f"{i+1}|TMT6plex")
    return "|".join(parts) if parts else "-"


def predict_ms2pip(pep, mods, z, model):
    """Return MS2PIP-predicted ``(mz, intensity, labels)`` remapped onto our b/y m/z.

    Returns ``(None, None, None)`` if MS2PIP fails or yields nothing usable.
    """
    try:
        _, inten, ann = _sp().predict(pep, ms2pip_modstr(pep, mods), z, model=model)
    except Exception:
        return None, None, None
    frag = {(t, n, c): mz for t, n, c, mz in by_ions(pep, mods, 1)}   # singly-charged, with mods
    mz, it, labs = [], [], []
    for a, i in zip(ann, inten):
        m = re.match(r"([by])(\d+)", str(a))
        if not m:
            continue
        key = (m.group(1), int(m.group(2)), 1)
        if key in frag:
            mz.append(frag[key])
            it.append(float(i))
            labs.append(str(a))
    if not mz:
        return None, None, None
    return np.asarray(mz), np.asarray(it), labs


def mirror_corr(omz, oit, pmz, pit, tol=0.05):
    """Pearson r between observed and MS2PIP-predicted intensities over predicted ions."""
    ov, pv = [], []
    for m_, i_ in zip(pmz, pit):
        d = np.abs(omz - m_)
        j = int(np.argmin(d))
        ov.append(oit[j] if d[j] <= tol else 0.0)
        pv.append(i_)
    if len(pv) < 3 or np.std(ov) == 0:
        return float("nan")
    return float(np.corrcoef(ov, pv)[0, 1])
