"""MS2PIP fragment-intensity prediction and observed-vs-predicted correlation.

Key gotcha preserved from the original: MS2PIP's returned fragment m/z ignore our
modification masses, so we remap its per-ion *intensities* onto our own correctly
computed b/y m/z (via the ion labels). This is what fixed negative TMT correlations.

MS2PIP 3 (``ms2pip.single_prediction``) is used when it is installed, else MS2PIP 4, which
takes the peptidoform as ProForma with every modification as a mass shift. MS2PIP 4 loads
its model files on every single prediction (seconds for the XGBoost models), so batch
callers first call :func:`predict_ms2pip_batch`, which predicts all PSMs in one call per
model and caches them for :func:`predict_ms2pip`. ``$PSMPLOT_MS2PIP_MODEL_DIR`` overrides
MS2PIP 4's model directory (default ``~/.ms2pip``).
"""
import os
import re
import numpy as np
from .constants import CAM, OXM, TMT, TMTPRO
from .annotate import by_ions

_SP = None
_V3 = None
_CACHE = {}      # MS2PIP 4: (proforma, model) -> {"b": intensities, "y": intensities}, linear


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


def _has_v3():
    """True when MS2PIP 3 (``ms2pip.single_prediction``) is installed; MS2PIP 4 is used otherwise."""
    global _V3
    if _V3 is None:
        try:
            import ms2pip.single_prediction  # noqa: F401
            _V3 = True
        except ImportError:
            _V3 = False
    return _V3


def proforma(pep, mods, z):
    """ProForma string with every modification as a mass shift, e.g. ``[+229.162932]-PEPC[+57.021464]K/2``."""
    s = "".join(a + (f"[{mods[i]:+.6f}]" if mods.get(i) else "") for i, a in enumerate(pep))
    if mods.get("N-term"):
        s = f"[{mods['N-term']:+.6f}]-" + s
    if mods.get("C-term"):
        s += f"-[{mods['C-term']:+.6f}]"
    return f"{s}/{int(z)}"


def _linear(pred):
    """MS2PIP 4 intensities per ion type, log2(I + 0.001), back to linear."""
    return {k.lower(): np.clip(2 ** np.asarray(v, float) - 0.001, 0, None) for k, v in pred.items()}


def _model_dir():
    return os.environ.get("PSMPLOT_MS2PIP_MODEL_DIR") or None


def predict_ms2pip_batch(items, processes=1):
    """Predict many PSMs with MS2PIP 4, one call per model, and cache them for :func:`predict_ms2pip`.

    ``items``: iterable of ``(pep, mods, z, model)``. Returns the number of peptidoforms predicted; a
    model whose batch fails is left to single predictions. No-op under MS2PIP 3. Callers that use
    ``processes`` > 1 need the usual ``if __name__ == "__main__":`` guard.
    """
    if _has_v3():
        return 0
    from ms2pip import predict_batch
    from psm_utils import PSM, PSMList, Peptidoform
    todo = {}
    for pep, mods, z, model in items:
        p = proforma(pep, mods, z)
        if (p, model) not in _CACHE:
            todo.setdefault(model, set()).add(p)
    n = 0
    for model, pfs in todo.items():
        pfs = sorted(pfs)
        try:
            res = predict_batch(PSMList(psm_list=[PSM(peptidoform=Peptidoform(p), spectrum_id=i)
                                                  for i, p in enumerate(pfs)]),
                                model=model, model_dir=_model_dir(), processes=processes)
        except Exception:
            continue
        for p, r in zip(pfs, res):
            if r.predicted_intensity is not None:
                _CACHE[(p, model)] = _linear(r.predicted_intensity)
                n += 1
    return n


def _predict_v4(pep, mods, z, model):
    key = (proforma(pep, mods, z), model)
    if key not in _CACHE:
        from ms2pip import predict_single
        _CACHE[key] = _linear(predict_single(key[0], model=model, model_dir=_model_dir()).predicted_intensity)
    return _CACHE[key]


def predict_ms2pip(pep, mods, z, model):
    """Return MS2PIP-predicted ``(mz, intensity, labels)`` remapped onto our b/y m/z.

    MS2PIP 3 when installed, else MS2PIP 4 (from the :func:`predict_ms2pip_batch` cache, or a single
    prediction). Returns ``(None, None, None)`` if MS2PIP fails or yields nothing usable.
    """
    try:
        if _has_v3():
            _, inten, ann = _sp().predict(pep, ms2pip_modstr(pep, mods), z, model=model)
            pairs = list(zip(ann, inten))
        else:
            pred = _predict_v4(pep, mods, z, model)
            pairs = [(f"{t}{j + 1}", v) for t in ("b", "y") for j, v in enumerate(pred.get(t, []))]
    except Exception:
        return None, None, None
    frag = {(t, n, c): mz for t, n, c, mz in by_ions(pep, mods, 1)}   # singly-charged, with mods
    mz, it, labs = [], [], []
    for a, i in pairs:
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
