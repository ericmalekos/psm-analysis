"""Spectrum loading for MGF and mzML, with per-file caching and scan indexing.

``load_spectrum(path, scan)`` returns a dict::

    {"mz": ndarray, "it": ndarray, "pmz": float|None, "pz": int|None, "title": str}

MGF files are read via ``IndexedMGF`` and a scan->title map is built once per file.
mzML files (optionally gzipped) are decompressed once into ``CACHE_DIR``.
"""
import os
import re
import gzip
import json
import tempfile
import numpy as np
from pyteomics import mgf, mzml

# Where gzipped mzML files are decompressed. Override with $PSMPLOT_CACHE.
CACHE_DIR = os.environ.get("PSMPLOT_CACHE", tempfile.gettempdir())

_MGF = {}


def _load_mgf(path, scan):
    if path not in _MGF:
        f = mgf.IndexedMGF(path)
        s2t = {}
        for k in f.index:
            m = re.search(r"\.(\d+)\.\d+\.\d+", k) or re.search(r"scan[=:](\d+)", k)
            if m:
                s2t.setdefault(m.group(1), k)
        _MGF[path] = (f, s2t)
    f, s2t = _MGF[path]
    t = s2t.get(str(scan))
    if t is None:
        return None
    s = f[t]
    ch = s["params"].get("charge")
    return dict(mz=np.asarray(s["m/z array"], float), it=np.asarray(s["intensity array"], float),
                pmz=float(s["params"]["pepmass"][0]), pz=int(ch[0]) if ch else None,
                title=s["params"].get("title", ""))


_MZML = {}


def _load_mzml(path, scan):
    if path not in _MZML:
        local = path
        if path.endswith(".gz"):
            os.makedirs(CACHE_DIR, exist_ok=True)
            local = os.path.join(CACHE_DIR, os.path.basename(path)[:-3])
            if not os.path.exists(local):
                with gzip.open(path, "rb") as fi, open(local, "wb") as fo:
                    fo.write(fi.read())
        _MZML[path] = mzml.MzML(local)
    f = _MZML[path]
    s = None
    for sid in (f"controllerType=0 controllerNumber=1 scan={scan}", f"scan={scan}", f"index={int(scan)-1}"):
        try:
            s = f.get_by_id(sid)
            break
        except Exception:
            s = None
    if s is None:
        for s in f:
            if s.get("id", "").endswith(f"scan={scan}"):
                break
        else:
            return None
    pmz = pz = None
    try:
        pre = s["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]
        pmz = float(pre.get("selected ion m/z"))
        pz = int(pre.get("charge state")) if "charge state" in pre else None
    except Exception:
        pass
    return dict(mz=np.asarray(s["m/z array"], float), it=np.asarray(s["intensity array"], float),
                pmz=pmz, pz=pz, title=s.get("id", ""))


def _load_proxi(path, scan):
    """Load a scan from a PROXI/USI JSON (e.g. ThermoRawFileParser `query` output).

    The file is a list of spectrum objects, each with ``mzs``/``intensities`` and an
    ``attributes`` list (scan number, charge state, selected ion m/z). Selects the entry
    whose scan number matches ``scan`` (or the sole/first entry if ``scan`` is None).
    """
    arr = json.load(open(path))
    if isinstance(arr, dict):
        arr = [arr]

    def scannum(s):
        for x in s.get("attributes", []):
            if x.get("name") == "scan number":
                return str(x.get("value"))
        return None

    s = None
    if scan is not None:
        for x in arr:
            if scannum(x) == str(scan):
                s = x
                break
    if s is None:
        s = arr[0] if arr else None
    if s is None:
        return None
    a = {x["name"]: x.get("value") for x in s.get("attributes", []) if "name" in x}
    pz = a.get("charge state")
    pmz = a.get("selected ion m/z")
    filt = str(a.get("filter string") or "")
    analyzer = "FTMS" if "FTMS" in filt else "ITMS" if "ITMS" in filt else "?"
    return dict(mz=np.asarray(s["mzs"], float), it=np.asarray(s["intensities"], float),
                pmz=float(pmz) if pmz is not None else None, pz=int(pz) if pz else None,
                analyzer=analyzer, filter=filt, title=f"scan {a.get('scan number')}")


def load_spectrum(path, scan):
    """Load one scan from an MGF, mzML(.gz), or PROXI/USI JSON file. None if the scan is absent."""
    p = path.lower()
    if p.endswith(".json"):
        return _load_proxi(path, scan)
    if p.endswith((".mzml", ".mzml.gz")):
        return _load_mzml(path, scan)
    return _load_mgf(path, scan)
