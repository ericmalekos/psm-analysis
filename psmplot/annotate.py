"""Fragment chemistry: modifications, b/y ions, adaptive tolerance, and peak annotation.

This is the "what matches what" layer. It has no matplotlib dependency, so the
annotation can be reused for grading/QC independently of plotting.
"""
import numpy as np
from .constants import (AA, PROTON, H2O, ISO, C13, CO, CAM, OXM, TMT, TMTPRO, TMT_REPORTERS,
                        TMTPRO_REPORTERS)


# ---------------- modifications ----------------
def tmt_kind(tmt, mods=None):
    """Resolve the TMT label: ``None`` (unlabelled), ``"tmt"`` (6/10/11-plex) or ``"tmtpro"`` (16/18-plex).

    ``tmt`` may be False/None, ``"tmt"``, ``"tmtpro"``, or True: TMT 6/10/11-plex, or TMTpro when
    ``mods`` already carry the TMTpro mass.
    """
    if not tmt:
        return None
    if isinstance(tmt, str):
        kind = tmt.lower()
        if kind not in ("tmt", "tmtpro"):
            raise ValueError(f"tmt must be 'tmt' or 'tmtpro', not {tmt!r}")
        return kind
    if mods and any(abs(d - TMTPRO) < 0.01 for d in mods.values()):
        return "tmtpro"
    return "tmt"


def base_mods(pep, tmt=False):
    """Deterministic modification dict ``{resIdx -> massdelta, 'N-term' -> ...}``.

    Fixed Carbamidomethyl on every C; the TMT label on the N-term and every K when ``tmt``
    (True or ``"tmt"``: TMT 6/10/11-plex; ``"tmtpro"``: TMTpro 16/18-plex).
    Variable Oxidation-M is added later by :func:`infer_oxidations`.
    """
    m = {i: CAM for i, a in enumerate(pep) if a == "C"}
    kind = tmt_kind(tmt)
    if kind:
        label = TMTPRO if kind == "tmtpro" else TMT
        m["N-term"] = label
        for i, a in enumerate(pep):
            if a == "K":
                m[i] = m.get(i, 0.0) + label
    return m


def neutral_mass(pep, mods):
    """Monoisotopic neutral mass M of a peptidoform."""
    return sum(AA[a] for a in pep) + H2O + sum(mods.values())


def infer_oxidations(pep, mods, obs_pmz, pz):
    """Add k Oxidation-M to best match the observed precursor mass.

    Returns ``(mods, n_ox, ppm)``. Tries 0..(#Met) oxidations and a few isotope
    offsets so an off-isotope precursor pick does not defeat the match.
    """
    if pz is None:
        return mods, 0, float("nan")
    obs_neutral = obs_pmz * pz - pz * PROTON
    n_m = pep.count("M")
    best = (mods, 0, 1e9)
    for k in range(0, n_m + 1):
        for iso in (0, 1, 2, -1):
            test = dict(mods)
            placed = 0
            for i, a in enumerate(pep):
                if a == "M" and placed < k:
                    test[i] = test.get(i, 0.0) + OXM
                    placed += 1
            nm = neutral_mass(pep, test) + iso * ISO
            ppm = (nm - obs_neutral) / obs_neutral * 1e6
            if abs(ppm) < abs(best[2]):
                best = (test, k, ppm)
    return best


# ---------------- fragments + annotation ----------------
def by_ions(pep, mods, maxch=2):
    """All b and y ions up to charge ``maxch``. Returns [(type, num, charge, m/z), ...]."""
    res = [AA[a] + mods.get(i, 0.0) for i, a in enumerate(pep)]
    nt = mods.get("N-term", 0.0)
    ct = mods.get("C-term", 0.0)
    n = len(pep)
    out = []
    for i in range(1, n):
        bs = sum(res[:i]) + nt
        ys = sum(res[i:]) + ct + H2O
        for c in range(1, maxch + 1):
            out.append(("b", i, c, (bs + c * PROTON) / c))
            out.append(("y", n - i, c, (ys + c * PROTON) / c))
    return out


def detect_tol(mz, ions):
    """Adaptive fragment tolerance: high-res iff >=5 singly-charged b/y match within 0.02 Da.

    Returns ``("Da", 0.025)`` for high-res Orbitrap MS2, else ``("Da", 0.4)`` for
    low-res ion-trap MS2. (Most datasets here are high-res-MS1 + low-res-CID-MS2.)
    """
    e = np.array([np.min(np.abs(mz - imz)) for t, n, c, imz in ions if c == 1])
    return ("Da", 0.025) if int((e < 0.02).sum()) >= 5 else ("Da", 0.4)


def clean(spec):
    """Drop the intact precursor (and immediate isotopes) so fragment signal dominates."""
    mz, it = spec["mz"], spec["it"]
    p = spec.get("pmz")
    keep = np.ones(len(mz), bool)
    if p:
        keep &= ~((mz > p - 1.6) & (mz < p + 3.2))
    return dict(spec, mz=mz[keep], it=it[keep])


def tmt_ions(kind, pmz=None, z=None):
    """m/z of the TMT ions that carry no sequence information, as ``[(mz, label), ...]``.

    The reporter ions (``TMTrep``), the intact tag ion (``TMTtag``: label + proton, and the tag +
    H2O) and, for a precursor of charge ``z`` >= 2 at ``pmz``, the complement ions (``TMTc``): the
    precursor less one reporter ion and CO, at charge z-1. Each with its first 13C isotope. Empty
    when ``kind`` is None.
    """
    if not kind:
        return []
    reporters = TMTPRO_REPORTERS if kind == "tmtpro" else TMT_REPORTERS
    tag = (TMTPRO if kind == "tmtpro" else TMT) + PROTON
    out = []
    for iso in (0.0, C13):
        out += [(r + iso, "TMTrep") for r in reporters]
        out += [(tag + iso, "TMTtag"), (tag + H2O + iso, "TMTtag")]
    if pmz and z and z >= 2:
        total = pmz * z
        for r in reporters:
            for iso in (0.0, C13):
                out.append(((total - r - CO + iso) / (z - 1), "TMTc"))
    return out


def drop_tmt(spec, kind, z=None, tol=0.025, keep_mz=None):
    """Remove TMT reporter, tag and complement peaks (within ``tol`` Da) so the sequence ions dominate,
    as :func:`clean` does for the precursor. A peak that also lies within ``tol`` of one of ``keep_mz``
    (the peptide's own b/y ions) is kept: complement ions in particular can coincide with a sequence ion.
    Complement ions are removed only at high resolution (``tol`` < 0.1 Da): their cluster (one ion per
    reporter, about 1/(z-1) apart) would cover several m/z units at a low-resolution tolerance and take
    neighbouring fragment ions with it. Returns ``spec`` unchanged when ``kind`` is None."""
    if not kind or not len(spec["mz"]):
        return spec
    mz = spec["mz"]
    ions = np.array([m for m, _ in tmt_ions(kind, spec.get("pmz") if tol < 0.1 else None, z)])
    drop = np.min(np.abs(mz[:, None] - ions[None, :]), axis=1) <= tol
    if keep_mz is not None and len(keep_mz):
        drop &= np.min(np.abs(mz[:, None] - np.asarray(keep_mz, float)[None, :]), axis=1) > tol
    return dict(spec, mz=mz[~drop], it=spec["it"][~drop])


def annotate(spec, pep, mods, maxch=2, tol=None, min_rel_int=0.02):
    """Match b/y ions to peaks above a noise floor (``min_rel_int`` of the base peak).

    Returns dict with ``tol``, ``matched`` [(peak_idx, ion, intensity, abs_err), ...],
    ``ions`` (all theoretical), ``covfrac`` (annotated fraction of TIC) and ``distinct``.
    """
    mz, it = spec["mz"], spec["it"]
    ions = by_ions(pep, mods, maxch)
    if tol is None:
        tol = detect_tol(mz, ions)
    mode, val = tol
    floor = min_rel_int * it.max() if len(it) else 0
    best = {}                                            # peak idx -> (ion, intensity, err)
    for t, num, c, imz in ions:
        d = np.abs(mz - imz)
        j = int(np.argmin(d))
        if it[j] < floor:
            continue
        if d[j] <= val and (j not in best or d[j] < best[j][2]):
            best[j] = ((t, num, c, imz), it[j], d[j])
    matched = [(j, ion, inten, d) for j, (ion, inten, d) in best.items()]
    covfrac = sum(m[2] for m in matched) / it.sum() if it.sum() else 0
    distinct = len({(ion[0], ion[1]) for _, ion, _, _ in matched})
    return dict(tol=tol, matched=matched, ions=ions, covfrac=covfrac, distinct=distinct)


def longest_series(pep, matched):
    """Length of the longest consecutive run of b ions and of y ions."""
    n = len(pep)
    byb = {ion[1] for _, ion, _, _ in matched if ion[0] == "b"}
    byy = {ion[1] for _, ion, _, _ in matched if ion[0] == "y"}

    def run(s):
        best = cur = 0
        for i in range(1, n):
            cur = cur + 1 if i in s else 0
            best = max(best, cur)
        return best

    return run(byb), run(byy)
