"""Grade a PSM spectrum against the Deutsch/TransCODE 2026 manual-validation rubric.

Classifies a spectrum as **excellent / good / caution / low_information** using objective
proxies for the paper's manual criteria. This grades the SPECTRUM, not the identification
(all reported PSMs already pass q<0.01 AND PEP<0.1) — it is the visual-corroboration check.

The "explained/benign" ion model (what a manual reviewer treats as legitimately present
without being a sequence ion) covers b/y/a at charge 1..min(2,z), H2O/NH3 losses, internal
by/ay ions, immonium ions, and the whole precursor family (charge-reduced + 13C isotopes +
multi-losses). This is what stops a precursor-loss or co-isolation base peak from being
mis-read as "unexplained". Ported from PSM_Verification grade_psms_rubric.py; automated
approximation of a manual procedure that also uses USIs and human judgement.
"""
import numpy as np
from .constants import AA, PROTON, H2O
from .annotate import clean, detect_tol, by_ions, annotate

CO = 27.994915
NH3 = 17.026549
ISO = 1.00336

GRADES = ("excellent", "good", "caution", "low_information")


def explained(pep, mods, z, pmz):
    """Return (mzs, labels) for every ion m/z a manual reviewer treats as benign/explained."""
    res = [AA[a] + mods.get(i, 0.0) for i, a in enumerate(pep)]
    nt = mods.get("N-term", 0.0)
    ct = mods.get("C-term", 0.0)
    L = len(pep)
    out = []
    mc = min(2, max(1, z))
    for i in range(1, L):
        bs = sum(res[:i]) + nt
        ys = sum(res[i:]) + ct + H2O
        for c in range(1, mc + 1):
            out += [((bs + c * PROTON) / c, f"b{i}^{c}"), ((ys + c * PROTON) / c, f"y{L-i}^{c}"),
                    ((bs - CO + c * PROTON) / c, f"a{i}^{c}")]
            out += [((bs - H2O + c * PROTON) / c, f"b{i}-H2O"), ((bs - NH3 + c * PROTON) / c, f"b{i}-NH3"),
                    ((ys - H2O + c * PROTON) / c, f"y{L-i}-H2O"), ((ys - NH3 + c * PROTON) / c, f"y{L-i}-NH3")]
    for i in range(1, L - 1):                                   # internal fragments
        for j in range(i + 1, L):
            m = sum(res[i:j])
            for c in (1, 2):
                out += [((m + c * PROTON) / c, f"int[{i}:{j}]by^{c}"), ((m - CO + c * PROTON) / c, f"int[{i}:{j}]ay^{c}")]
    for a in set(pep):                                          # immonium
        out.append((AA[a] - 26.00307, f"imm{a}"))
    M = pmz * z - z * PROTON                                    # precursor family
    LOSSES = ((0, ""), (-H2O, "-H2O"), (-NH3, "-NH3"), (-2 * H2O, "-2H2O"),
              (-H2O - NH3, "-H2O-NH3"), (-2 * NH3, "-2NH3"), (-3 * H2O, "-3H2O"))
    for c in range(1, z + 1):
        for iso, il in ((0, ""), (ISO, "+1"), (2 * ISO, "+2"), (-ISO, "-1")):
            for loss, ll in LOSSES:
                out.append(((M + iso + loss + c * PROTON) / c, f"prec{c}{il}{ll}"))
    out = [(m, l) for m, l in out if m > 0]
    out.sort()
    return np.array([m for m, l in out]), [l for m, l in out]


def unexplained(cs, emz, elab, tol, floor_frac=0.02):
    """Return (base_explained, max_unexplained_rel, n_unexplained, base_assignment, base_mz)."""
    mz, it = cs["mz"], cs["it"]
    if not len(it):
        return False, 0.0, 0, "none", "(empty)"
    base = it.max()
    floor = floor_frac * base

    def assign(m):
        lo = np.searchsorted(emz, m - tol)
        hi = np.searchsorted(emz, m + tol)
        if hi > lo:
            return elab[lo + int(np.argmin(np.abs(emz[lo:hi] - m)))]
        return None

    unex = [it[j] / base for j in range(len(mz)) if it[j] >= floor and assign(mz[j]) is None]
    bidx = int(np.argmax(it))
    ba = assign(mz[bidx])
    return (ba is not None), (max(unex) if unex else 0.0), len(unex), (ba or "UNEXPLAINED"), f"{mz[bidx]:.3f}"


def grade(spec, pep, mods, z, nox, ppm, tmt=False, tol=None):
    """Grade one PSM spectrum. Returns a dict with ``grade`` plus the component metrics.

    ``mods``/``nox``/``ppm`` are the resolved peptidoform values (from ``infer_oxidations``)
    — pass the same ones used for plotting so the grade matches the figure. ``tol`` defaults
    to the empirical per-spectrum tolerance (``detect_tol``).
    """
    cs = clean(spec)
    mc = min(2, max(1, z))
    if tol is None:
        tol = detect_tol(cs["mz"], by_ions(pep, mods, mc))
    ann = annotate(cs, pep, mods, maxch=mc, tol=tol)
    L = len(pep)
    sites = set(range(1, L))
    covered = set()
    bnums = set()
    for j, (t, num, c, imz), inten, d in ann["matched"]:
        if t == "b":
            covered.add(num)
            bnums.add(num)
        else:
            covered.add(L - num)
    covered &= sites
    cov = len(covered) / max(len(sites), 1)
    gap = g = 0
    for s in range(1, L):
        if s in covered:
            g = 0
        else:
            g += 1
            gap = max(gap, g)
    nb = len(bnums)
    ny = len({num for j, (t, num, c, i), ii, d in ann["matched"] if t == "y"})
    b2 = (2 in bnums)
    emz, elab = explained(pep, mods, spec.get("pz") or z, spec["pmz"])
    base_exp, max_unexp, n_unexp, base_assign, base_mz = unexplained(cs, emz, elab, tol[1])
    # per-criterion pass/fail (the 7 conjuncts of "excellent"), plus the two down-grade triggers
    cov_ok = cov >= 0.90                                    # near-complete sequence coverage
    gap_ok = gap <= 1                                       # no big break in the ladder
    series_ok = b2 or ny >= L - 1                           # both termini / a (near-)complete y-series
    ppm_ok = abs(ppm) <= 12 or (tmt and abs(ppm) <= 13)     # precursor mass accuracy (TMT ~-12 ppm offset)
    nox_ok = nox <= 2                                       # not over-modified
    base_ok = base_exp                                      # the base peak is an explained/benign ion
    unexp_ok = max_unexp < 0.5                              # no dominant unexplained peak
    caution_trigger = (max_unexp >= 0.6 and cov < 0.85)     # dominant stray peak AND incomplete ladder
    lowinfo_trigger = (cov < 0.55)                          # too few sequence ions
    if caution_trigger:
        grade_ = "caution"
    elif lowinfo_trigger:
        grade_ = "low_information"
    elif cov_ok and gap_ok and series_ok and ppm_ok and nox_ok and base_ok and unexp_ok:
        grade_ = "excellent"
    else:
        grade_ = "good"
    return dict(grade=grade_, cov=round(cov, 2), gap=gap, nb=nb, ny=ny, L=L, b2=b2,
                base_exp=base_exp, max_unexp=round(max_unexp, 2), n_unexp=n_unexp, ppm=round(ppm, 1),
                nox=nox, tol=tol[1], base_assign=base_assign, base_mz=base_mz,
                cov_ok=cov_ok, gap_ok=gap_ok, series_ok=series_ok, ppm_ok=ppm_ok, nox_ok=nox_ok,
                base_ok=base_ok, unexp_ok=unexp_ok,
                caution_trigger=caution_trigger, lowinfo_trigger=lowinfo_trigger)
