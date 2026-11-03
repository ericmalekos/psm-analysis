"""psmplot — annotate and mirror-plot peptide-spectrum matches (PSMs) for manual verification.

Typical use::

    import psmplot
    psmplot.render_psm(
        peptide="LMDLSHTNSVAELQETPGK", charge=3,
        spectrum_path=".../run.mgf", scan=30662,
        model="CID", out_png="fig.png",
        title="Gvin-ps7  LMDLSHTNSVAELQETPGK  3+")

Restyle by editing ``psmplot/theme.py`` or passing ``theme=Theme(...)``.
"""
from .constants import AA, PROTON, H2O, CAM, OXM, TMT, TMTPRO
from .spectra import load_spectrum
from .annotate import (base_mods, neutral_mass, infer_oxidations, by_ions,
                       detect_tol, clean, annotate, longest_series, tmt_kind, tmt_ions, drop_tmt)
from .prediction import predict_ms2pip, mirror_corr, ms2pip_modstr
from .plot import plot_psm
from .theme import Theme, DEFAULT, MINIMAL
from .datasets import model_for, is_tmt
from .grade import grade, GRADES

__all__ = [
    "render_psm", "load_spectrum", "plot_psm", "Theme", "DEFAULT", "MINIMAL",
    "base_mods", "infer_oxidations", "annotate", "clean", "by_ions",
    "detect_tol", "longest_series", "predict_ms2pip", "mirror_corr",
    "neutral_mass", "model_for", "is_tmt", "grade", "GRADES",
    "tmt_kind", "tmt_ions", "drop_tmt",
]


def render_psm(peptide, charge, spectrum_path, scan, out_png,
               model="HCD2019", tmt=False, title=None, subtitle="", theme=None, badge=False):
    """Resolve mods + charge for one PSM and render it. Returns a QC dict.

    Loads the scan, applies fixed Carbamidomethyl (and the TMT label if ``tmt``: True or
    ``"tmt"`` for TMT 6/10/11-plex, ``"tmtpro"`` for TMTpro 16/18-plex), infers variable
    Oxidation-M from the precursor mass, then draws the annotated + mirror figure. For TMT the
    reporter, tag and complement ions are removed before annotation and grading.
    Raises ``ValueError`` if the scan cannot be loaded.

    The returned dict carries: ``nb, ny, bmax, ymax`` (b/y counts + longest runs),
    ``covfrac`` (annotated fraction of fragment TIC), ``n_ox``, ``ppm`` (precursor
    mass error), ``highres`` (bool), ``tol`` and ``ms2pip_r``.
    """
    charge = int(charge)
    spec = load_spectrum(spectrum_path, scan)
    if spec is None or not len(spec["mz"]):
        raise ValueError(f"scan {scan} not found (or empty) in {spectrum_path}")
    mods = base_mods(peptide, tmt=tmt)
    mods, n_ox, ppm = infer_oxidations(peptide, mods, spec["pmz"], spec["pz"] or charge)
    grade_info = grade(spec, peptide, mods, charge, n_ox, ppm, tmt=tmt) if badge else None
    if title is None:
        title = f"{peptide}   {charge}+   (scan {scan})"
    ann = plot_psm(spec, peptide, mods, charge, model, title, out_png, subtitle=subtitle,
                   theme=theme, badge=(grade_info["grade"] if grade_info else None), tmt=tmt)
    # observed-vs-predicted correlation for the QC record (mirror already drew it)
    sc = drop_tmt(clean(spec), tmt_kind(tmt, mods), spec["pz"] or charge, ann["tol"][1],
                  keep_mz=[m for _, _, _, m in by_ions(peptide, mods, min(2, max(1, charge)))])
    pm, pit, _ = predict_ms2pip(peptide, mods, charge, model)
    r = (mirror_corr(sc["mz"], sc["it"] / sc["it"].max(), pm, pit / pit.max(),
                     tol=max(0.05, ann["tol"][1]))
         if pm is not None and len(pit) else float("nan"))
    return dict(peptide=peptide, charge=charge, scan=scan, model=model,
                nb=ann["nb"], ny=ann["ny"], bmax=ann["bmax"], ymax=ann["ymax"],
                covfrac=ann["covfrac"], n_ox=n_ox, ppm=ppm,
                highres=ann["tol"][1] < 0.1, tol=ann["tol"][1], ms2pip_r=r,
                grade=(grade_info["grade"] if grade_info else None), out_png=out_png)
