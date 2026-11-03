"""Figure construction: fragment-coverage map + annotated observed spectrum + observed-vs-MS2PIP mirror.

Three stacked panels:
  1. a fragment-coverage map (peptide sequence with b/y cleavage flags),
  2. the observed MS2 spectrum with b/y ions annotated,
  3. the same spectrum mirrored against the MS2PIP prediction.

All visual values are read from a :class:`~psmplot.theme.Theme`; this module contains
layout logic only. ``plot_psm`` is the low-level entry point (spectrum + resolved mods);
most callers use ``psmplot.render_psm``, which resolves mods and charge for you.
"""
import os
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .annotate import annotate, clean, longest_series, detect_tol, by_ions, tmt_kind, drop_tmt
from .prediction import predict_ms2pip, mirror_corr
from .theme import DEFAULT


def _ion_label(ion):
    """Human ion label, e.g. b7, y3, y5 with a superscript 2+ for doubly charged."""
    t, num, c, _ = ion
    return f"{t}{num}" + ("" if c == 1 else "²⁺" if c == 2 else f"({c}+)")


def _coverage_map(ax, pep, mods, ann, th):
    """Draw the peptide as a fragment-coverage map: sequence letters with b/y cleavage flags.

    Convention: a cleavage between residues i and i+1 yields b_(i+1) (N-terminal) and
    y_(n-i-1) (C-terminal). An observed b-ion draws a blue flag *below* the sequence
    with its foot pointing left (toward the N-terminus it contains); an observed y-ion
    draws a red flag *above* with its foot pointing right (toward the C-terminus).
    """
    n = len(pep)
    have_b = {ion[1] for _, ion, _, _ in ann["matched"] if ion[0] == "b"}
    have_y = {ion[1] for _, ion, _, _ in ann["matched"] if ion[0] == "y"}
    ax.set_xlim(-0.7, n - 0.3)
    ax.set_ylim(0, 1)
    ax.axis("off")
    base = 0.5
    for i, a in enumerate(pep):
        modded = i in mods and a in "CM"
        ax.text(i, base, a, fontsize=th.map_fontsize, ha="center", va="center",
                family="monospace", weight="bold",
                color=th.mod_residue_color if modded else "#222222")
    top_y = base + th.map_gap + th.map_tick_len
    bot_y = base - th.map_gap - th.map_tick_len
    for k in range(1, n):                      # cleavage after residue k -> b_k / y_(n-k)
        xb = k - 0.5
        if k in have_b:                        # blue flag below, foot toward N-terminus (left)
            ax.plot([xb, xb], [base - th.map_gap, bot_y], color=th.b_color, lw=th.map_lw, solid_capstyle="round")
            ax.plot([xb, xb - th.map_flag_len], [bot_y, bot_y], color=th.b_color, lw=th.map_lw, solid_capstyle="round")
        if (n - k) in have_y:                  # red flag above, foot toward C-terminus (right)
            ax.plot([xb, xb], [base + th.map_gap, top_y], color=th.y_color, lw=th.map_lw, solid_capstyle="round")
            ax.plot([xb, xb + th.map_flag_len], [top_y, top_y], color=th.y_color, lw=th.map_lw, solid_capstyle="round")


def plot_psm(spec, pep, mods, z, model, title, out_png, subtitle="", tol=None, theme=None, badge=None,
             tmt=None):
    """Render one PSM to ``out_png`` (+ a .pdf if the theme asks). Returns the annotation dict.

    ``badge`` (optional) is a rubric grade string ("excellent"/"good"/"caution"/
    "low_information"); when given, a colored grade badge is drawn in the spectrum panel.
    ``tmt`` (True, ``"tmt"`` or ``"tmtpro"``) removes the TMT reporter, tag and complement ions
    from the plotted spectrum, as the intact precursor is removed.
    """
    th = theme or DEFAULT
    spec = clean(spec)                                    # drop the intact precursor
    kind = tmt_kind(tmt, mods)
    if kind:                                              # and the TMT reporter, tag and complement ions
        if tol is None:
            tol = detect_tol(spec["mz"], by_ions(pep, mods, min(2, max(1, z))))
        spec = drop_tmt(spec, kind, spec.get("pz") or z, tol[1],
                        keep_mz=[m for _, _, _, m in by_ions(pep, mods, min(2, max(1, z)))])
    if not len(spec["mz"]):
        raise ValueError("no fragment peaks left after removing the precursor" + (" and TMT ions" if kind else ""))
    ann = annotate(spec, pep, mods, maxch=min(2, max(1, z)), tol=tol,
                   min_rel_int=th.annotate_min_rel_int)
    mz, it = spec["mz"], spec["it"]
    itn = it / it.max()
    bmax, ymax = longest_series(pep, ann["matched"])
    nb = len({ion[1] for _, ion, _, _ in ann["matched"] if ion[0] == "b"})
    ny = len({ion[1] for _, ion, _, _ in ann["matched"] if ion[0] == "y"})
    tolstr = f"{ann['tol'][1]} Da ({'high' if ann['tol'][1] < 0.1 else 'low'}-res)"
    qc = (f"{nb} b-ions, {ny} y-ions   ·   longest run b{bmax}/y{ymax}   ·   "
          f"{ann['covfrac'] * 100:.0f}% of fragment TIC annotated   ·   tol {tolstr}"
          + (f"   ·   {'TMTpro' if kind == 'tmtpro' else 'TMT'} ions removed" if kind else ""))

    fig = plt.figure(figsize=th.figsize)
    gs = fig.add_gridspec(3, 1, height_ratios=list(th.height_ratios), hspace=th.hspace)
    ax_map = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])
    matched_idx = {j: ion for j, ion, _, _ in ann["matched"]}

    # ---- header: title, the scan line under it, then the QC line (stacked from the top) ----
    # Callers pass the scan/run inside a trailing "(… scan N)" on the title; lift it onto its own
    # smaller line directly under the title rather than leaving it in the (now smaller) title.
    scan_line = ""
    if title:
        m = re.search(r"\s*\(([^)]*scan[^)]*)\)\s*$", title, re.I)
        if m:
            scan_line, title = m.group(1).strip(), title[:m.start()].rstrip()
    y, drew = 0.990, False
    if th.show_title and title:
        fig.text(0.5, y, title, ha="center", va="top", fontsize=th.title_fontsize,
                 weight=th.title_weight, linespacing=th.title_linespacing)
        y -= 0.024
        drew = True
    if th.show_title and scan_line:
        fig.text(0.5, y, scan_line, ha="center", va="top",
                 fontsize=th.title_sub_fontsize, color=th.title_sub_color)
        y -= 0.019
        drew = True
    if th.show_qc:
        fig.text(0.5, y, qc, ha="center", va="top", fontsize=th.qc_fontsize, color=th.qc_color)
        y -= 0.018
        drew = True
    top = (y - 0.004) if drew else 0.975

    # ---- row 0: fragment-coverage map ----
    _coverage_map(ax_map, pep, mods, ann, th)

    # ---- row 1: annotated observed spectrum ----
    xlim = (min(th.xlim_left_floor, mz.min()) - th.xlim_pad_left, mz.max() * th.xlim_pad_right)
    ax1.vlines(mz, 0, itn, color=th.unmatched_color, lw=th.peak_lw_unmatched, zorder=1)
    for j, ion in matched_idx.items():
        col = th.b_color if ion[0] == "b" else th.y_color
        ax1.vlines(mz[j], 0, itn[j], color=col, lw=th.peak_lw_matched, zorder=3)
        if itn[j] > th.frag_label_min_rel_int:
            ax1.text(mz[j], itn[j] + th.frag_label_offset, _ion_label(ion), color=col,
                     fontsize=th.frag_label_fontsize, ha="center", va="bottom",
                     rotation=th.frag_label_rotation, zorder=4)
    ax1.set_ylim(0, th.top_ymax)
    ax1.set_xlim(*xlim)
    if th.show_axis_labels:
        ax1.set_ylabel("relative intensity")
    if th.show_legend:
        ax1.text(0.995, 0.99, "b ions", transform=ax1.transAxes, ha="right", va="top",
                 color=th.b_color, fontsize=th.legend_fontsize, weight="bold")
        ax1.text(0.995, 0.93, "y ions", transform=ax1.transAxes, ha="right", va="top",
                 color=th.y_color, fontsize=th.legend_fontsize, weight="bold")
    if badge:
        bcol = th.badge_colors.get(badge, "#888888")
        ax1.text(0.012, 0.985, badge.replace("_", " ").upper(), transform=ax1.transAxes,
                 ha="left", va="top", fontsize=th.badge_fontsize, weight="bold", color="white",
                 bbox=dict(boxstyle="round,pad=0.35", facecolor=bcol, edgecolor="none"))

    # ---- row 2: mirror vs MS2PIP ----
    pmz, pit, pann = predict_ms2pip(pep, mods, z, model)
    ax2.axhline(0, color="k", lw=th.peak_lw_unmatched)
    ax2.vlines(mz, 0, itn, color=th.unmatched_color, lw=th.peak_lw_unmatched)
    for j, ion in matched_idx.items():
        col = th.b_color if ion[0] == "b" else th.y_color
        ax2.vlines(mz[j], 0, itn[j], color=col, lw=th.peak_lw_matched)
        if th.mirror_label_observed and itn[j] > th.frag_label_min_rel_int:
            ax2.text(mz[j], itn[j] + th.frag_label_offset, _ion_label(ion), color=col,
                     fontsize=th.frag_label_fontsize, ha="center", va="bottom",
                     rotation=th.frag_label_rotation)
    if pmz is not None and len(pit):
        pn = pit / pit.max()
        for m_, i_, a_ in zip(pmz, pn, pann):
            col = th.b_color if a_.startswith("b") else th.y_color if a_.startswith("y") else th.unmatched_color
            ax2.vlines(m_, 0, -i_, color=col, lw=th.peak_lw_matched)
            if th.mirror_label_predicted and i_ > th.frag_label_min_rel_int:
                ax2.text(m_, -i_ - th.frag_label_offset, a_, color=col,
                         fontsize=th.frag_label_fontsize, ha="center", va="top",
                         rotation=th.frag_label_rotation)
        if th.show_corr:
            r = mirror_corr(mz, itn, pmz, pn, tol=max(0.05, ann["tol"][1]))
            ax2.text(0.995, 0.06, f"MS2PIP ({model}) Pearson r = {r:.2f}", transform=ax2.transAxes,
                     ha="right", va="bottom", fontsize=th.corr_fontsize, color=th.corr_text_color,
                     weight="bold")
    elif th.show_corr:
        ax2.text(0.5, -0.5, "MS2PIP prediction unavailable", ha="center", color="#999")
    ax2.set_ylim(-th.mirror_ymax, th.mirror_ymax)
    ax2.set_xlim(*xlim)
    if th.show_axis_labels:
        ax2.set_xlabel("m/z")
        ax2.set_ylabel("relative intensity")
    if th.show_mirror_side_labels:
        ax2.text(0.005, 0.97, th.observed_label, transform=ax2.transAxes, ha="left", va="top",
                 fontsize=th.legend_fontsize, weight="bold", color="#444")
        ax2.text(0.005, 0.03, th.predicted_label, transform=ax2.transAxes, ha="left", va="bottom",
                 fontsize=th.legend_fontsize, weight="bold", color="#444")

    for ax in (ax1, ax2):
        ax.spines[["top", "right"]].set_visible(False)
    if th.show_subtitle and subtitle:
        bottom = 0.075
        fig.text(0.5, 0.015, subtitle, ha="center", va="bottom",
                 fontsize=th.subtitle_fontsize, color=th.subtitle_color)
    else:
        bottom = 0.062
    fig.subplots_adjust(left=0.11, right=0.97, top=top, bottom=bottom)
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    plt.savefig(out_png, dpi=th.dpi)
    if th.save_pdf and out_png.lower().endswith(".png"):
        plt.savefig(out_png[:-4] + ".pdf")
    plt.close()
    ann["nb"], ann["ny"], ann["bmax"], ann["ymax"] = nb, ny, bmax, ymax
    return ann
