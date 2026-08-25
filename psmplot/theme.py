"""All visual choices in one place.

Everything that controls how a PSM figure *looks* — palette, figure size, fonts,
axis headroom, the fragment-coverage map geometry, label thresholds — is a field on
:class:`Theme`. To restyle the plots, edit the defaults here or pass a customized
``Theme`` to ``plot_psm(...)`` / ``render_psm(...)``. The plotting code in ``plot.py``
reads every magic number from this object, so nothing visual is hard-coded there.
"""
from dataclasses import dataclass, replace, field
from typing import Tuple


@dataclass
class Theme:
    # ---------- palette ----------
    b_color: str = "#1f6fd6"           # b-ion blue (N-terminal fragments)
    y_color: str = "#d1341f"           # y-ion red  (C-terminal fragments)
    unmatched_color: str = "#b8b8b8"   # noise / unannotated peaks
    corr_text_color: str = "#222222"
    subtitle_color: str = "#555555"
    mod_residue_color: str = "#7b3fbf"  # modified residue letter in the coverage map

    # ---------- figure ----------
    figsize: Tuple[float, float] = (7.0, 9.0)   # compact portrait; compressed horizontally
    dpi: int = 150
    save_pdf: bool = True                        # also write a vector .pdf next to the .png
    # three stacked rows: [fragment-coverage map, observed spectrum, mirror]
    height_ratios: Tuple[float, float, float] = (1.5, 6.0, 6.0)
    hspace: float = 0.32

    # ---------- peaks ----------
    peak_lw_unmatched: float = 0.6
    peak_lw_matched: float = 1.2

    # ---------- fragment labels (spectrum panel) ----------
    frag_label_fontsize: float = 6.3
    frag_label_min_rel_int: float = 0.03   # only label peaks taller than this (rel. to base)
    frag_label_offset: float = 0.012       # vertical gap between peak top and its label
    frag_label_rotation: float = 90.0

    # ---------- header text ----------
    title_fontsize: float = 9.0            # smaller than before (was 11)
    title_weight: str = "bold"
    title_linespacing: float = 1.4
    title_sub_fontsize: float = 8.0        # the scan line, small, directly under the title
    title_sub_color: str = "#555555"
    qc_fontsize: float = 8.0               # the QC summary line under the scan line
    qc_color: str = "#333333"

    # ---------- axis limits ----------
    top_ymax: float = 1.20            # spectrum-panel headroom for rotated peak labels
    mirror_ymax: float = 1.30         # mirror-panel half-height (room for ion labels on the predicted side)
    xlim_left_floor: float = 150.0    # left edge = min(this, data.min) - xlim_pad_left
    xlim_pad_left: float = 20.0
    xlim_pad_right: float = 1.02      # right edge = data.max * this

    # ---------- fragment-coverage map (top strip) ----------
    map_fontsize: float = 11.0        # sequence letters
    map_lw: float = 1.5               # flag line width
    map_flag_len: float = 0.36        # horizontal foot length (in residue-width units)
    map_tick_len: float = 0.24        # how far the vertical flag extends from the baseline
    map_gap: float = 0.06             # gap between the letter row and a flag's vertical start

    # ---------- legend + inline annotations ----------
    legend_fontsize: float = 8.5
    corr_fontsize: float = 9.0
    subtitle_fontsize: float = 7.5
    observed_label: str = "observed"          # upper half of the mirror
    predicted_label: str = "predicted (MS2PIP)"  # lower half; set to "theoretical" if you prefer
    mirror_label_predicted: bool = True       # annotate the predicted (lower) peaks with ion labels
    mirror_label_observed: bool = False       # also annotate the observed (upper) peaks (duplicates panel 1)

    # ---------- text visibility (independent toggles; MINIMAL flips the "chrome" ones off) ----------
    show_title: bool = True
    show_subtitle: bool = True
    show_qc: bool = True                  # the QC summary line under the title
    show_legend: bool = True              # the "b ions" / "y ions" key
    show_mirror_side_labels: bool = True  # the "observed" / "predicted (MS2PIP)" half-labels
    show_corr: bool = True                # the Pearson-r annotation
    show_axis_labels: bool = True         # x/y axis titles (tick numbers are always shown)

    # ---------- rubric grade badge (drawn only when a grade is supplied) ----------
    badge_fontsize: float = 8.5
    badge_colors: dict = field(default_factory=lambda: {
        "excellent": "#1a8f3c", "good": "#2b6cb0",
        "caution": "#c05621", "low_information": "#718096"})

    # ---------- annotation noise floor ----------
    annotate_min_rel_int: float = 0.02  # peaks below this (rel. to base) are never matched

    def variant(self, **kw):
        """Return a copy with some fields overridden: ``theme.variant(dpi=300)``."""
        return replace(self, **kw)


DEFAULT = Theme()

# Minimal-text preset: strip the title, subtitle, b/y legend key, the
# observed/predicted side-labels, and the Pearson-r annotation; KEEP the b/y peak
# labels, the QC line, and the axis labels + ticks. For dropping panels into a
# figure where you add your own captions. (To keep the r: MINIMAL.variant(show_corr=True).)
MINIMAL = DEFAULT.variant(show_title=False, show_subtitle=False, show_legend=False,
                          show_mirror_side_labels=False, show_corr=False)
