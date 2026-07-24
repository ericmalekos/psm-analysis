#!/usr/bin/env python3
"""Static (SVG + PNG) structure evidence map — the print/paper version of the interactive viewer.

Renders the predicted fold as a depth-cued Cα "worm" (spline-smoothed backbone tube), each residue
tinted by the strongest PSM covering it — hue = spectral grade, shade = genome uniqueness, grey = no
support. This is the standard way to colour a structure by a per-residue property (like a B-factor
putty), and it exports as clean vector SVG (scalable/editable for the figure) plus a 300-dpi PNG.

Category logic is imported verbatim from build_evidence_map.py, so this matches the interactive map
residue-for-residue and takes the same one-command AlphaFold3 swap:
    render_evidence_static.py SEMR=af3/semr.cif SYNIR=af3/synir.cif

Outputs per locus:  <KEY>_evidence.svg  and  <KEY>_evidence.png   (two orientations + legend).
Run with the ms2rescore python (Biopython + matplotlib).
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
from Bio.PDB import PDBParser, MMCIFParser

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_evidence_map as bem

HERE = os.path.dirname(os.path.abspath(__file__))
WHITE = np.array([1.0, 1.0, 1.0])
GRADE_ORDER = ["excellent", "good", "caution", "low_information"]


def ca_xyz(path):
    fmt = "cif" if path.lower().endswith((".cif", ".mmcif")) else "pdb"
    parser = MMCIFParser(QUIET=True) if fmt == "cif" else PDBParser(QUIET=True)
    st = parser.get_structure("m", path)
    model = next(iter(st))
    chain = max(model, key=lambda c: sum(1 for _ in c))
    return {res.id[1]: np.array(res["CA"].get_coord(), float) for res in chain if "CA" in res}


def catmull_rom(P, n=12):
    """Smooth a polyline P (N,3) through its points. Returns (dens_pts, src_residue_index)."""
    P = np.asarray(P, float)
    N = len(P)
    if N < 3:
        return P.copy(), np.arange(N)
    ext = np.vstack([P[0], P, P[-1]])          # pad both ends so every real segment is interpolated
    out, src = [], []
    for i in range(1, len(ext) - 2):
        p0, p1, p2, p3 = ext[i - 1], ext[i], ext[i + 1], ext[i + 2]
        for t in np.linspace(0.0, 1.0, n, endpoint=False):
            t2, t3 = t * t, t * t * t
            pt = 0.5 * (2 * p1 + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2
                        + (-p0 + 3 * p1 - 3 * p2 + p3) * t3)
            out.append(pt); src.append(i - 1)   # p1 == P[i-1]
    out.append(P[-1]); src.append(N - 1)
    return np.array(out), np.array(src)


def orient(P):
    """PCA frame: PC1=long axis (x), PC2 (y), PC3=depth. Returns (mean, rotation Vt)."""
    mu = P.mean(0)
    _, _, Vt = np.linalg.svd(P - mu, full_matrices=False)
    if np.linalg.det(Vt) < 0:                   # keep right-handed so depth sign is consistent
        Vt[-1] *= -1
    return mu, Vt


def draw_worm(ax, proj, srcidx, res_rgb, ax_i, ax_j, ax_k, sign, lw):
    """One projected view: depth-sorted, depth-shaded per-residue coloured tube."""
    xy = proj[:, [ax_i, ax_j]]
    depth = sign * proj[:, ax_k]
    seg = np.stack([xy[:-1], xy[1:]], axis=1)
    scol = np.array([res_rgb[srcidx[k]] for k in range(len(xy) - 1)])
    sdep = 0.5 * (depth[:-1] + depth[1:])
    if sdep.ptp() > 0:
        t = (sdep - sdep.min()) / sdep.ptp()
    else:
        t = np.ones_like(sdep)
    w = (0.72 + 0.28 * t)[:, None]              # near = full colour, far = faded toward white
    shaded = scol * w + WHITE * (1 - w)
    order = np.argsort(sdep)                     # far first, near drawn on top
    lwv = (lw * (0.85 + 0.30 * t))[order]        # near slightly fatter -> depth
    lc = LineCollection(seg[order], colors=shaded[order], linewidths=lwv,
                        capstyle="round", joinstyle="round", antialiased=True)
    ax.add_collection(lc)
    m = 0.06 * max(xy[:, 0].ptp(), xy[:, 1].ptp())
    ax.set_xlim(xy[:, 0].min() - m, xy[:, 0].max() + m)
    ax.set_ylim(xy[:, 1].min() - m, xy[:, 1].max() + m)
    ax.set_aspect("equal"); ax.axis("off")
    return xy


def render(key, cfg, resinfo, coords):
    resid_cat = {ri["resid"]: ri["cat"] for ri in resinfo.values()}
    ordered = [r for r in sorted(coords) if r in resid_cat]
    P = np.array([coords[r] for r in ordered])
    res_rgb = np.array([to_rgb(bem.COLORS[resid_cat[r]]) if resid_cat[r] else to_rgb(bem.UNCOVERED)
                        for r in ordered])
    dens, src = catmull_rom(P, n=12)
    mu, Vt = orient(P)
    proj = (dens - mu) @ Vt.T

    length = len(cfg["seq"])
    covered = sum(1 for c in resid_cat.values() if c)
    counts = {}
    for c in resid_cat.values():
        if c:
            counts[c] = counts.get(c, 0) + 1

    fig = plt.figure(figsize=(11.2, 6.4), facecolor="white")
    fig.suptitle(f"{cfg['title']}   ·   {cfg['gene']}", x=0.5, y=0.965,
                 fontsize=15, fontweight="bold", ha="center", color="#12191c")
    fig.text(0.5, 0.915, f"predicted fold coloured by PSM support   ·   {cfg['method']}   ·   "
             f"{covered}/{length} residues supported ({round(100*covered/length)}%)",
             ha="center", fontsize=10.5, color="#5c6b6b")

    lw = 6.2
    axA = fig.add_axes([0.02, 0.20, 0.47, 0.66])
    xyA = draw_worm(axA, proj, src, res_rgb, 0, 1, 2, +1, lw)
    axA.text(0.5, 1.02, "front", transform=axA.transAxes, ha="center", fontsize=10, color="#5c6b6b")
    axB = fig.add_axes([0.51, 0.20, 0.47, 0.66])
    draw_worm(axB, proj, src, res_rgb, 0, 2, 1, -1, lw)
    axB.text(0.5, 1.02, "top (rotated 90°)", transform=axB.transAxes, ha="center", fontsize=10, color="#5c6b6b")

    # N / C termini on the front view
    for idx, lab in [(0, "N"), (len(P) - 1, "C")]:
        p = (P[idx] - mu) @ Vt.T
        axA.annotate(lab, (p[0], p[1]), fontsize=11, fontweight="bold", color="#12191c",
                     ha="center", va="center",
                     bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="#12191c", lw=1.1))

    # legend: categories present (canonical order) + uncovered, with residue counts
    handles, labels = [], []
    for g in GRADE_ORDER:
        for u in ("unique", "multiple"):
            ck = f"{g}|{u}"
            if counts.get(ck):
                handles.append(Line2D([0], [0], color=bem.COLORS[ck], lw=7, solid_capstyle="round"))
                labels.append(f"{bem.GRADE_LABEL[g]} · {'unique' if u=='unique' else 'multiple loci'}  ({counts[ck]})")
    handles.append(Line2D([0], [0], color=bem.UNCOVERED, lw=7, solid_capstyle="round"))
    labels.append(f"no PSM support  ({length - covered})")
    leg = fig.legend(handles, labels, loc="lower center", ncol=min(4, len(handles)),
                     frameon=False, fontsize=9.5, handlelength=1.6, columnspacing=1.8,
                     bbox_to_anchor=(0.5, 0.015))
    for txt in leg.get_texts():
        txt.set_color("#12191c")
    fig.text(0.5, 0.11, "hue = spectral grade   ·   dark = unique to locus / light = also at other genomic loci",
             ha="center", fontsize=9, color="#5c6b6b", style="italic")

    svg = os.path.join(HERE, f"{key}_evidence.svg")
    png = os.path.join(HERE, f"{key}_evidence.png")
    fig.savefig(svg, facecolor="white", bbox_inches="tight", pad_inches=0.2)
    fig.savefig(png, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)
    print(f"  [{key}] {covered}/{length} supported -> {os.path.basename(svg)}, {os.path.basename(png)}")
    return svg, png


def render_refined(key, cfg, resinfo, coords):
    """Single clean orientation, PyMOL-free: a Cα PUTTY — tube width AND opacity encode pLDDT, so the
    disordered low-confidence loop recedes and the structured, evidence-bearing core dominates the frame."""
    resid_cat = {ri["resid"]: ri["cat"] for ri in resinfo.values()}
    resid_pl = {ri["resid"]: ri["plddt"] for ri in resinfo.values()}
    ordered = [r for r in sorted(coords) if r in resid_cat]
    P = np.array([coords[r] for r in ordered])
    res_rgb = np.array([to_rgb(bem.COLORS[resid_cat[r]]) if resid_cat[r] else to_rgb(bem.UNCOVERED)
                        for r in ordered])
    pl = np.array([resid_pl[r] for r in ordered], float)
    plo, phi = np.percentile(pl, 5), np.percentile(pl, 95)
    conf = np.clip((pl - plo) / max(phi - plo, 1e-6), 0, 1)
    dens, src = catmull_rom(P, n=12)
    mu, Vt = orient(P)
    proj = (dens - mu) @ Vt.T
    xy = proj[:, [0, 1]]
    depth = proj[:, 2]
    seg = np.stack([xy[:-1], xy[1:]], axis=1)
    ssrc = src[:-1]
    base = res_rgb[ssrc]
    cseg = conf[ssrc]
    sdep = 0.5 * (depth[:-1] + depth[1:])
    t = (sdep - sdep.min()) / (sdep.ptp() or 1)
    wd = (0.74 + 0.26 * t)[:, None]
    c1 = base * wd + WHITE * (1 - wd)                 # depth shade
    wc = (0.5 + 0.5 * cseg)[:, None]
    shaded = c1 * wc + WHITE * (1 - wc)               # low pLDDT -> paler
    lw = (1.3 + 7.4 * cseg) * (0.85 + 0.30 * t)       # low pLDDT -> thinner
    order = np.argsort(sdep)

    length = len(cfg["seq"])
    covered = sum(1 for c in resid_cat.values() if c)
    counts = {}
    for c in resid_cat.values():
        if c:
            counts[c] = counts.get(c, 0) + 1

    fig = plt.figure(figsize=(7.6, 8.0), facecolor="white")
    fig.suptitle(f"{cfg['title']}   ·   {cfg['gene']}", x=0.5, y=0.975, fontsize=14.5,
                 fontweight="bold", color="#12191c")
    fig.text(0.5, 0.935, f"Cα putty — tube width & opacity = pLDDT   ·   {cfg['method']}   ·   "
             f"{covered}/{length} supported ({round(100*covered/length)}%)",
             ha="center", fontsize=10, color="#5c6b6b")
    ax = fig.add_axes([0.03, 0.16, 0.94, 0.74])
    ax.add_collection(LineCollection(seg[order], colors=shaded[order], linewidths=lw[order],
                                     capstyle="round", joinstyle="round", antialiased=True))
    m = 0.06 * max(xy[:, 0].ptp(), xy[:, 1].ptp())
    ax.set_xlim(xy[:, 0].min() - m, xy[:, 0].max() + m)
    ax.set_ylim(xy[:, 1].min() - m, xy[:, 1].max() + m)
    ax.set_aspect("equal"); ax.axis("off")
    for idx, lab in [(0, "N"), (len(P) - 1, "C")]:
        p = (P[idx] - mu) @ Vt.T
        ax.annotate(lab, (p[0], p[1]), fontsize=11, fontweight="bold", color="#12191c",
                    ha="center", va="center",
                    bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="#12191c", lw=1.1))
    handles, labels = [], []
    for g in GRADE_ORDER:
        for u in ("unique", "multiple"):
            ck = f"{g}|{u}"
            if counts.get(ck):
                handles.append(Line2D([0], [0], color=bem.COLORS[ck], lw=7, solid_capstyle="round"))
                labels.append(f"{bem.GRADE_LABEL[g]} · {'unique' if u=='unique' else 'multiple loci'}  ({counts[ck]})")
    handles.append(Line2D([0], [0], color=bem.UNCOVERED, lw=7, solid_capstyle="round"))
    labels.append(f"no PSM support  ({length - covered})")
    leg = fig.legend(handles, labels, loc="lower center", ncol=min(4, len(handles)), frameon=False,
                     fontsize=9.5, handlelength=1.6, columnspacing=1.8, bbox_to_anchor=(0.5, 0.02))
    for txt in leg.get_texts():
        txt.set_color("#12191c")
    svg = os.path.join(HERE, f"{key}_putty.svg")
    png = os.path.join(HERE, f"{key}_putty.png")
    fig.savefig(svg, facecolor="white", bbox_inches="tight", pad_inches=0.2)
    fig.savefig(png, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)
    print(f"  [{key}] putty -> {os.path.basename(svg)}, {os.path.basename(png)}")


def main():
    overrides = {}
    for a in sys.argv[1:]:
        if "=" in a:
            k, p = a.split("=", 1)
            overrides[k.strip().upper()] = p.strip()
    evidence = bem.load_evidence()
    done = 0
    for key, cfg in bem.LOCI.items():
        cfg = dict(cfg)
        if key in overrides:
            cfg["model"] = overrides[key]
            if overrides[key].lower().endswith((".cif", ".mmcif")):
                cfg["method"] = "AlphaFold3"
        if not cfg["model"] or not os.path.exists(cfg["model"]):
            print(f"[{key}] no model yet — skipping (provide {key}=/path/to/model.cif)")
            continue
        if key not in evidence:
            print(f"[{key}] no peptides — skipping"); continue
        _, chain_res, _, _ = bem.parse_structure(cfg["model"])
        model_seq = "".join(a for _, a, _ in chain_res)
        offset = bem.align_offset(model_seq, cfg["seq"], key)
        resinfo, _ = bem.residue_categories(evidence[key], cfg["seq"], chain_res, offset)
        coords = ca_xyz(cfg["model"])
        render(key, cfg, resinfo, coords)
        render_refined(key, cfg, resinfo, coords)
        done += 1
    if not done:
        print("Nothing rendered.")


if __name__ == "__main__":
    main()
