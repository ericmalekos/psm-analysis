#!/usr/bin/env python3
"""Render a predicted structure coloured by PSM support — publication CARTOON (+ optional surface and
2-D secondary-structure track). Every path and important option is a CLI flag, so this is the reusable
generator for the AlphaFold3 models (drop-in: point --model at the .cif).

Colouring: each residue takes its strongest peptide's 2-D category — spectral grade (hue) x genome
uniqueness (dark=unique / light=multiple loci); grey = no confident PSM. Evidence comes from the combined
PSM table; category logic is shared with build_evidence_map.py. PyMOL ray-traces on CPU (no GPU/display).

Examples
--------
  # SEMR ESMFold preview (registry default model), cartoon only:
  render_pymol.py

  # AlphaFold3 models, all three representations, into a chosen folder:
  render_pymol.py --model SEMR=af3/semr.cif --model SYNIR=af3/synir.cif \
                  --representations cartoon,surface,ss-track --outdir figures/

  # tune: keep every residue, hi-res ray, custom pymol:
  render_pymol.py SEMR=af3/semr.cif --no-trim --ray 3000x2400 --pymol /opt/pymol/bin/pymol
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_evidence_map as bem

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PYMOL = "/public/groups/carpenterlab/people/emalekos/conda/pymol/bin/pymol"
DEFAULT_EVIDENCE = os.path.normpath(bem.EVIDENCE_TSV)
GRADE_ORDER = ["excellent", "good", "caution", "low_information"]
REPS = ("cartoon", "surface", "ss-track")


class PipelineError(Exception):
    """A per-locus failure that should be reported but not crash the whole run."""


# ----------------------------------------------------------------- helpers
def hex_rgb(h):
    return [int(h[1:3], 16) / 255, int(h[3:5], 16) / 255, int(h[5:7], 16) / 255]


def ranges(ints):
    """[1,2,3,7,8] -> '1-3+7-8' (compact PyMOL resi selection)."""
    ints = sorted(ints)
    out, i = [], 0
    while i < len(ints):
        j = i
        while j + 1 < len(ints) and ints[j + 1] == ints[j] + 1:
            j += 1
        out.append(f"{ints[i]}-{ints[j]}" if j > i else f"{ints[i]}")
        i = j + 1
    return "+".join(out)


def resolve_pymol(arg):
    for cand in (arg, DEFAULT_PYMOL, shutil.which("pymol")):
        if cand and os.path.exists(cand) and os.access(cand, os.X_OK):
            return cand
    raise SystemExit(
        f"error: PyMOL executable not found (tried: {arg or '-'}, {DEFAULT_PYMOL}, $PATH).\n"
        f"       install with:  micromamba create -p <prefix>/pymol -c conda-forge pymol-open-source\n"
        f"       then pass --pymol <prefix>/pymol/bin/pymol")


def parse_model_args(tokens):
    """['SEMR=a.cif', 'SYNIR=b.cif'] -> {'SEMR':'a.cif','SYNIR':'b.cif'}."""
    out = {}
    for t in tokens or []:
        if "=" not in t:
            raise SystemExit(f"error: --model/positional must be LOCUS=PATH, got {t!r}")
        k, p = t.split("=", 1)
        out[k.strip().upper()] = os.path.abspath(os.path.expanduser(p.strip()))
    return out


def trim_threshold(args, plddts):
    if args.no_trim:
        return None
    if args.trim_plddt is not None:
        return args.trim_plddt
    return 0.5 if max(plddts) <= 1.5 else 50.0          # auto: ESMFold 0-1 vs AlphaFold 0-100


# ----------------------------------------------------------------- pymol script
def build_pml(key, cfg, resinfo, palette, outdir, thresh, ray_w, ray_h, rt_mode, reps):
    runs = bem.contiguous_runs(sorted((ri["resid"], ri["cat"]) for ri in resinfo.values()))
    n_hidden = 0
    want_3d = ("cartoon" in reps) or ("surface" in reps)
    model = os.path.abspath(cfg["model"])
    ss_tsv = os.path.join(outdir, f"{key}_ss.tsv")
    cart_png = os.path.join(outdir, f"{key}_cartoon_raw.png")
    surf_png = os.path.join(outdir, f"{key}_surface_raw.png")
    L = ["python", f"cmd.load({model!r}, 'obj')", "python end",
         "hide everything", "bg_color white", "set ray_opaque_background, 0",
         "set cartoon_fancy_helices, 1", "set cartoon_smooth_loops, 1",
         "set cartoon_highlight_color, grey60", "set cartoon_transparency, 0",
         "set ambient, 0.42", "set direct, 0.55", "set specular, 0.25", "set ray_shadows, 1",
         "set antialias, 2", f"set ray_trace_mode, {rt_mode}", "dss"]
    # SS dump first (full length, before any removal)
    L += ["python", "ssd = []",
          "cmd.iterate('obj and name CA', 'ssd.append((int(resi), ss))', space={'ssd': ssd})",
          f"open({ss_tsv!r},'w').write('\\n'.join('%d\\t%s' % (r, s or 'L') for r, s in sorted(ssd)))",
          "python end"]
    if want_3d:
        L += ["show cartoon", "color grey70, obj"]
        for cat in palette:
            sel = "+".join(f"{s}-{e}" for cc, s, e in runs if cc == cat)
            if sel:
                r, g, b = hex_rgb(palette[cat])
                nm = "c_" + cat.replace("|", "_")
                L += [f"set_color {nm}, [{r:.3f},{g:.3f},{b:.3f}]", f"color {nm}, obj and resi {sel}"]
        if cfg.get("boundary"):
            L += [f"color black, obj and resi {cfg['boundary']}-{cfg['boundary']+1}"]
        if thresh is not None:
            hide = sorted(ri["resid"] for ri in resinfo.values()
                          if ri["cat"] is None and ri["plddt"] < thresh)
            n_hidden = len(hide)
            if hide:
                L += [f"remove obj and resi {ranges(hide)}"]
        L += ["orient obj"]
        if "cartoon" in reps:
            L += ["python", f"cmd.ray({ray_w}, {ray_h})", f"cmd.png({cart_png!r}, dpi=300)", "python end"]
        if "surface" in reps:
            L += ["hide cartoon", "set surface_quality, 1", "set transparency, 0", "show surface, obj",
                  "python", f"cmd.ray({ray_w}, {ray_h})", f"cmd.png({surf_png!r}, dpi=300)", "python end"]
    pml = os.path.join(outdir, f"{key}_pymol.pml")
    with open(pml, "w") as fo:
        fo.write("\n".join(L) + "\n")
    return pml, n_hidden, ss_tsv, cart_png, surf_png


def run_pymol(pymol, pml, timeout):
    try:
        r = subprocess.run([pymol, "-cq", pml], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise PipelineError(f"PyMOL timed out after {timeout}s")
    if r.returncode != 0:
        raise PipelineError("PyMOL exited %d:\n%s" % (r.returncode, (r.stderr or r.stdout)[-1200:]))
    return r


# ----------------------------------------------------------------- figures
def cat_legend(resinfo, length, covered, palette, ss=False):
    counts = {}
    for ri in resinfo.values():
        if ri["cat"]:
            counts[ri["cat"]] = counts.get(ri["cat"], 0) + 1
    handles, labels = [], []
    for g in GRADE_ORDER:
        for u in ("unique", "multiple"):
            ck = f"{g}|{u}"
            if counts.get(ck):
                handles.append(Rectangle((0, 0), 1, 1, facecolor=palette[ck]))
                labels.append(f"{bem.GRADE_LABEL[g]} · {'unique' if u=='unique' else 'multiple loci'} ({counts[ck]})")
    handles.append(Rectangle((0, 0), 1, 1, facecolor=bem.UNCOVERED))
    labels.append(f"no PSM support ({length - covered})")
    if ss:
        handles += [Line2D([0], [0], color="#39403f", lw=6), Line2D([0], [0], color="#39403f", lw=3),
                    Line2D([0], [0], color="#39403f", lw=1)]
        labels += ["helix", "strand", "loop"]
    return handles, labels


def compose_labeled(raw, out, cfg, resinfo, kind, n_hidden, palette, dpi, transparent):
    img = plt.imread(raw)
    asp = img.shape[0] / img.shape[1]
    fw = 8.4
    ih = fw * asp
    Ht = ih + 1.75
    fig = plt.figure(figsize=(fw, Ht), facecolor="none" if transparent else "white")
    length = len(cfg["seq"])
    covered = sum(1 for ri in resinfo.values() if ri["cat"])
    ax = fig.add_axes([0.02, 1.0 / Ht, 0.96, ih / Ht])
    ax.imshow(img)
    ax.axis("off")
    fig.text(0.5, 1 - 0.30 / Ht, f"{cfg['title']}   ·   {cfg['gene']}", ha="center",
             fontsize=14, fontweight="bold", color="#12191c")
    sub = (f"{kind} coloured by PSM support   ·   {cfg['method']}   ·   "
           f"{covered}/{length} supported ({round(100*covered/length)}%)")
    if n_hidden:
        sub += f"   ·   {n_hidden} disordered unsupported residues hidden"
    fig.text(0.5, 1 - 0.60 / Ht, sub, ha="center", fontsize=9.5, color="#5c6b6b")
    handles, labels = cat_legend(resinfo, length, covered, palette)
    leg = fig.legend(handles, labels, loc="lower center", ncol=min(4, len(handles)), frameon=False,
                     fontsize=9, handlelength=1.4, columnspacing=1.6, bbox_to_anchor=(0.5, 0.05 / Ht))
    for t in leg.get_texts():
        t.set_color("#12191c")
    fig.savefig(out, dpi=dpi, facecolor="none" if transparent else "white", transparent=transparent)
    plt.close(fig)


def draw_ss_track(key, cfg, resinfo, ss_by_resid, palette, outdir, dpi):
    resid_cat = {ri["resid"]: ri["cat"] for ri in resinfo.values()}
    resid_orf = {ri["resid"]: pos for pos, ri in resinfo.items()}
    resids = sorted(resid_cat)
    xs = {r: resid_orf[r] for r in resids}

    def sstype(r):
        s = ss_by_resid.get(r, "L")
        return "H" if s == "H" else ("S" if s in ("S", "E") else "L")

    HGT = {"H": 1.0, "S": 0.60, "L": 0.18}
    runs = []
    for r in resids:
        st = sstype(r)
        if runs and runs[-1][0] == st and xs[r] == runs[-1][2] + 1:
            runs[-1][2] = xs[r]; runs[-1][3].append(r)
        else:
            runs.append([st, xs[r], xs[r], [r]])

    length = len(cfg["seq"])
    covered = sum(1 for c in resid_cat.values() if c)
    x0, x1 = min(xs.values()), max(xs.values())
    fig = plt.figure(figsize=(max(8.5, (x1 - x0) / 26.0), 2.75), facecolor="white")
    fig.suptitle(f"{cfg['title']}   ·   {cfg['gene']}", x=0.5, y=0.98, fontsize=13.5,
                 fontweight="bold", color="#12191c")
    fig.text(0.5, 0.885, f"secondary structure coloured by PSM support   ·   {cfg['method']}   ·   "
             f"{covered}/{length} residues supported ({round(100*covered/length)}%)",
             ha="center", fontsize=9.5, color="#5c6b6b")
    ax = fig.add_axes([0.035, 0.30, 0.94, 0.46])

    def col(r):
        c = resid_cat[r]
        return palette[c] if c else bem.UNCOVERED

    for st, a, b, rr in runs:
        h = HGT[st]
        for r in rr:
            ax.add_patch(Rectangle((xs[r] - 0.5, -h / 2), 1.0, h, facecolor=col(r), edgecolor="none", zorder=2))
        if st == "H":
            ax.add_patch(Rectangle((a - 0.5, -h / 2), (b - a + 1), h, fill=False, edgecolor="#39403f", lw=0.7, zorder=3))
        elif st == "S":
            ax.add_patch(Polygon([(b + 0.5, -HGT["H"] / 2), (b + 0.5 + 2.2, 0), (b + 0.5, HGT["H"] / 2)],
                                 closed=True, facecolor=col(rr[-1]), edgecolor="#39403f", lw=0.6, zorder=2))
            ax.add_patch(Rectangle((a - 0.5, -h / 2), (b - a + 1), h, fill=False, edgecolor="#39403f", lw=0.6, zorder=3))

    ax.set_xlim(x0 - 2, x1 + 5)
    ax.set_ylim(-0.85, 0.85)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#c9d3cf")
    ax.set_xticks(list(range(((x0 // 20) + 1) * 20, x1 + 1, 20)))
    ax.tick_params(axis="x", labelsize=8, colors="#5c6b6b", length=3)
    for lab, xx in [("N", x0), ("C", x1)]:
        ax.annotate(lab, (xx, 0), fontsize=10, fontweight="bold", ha="center", va="center",
                    bbox=dict(boxstyle="circle,pad=0.12", fc="white", ec="#12191c", lw=1.0), zorder=5)
    handles, labels = cat_legend(resinfo, length, covered, palette, ss=True)
    leg = fig.legend(handles, labels, loc="lower center", ncol=min(5, len(handles)), frameon=False,
                     fontsize=8.6, handlelength=1.5, columnspacing=1.5, bbox_to_anchor=(0.5, 0.0))
    for t in leg.get_texts():
        t.set_color("#12191c")
    outs = []
    for ext in ("svg", "png"):
        p = os.path.join(outdir, f"{key}_ss_track.{ext}")
        fig.savefig(p, facecolor="white", bbox_inches="tight", pad_inches=0.15, dpi=dpi)
        outs.append(p)
    plt.close(fig)
    return outs


# ----------------------------------------------------------------- per-locus driver
def process_locus(key, cfg, evidence, args, palette, pymol):
    if key not in evidence:
        raise PipelineError(f"no peptides for {key} in {args.evidence}")
    model = cfg["model"]
    if not model:
        raise PipelineError(f"no model — pass {key}=/path/to/model.(pdb|cif)")
    if not os.path.exists(model):
        raise PipelineError(f"model not found: {model}")
    if not model.lower().endswith((".pdb", ".cif", ".mmcif")):
        raise PipelineError(f"model must be .pdb/.cif: {model}")

    _, chain_res, _, _ = bem.parse_structure(model)
    if not chain_res:
        raise PipelineError(f"no residues parsed from {model}")
    model_seq = "".join(a for _, a, _ in chain_res)
    exact = (model_seq == cfg["seq"]) or (model_seq in cfg["seq"])
    if not exact and not args.allow_mismatch:
        raise PipelineError(
            f"model sequence does not match the {key} ORF ({len(model_seq)} vs {len(cfg['seq'])} aa). "
            f"Numbering-dependent colouring would be wrong. Re-fold the exact ORF, or pass "
            f"--allow-mismatch to force best-offset alignment.")
    offset = bem.align_offset(model_seq, cfg["seq"], key)
    resinfo, _ = bem.residue_categories(evidence[key], cfg["seq"], chain_res, offset)
    plddts = [ri["plddt"] for ri in resinfo.values()]
    thresh = trim_threshold(args, plddts)

    reps = args.representations
    ray_w, ray_h = args.ray
    pml, n_hidden, ss_tsv, cart_raw, surf_raw = build_pml(
        key, cfg, resinfo, palette, args.outdir, thresh, ray_w, ray_h, args.ray_trace_mode, reps)
    for stale in (cart_raw, surf_raw, ss_tsv):        # a leftover must never mask a failed render
        if os.path.exists(stale):
            os.remove(stale)
    run_pymol(pymol, pml, args.pymol_timeout)

    made = []
    if "cartoon" in reps:
        if not os.path.exists(cart_raw):
            raise PipelineError("PyMOL produced no cartoon PNG (see the .pml)")
        out = os.path.join(args.outdir, f"{key}_cartoon.png")
        compose_labeled(cart_raw, out, cfg, resinfo, "cartoon", n_hidden, palette, args.dpi, args.transparent)
        made.append(out)
    if "surface" in reps:
        if not os.path.exists(surf_raw):
            raise PipelineError("PyMOL produced no surface PNG")
        out = os.path.join(args.outdir, f"{key}_surface.png")
        compose_labeled(surf_raw, out, cfg, resinfo, "surface", n_hidden, palette, args.dpi, args.transparent)
        made.append(out)
    if "ss-track" in reps:
        if not os.path.exists(ss_tsv):
            raise PipelineError("PyMOL produced no SS dump")
        ss_by_resid = {int(a): s for a, s in (ln.rstrip("\n").split("\t")
                       for ln in open(ss_tsv) if ln.strip())}
        made += draw_ss_track(key, cfg, resinfo, ss_by_resid, palette, args.outdir, args.dpi)

    if not args.keep_raw:
        for p in (cart_raw, surf_raw, ss_tsv, pml):
            if os.path.exists(p):
                os.remove(p)
    covered = sum(1 for ri in resinfo.values() if ri["cat"])
    return dict(key=key, made=made, offset=offset, covered=covered, length=len(cfg["seq"]),
                hidden=n_hidden, exact=exact)


# ----------------------------------------------------------------- CLI
def build_parser():
    p = argparse.ArgumentParser(
        description="Render a predicted structure coloured by PSM support (cartoon / surface / SS-track).",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("models", nargs="*", metavar="LOCUS=PATH",
                   help="model file per locus, e.g. SEMR=af3/semr.cif (also accepted via --model)")
    p.add_argument("--model", action="append", default=[], metavar="LOCUS=PATH",
                   help="same as positional LOCUS=PATH; repeatable")
    p.add_argument("--evidence", default=DEFAULT_EVIDENCE, metavar="TSV",
                   help="combined PSM table (default: %(default)s)")
    p.add_argument("--pymol", default=None, metavar="BIN",
                   help=f"pymol executable (default: {DEFAULT_PYMOL} or $PATH)")
    p.add_argument("--outdir", default=SCRIPT_DIR, metavar="DIR", help="output directory (default: script dir)")
    p.add_argument("--loci", default=None, metavar="LIST",
                   help="comma-separated subset to render (default: every locus that has a model)")
    p.add_argument("--representations", default="cartoon", metavar="LIST",
                   help="comma-separated: cartoon,surface,ss-track (default: cartoon)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--trim-plddt", type=float, default=None, metavar="X",
                   help="hide UNSUPPORTED residues below this pLDDT from 3-D views "
                        "(default: auto 0.5 for 0-1 scale, 50 for 0-100)")
    g.add_argument("--no-trim", action="store_true", help="keep every residue in the 3-D views")
    p.add_argument("--ray", default="2000x1600", metavar="WxH", help="ray-trace resolution (default: %(default)s)")
    p.add_argument("--ray-trace-mode", type=int, default=0, choices=(0, 1, 2, 3),
                   help="PyMOL ray_trace_mode (0 smooth, 1 outline; default: 0)")
    p.add_argument("--dpi", type=int, default=300, help="output figure dpi (default: %(default)s)")
    p.add_argument("--transparent", action="store_true", help="transparent background on cartoon/surface PNGs")
    p.add_argument("--palette", default=None, metavar="JSON",
                   help='JSON file overriding category colours, e.g. {"excellent|unique":"#08306b"}')
    p.add_argument("--config", default=None, metavar="JSON",
                   help="JSON adding/overriding loci: {NAME:{seq,boundary,title,gene,note,method,landmarks}}")
    p.add_argument("--pymol-timeout", type=int, default=600, metavar="SEC", help="per-locus PyMOL timeout")
    p.add_argument("--keep-raw", action="store_true", help="keep intermediate .pml / raw PNG / SS dump")
    p.add_argument("--allow-mismatch", action="store_true",
                   help="proceed even if the model sequence does not match the ORF (best-offset align)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    # representations
    reps = [r.strip() for r in args.representations.split(",") if r.strip()]
    bad = [r for r in reps if r not in REPS]
    if bad:
        raise SystemExit(f"error: unknown representation(s) {bad}; choose from {', '.join(REPS)}")
    args.representations = reps

    # ray resolution
    try:
        w, h = args.ray.lower().split("x")
        args.ray = (int(w), int(h))
    except Exception:
        raise SystemExit(f"error: --ray must be WxH (e.g. 2000x1600), got {args.ray!r}")

    # paths
    args.evidence = os.path.abspath(os.path.expanduser(args.evidence))
    if not os.path.exists(args.evidence):
        raise SystemExit(f"error: --evidence not found: {args.evidence}")
    args.outdir = os.path.abspath(os.path.expanduser(args.outdir))
    os.makedirs(args.outdir, exist_ok=True)
    pymol = resolve_pymol(args.pymol)

    # palette
    palette = dict(bem.COLORS)
    if args.palette:
        with open(args.palette) as fh:
            over = json.load(fh)
        unknown = [k for k in over if k not in palette]
        if unknown:
            raise SystemExit(f"error: --palette has unknown categories {unknown}; valid: {list(palette)}")
        palette.update(over)

    # locus registry (+ optional config overrides), then model overrides
    loci = {k: dict(v) for k, v in bem.LOCI.items()}
    if args.config:
        with open(args.config) as fh:
            for name, over in json.load(fh).items():
                loci.setdefault(name.upper(), {}).update(over)
    overrides = parse_model_args(args.models + args.model)
    for name, path in overrides.items():
        if name not in loci:
            raise SystemExit(f"error: unknown locus {name!r}; known: {list(loci)} (add via --config)")
        loci[name]["model"] = path
        if path.lower().endswith((".cif", ".mmcif")):
            loci[name]["method"] = "AlphaFold3"

    # which loci
    if args.loci:
        want = [x.strip().upper() for x in args.loci.split(",") if x.strip()]
        missing = [x for x in want if x not in loci]
        if missing:
            raise SystemExit(f"error: --loci names unknown loci {missing}; known: {list(loci)}")
    else:
        want = [k for k, v in loci.items()
                if k in overrides or (v.get("model") and os.path.exists(v["model"]))]
    if not want:
        raise SystemExit("error: no loci to render — pass a model, e.g. SEMR=/path/model.cif")

    evidence = bem.load_evidence(args.evidence)

    print(f"pymol   : {pymol}")
    print(f"evidence: {args.evidence}")
    print(f"outdir  : {args.outdir}")
    print(f"reps    : {', '.join(args.representations)}   loci: {', '.join(want)}\n")

    ok, failed = [], []
    for key in want:
        try:
            res = process_locus(key, loci[key], evidence, args, palette, pymol)
        except PipelineError as e:
            failed.append((key, str(e)))
            print(f"[{key}] FAILED: {e}\n")
            continue
        ok.append(res)
        note = "" if res["exact"] else "  (forced align — sequence mismatch!)"
        print(f"[{key}] {res['covered']}/{res['length']} supported, offset {res['offset']}, "
              f"hid {res['hidden']} disordered{note}")
        for m in res["made"]:
            print(f"        {os.path.relpath(m, args.outdir)}")
        if args.verbose:
            print()

    print(f"\ndone: {len(ok)} rendered, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
