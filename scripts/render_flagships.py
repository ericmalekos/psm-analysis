#!/usr/bin/env python3
"""Regenerate the paper's flagship PSM figures through psmplot, in full and minimal styles.

The flagship selection lives in PSM_Verification/render_qc.tsv (target, rank, peptide,
charge, dataset, run, scan + Percolator score/q/PEP + unique/shared + ppm). That table
carries everything except the spectrum path, which we resolve from PSM_master_table.tsv.

For each flagship whose spectrum is on disk we render two figures:
  <outdir>/full/<target>/rank<k>_<peptide>_<dataset>.png     (title, subtitle, legend, r)
  <outdir>/minimal/<target>/rank<k>_<peptide>_<dataset>.png  (--minimal: ion labels + QC + axes)

Flagships whose spectra were purged from iProX (the *_macrophage variants, etc.) are
reported as needing re-download rather than silently skipped.
"""
import os
import sys
import csv
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import psmplot

PSMV = "/public/groups/carpenterlab/people/emalekos/sORFs/Proteomics/PSM_Verification"


def main():
    ap = argparse.ArgumentParser(description="Regenerate flagship PSM figures (full + minimal).")
    ap.add_argument("--qc", default=f"{PSMV}/render_qc.tsv", help="flagship selection TSV")
    ap.add_argument("--master", default=f"{PSMV}/PSM_master_table.tsv", help="table with spectra_path")
    ap.add_argument("--outdir", default=f"{PSMV}/figures_psmplot", help="output root")
    ap.add_argument("--styles", default="full,minimal", help="comma list: full,minimal")
    ap.add_argument("--badge", action="store_true",
                    help="also render a 'badged' series (full style + rubric grade badge)")
    args = ap.parse_args()

    qc = list(csv.DictReader(open(args.qc), delimiter="\t"))
    mt = list(csv.DictReader(open(args.master), delimiter="\t"))
    by3 = {(r["peptide"], r["scan"], r["dataset"]): r for r in mt}
    by2 = {}
    for r in mt:
        by2.setdefault((r["peptide"], r["scan"]), r)

    styles = [s.strip() for s in args.styles.split(",") if s.strip()]
    if args.badge and "badged" not in styles:
        styles.append("badged")
    theme_of = {"full": None, "minimal": psmplot.MINIMAL, "badged": None}

    done, missing = [], []
    for q in qc:
        m = by3.get((q["peptide"], q["scan"], q["dataset"])) or by2.get((q["peptide"], q["scan"]))
        path = m["spectra_path"] if m else ""
        if not path or not os.path.exists(path):
            missing.append(q)
            continue
        pep, z, ds, run = q["peptide"], int(q["charge"]), q["dataset"], q["run"]
        model, tmt = psmplot.model_for(ds, run), psmplot.is_tmt(ds)
        try:
            spec = psmplot.load_spectrum(path, q["scan"])
            if spec is None or not len(spec["mz"]):
                missing.append(q)
                continue
            mods = psmplot.base_mods(pep, tmt=tmt)
            mods, nox, ppm = psmplot.infer_oxidations(pep, mods, spec["pmz"], spec["pz"] or z)
        except Exception as e:
            print(f"  LOAD FAIL {q['target']} {pep}: {e}")
            missing.append(q)
            continue

        share = "unique to target" if q.get("unique_to_target") == "True" else f"shared w/ {q.get('shared_with','')}"
        try:
            scores = f"score {float(q['score']):.2f}, q={float(q['qvalue']):.1e}, PEP={float(q['pep']):.1e}"
        except Exception:
            scores = ""
        oxs = f" · {nox} Ox-Met" if nox else ""
        title = f"{q['target']}   {pep}   {z}+   ({run.split('.')[0]} scan {q['scan']})"
        subtitle = f"{ds}/{q.get('db','')} · {scores} · {share} · {ppm:+.1f} ppm{oxs}"
        g = psmplot.grade(spec, pep, mods, z, nox, ppm, tmt=tmt)
        gr = g["grade"]
        stem = f"rank{q['rank']}_{pep[:24]}_{ds}".replace("/", "_")
        for style in styles:
            out = os.path.join(args.outdir, style, q["target"], stem + ".png")
            psmplot.plot_psm(spec, pep, mods, z, model, title, out, subtitle=subtitle,
                             theme=theme_of[style], badge=(gr if style == "badged" else None))
        q["_grade"] = gr
        done.append(q)
        print(f"  [ok] {q['target']:22s} {pep[:26]:26s} {ds:14s} {model:10s} r={q.get('ms2pip_r','?'):>5} grade={gr}")

    print(f"\nrendered {len(done)}/{len(qc)} flagships in styles [{', '.join(styles)}] -> {args.outdir}")
    from collections import Counter
    print("grades:", dict(Counter(q["_grade"] for q in done)))
    exc = [q for q in done if q["_grade"] == "excellent"]
    print(f"excellent ({len(exc)}):" + ("" if exc else " none"))
    for q in exc:
        print(f"   {q['target']:22s} {q['peptide']}")
    if missing:
        print(f"\n{len(missing)} flagships need iProX re-download (spectra not on disk):")
        for q in missing:
            print(f"   {q['target']:22s} {q['peptide'][:26]:26s} z{q['charge']} scan {q['scan']:>7}  {q['dataset']}")


if __name__ == "__main__":
    main()
