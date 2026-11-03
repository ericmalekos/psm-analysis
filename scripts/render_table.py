#!/usr/bin/env python3
"""Batch-render every PSM in a TSV (one figure per row).

The TSV must have (at least) these columns:
    peptide  charge  scan  spectra_path   [dataset]  [run]  [target]  [tmt]

- ``dataset`` (optional) selects the MS2PIP model + TMT via psmplot/datasets.py.
- ``tmt`` (optional: ``tmt`` or ``tmtpro``) sets the TMT label for that row (overrides ``dataset``).
- ``target`` (optional) groups figures into subfolders and names them.
- rows whose ``spectra_path`` is missing/absent on disk are skipped.

This renders exactly the rows you give it; PSM *selection* (which peptides are the
flagship 3 per protein) is a separate, project-specific concern and is not done here.

Example:
    python scripts/render_table.py \\
        --table ../PSM_Verification/render_qc.tsv \\
        --outdir ./out
"""
import os
import sys
import csv
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import psmplot


def main():
    ap = argparse.ArgumentParser(description="Batch-render PSM figures from a TSV.")
    ap.add_argument("--table", required=True, help="TSV of PSMs to render")
    ap.add_argument("--outdir", default="./out", help="output directory")
    ap.add_argument("--limit", type=int, default=0, help="render at most N rows (0 = all)")
    ap.add_argument("--minimal", action="store_true",
                    help="strip title/subtitle/legend/side-labels; keep ion labels, QC line, axes")
    ap.add_argument("--processes", type=int, default=1, help="MS2PIP 4 batch prediction processes")
    args = ap.parse_args()
    theme = psmplot.MINIMAL if args.minimal else None

    rows = list(csv.DictReader(open(args.table), delimiter="\t"))
    if args.limit:
        rows = rows[:args.limit]

    def params(r):
        pep = r["peptide"]
        charge = int(r["charge"]) if str(r.get("charge", "")).strip() else 2
        dataset = r.get("dataset", "")
        run = r.get("run", "")
        model = psmplot.model_for(dataset, run) if dataset else "HCD2019"
        tmt = (r.get("tmt") or "").strip() or (psmplot.is_tmt(dataset) if dataset else False)
        return pep, charge, dataset, run, model, tmt

    # MS2PIP 4 loads its model files on every single prediction: predict all rows first, one call per model
    items = []
    for r in rows:
        path = r.get("spectra_path", "")
        if path and os.path.exists(path):
            pep, charge, dataset, run, model, tmt = params(r)
            try:
                _, mods, _, _ = psmplot.resolve_psm(pep, charge, path, r["scan"], tmt=tmt)
            except Exception:
                continue
            items.append((pep, mods, charge, model))
    psmplot.predict_ms2pip_batch(items, processes=args.processes)

    n_ok = n_skip = 0
    qc_out = []
    for r in rows:
        path = r.get("spectra_path", "")
        if not path or not os.path.exists(path):
            n_skip += 1
            continue
        pep, charge, dataset, run, model, tmt = params(r)
        target = r.get("target", "") or "PSM"

        sub = os.path.join(args.outdir, target)
        stem = f"{pep[:24]}_{dataset or 'run'}_scan{r['scan']}".replace("/", "_")
        out = os.path.join(sub, stem + ".png")
        title = f"{target}   {pep}   {charge}+   (scan {r['scan']})"
        subtitle = (f"{dataset} · {run.split('.')[0]}" if dataset or run else "")
        try:
            qc = psmplot.render_psm(peptide=pep, charge=charge, spectrum_path=path, scan=r["scan"],
                                    out_png=out, model=model, tmt=tmt, title=title, subtitle=subtitle,
                                    theme=theme)
        except Exception as e:
            print(f"  SKIP {target} {pep} scan {r['scan']}: {e}")
            n_skip += 1
            continue
        n_ok += 1
        qc.update(target=target, dataset=dataset, run=run)
        qc_out.append(qc)
        print(f"  {target:16s} {pep:26s} b{qc['nb']}/y{qc['ny']} "
              f"r={qc['ms2pip_r']:.2f} TIC={qc['covfrac']*100:.0f}% -> {os.path.relpath(out)}")

    if qc_out:
        fields = ["target", "peptide", "charge", "scan", "dataset", "run", "model",
                  "nb", "ny", "bmax", "ymax", "covfrac", "ms2pip_r", "ppm", "n_ox", "highres", "out_png"]
        with open(os.path.join(args.outdir, "render_qc.tsv"), "w") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for q in qc_out:
                w.writerow(q)
    print(f"\nrendered {n_ok} figures, skipped {n_skip} -> {args.outdir}")


if __name__ == "__main__":
    main()
