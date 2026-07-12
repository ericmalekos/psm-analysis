#!/usr/bin/env python3
"""Grade every flagship PSM and write grades.tsv — each rubric metric with its PASS/FAIL.

For each of the 29 flagship PSMs (render_qc.tsv), load the spectrum (on-disk MGF via
PSM_master_table.tsv, or the preserved macrophage PROXI JSON), grade it with psmplot.grade,
and emit a row carrying every metric and whether it passes its "excellent" threshold. The
final `grade` is: caution if a stray-peak trigger fires; else low_information if coverage is
too low; else excellent iff ALL seven *_pass columns are PASS; else good.

Thresholds: cov>=0.90 · gap<=1 · (b2 present OR n_y>=len-1) · |ppm|<=12 (<=13 TMT) · n_ox<=2 ·
base peak explained · max_unexplained<0.50.  (Down-grade triggers: caution = max_unexp>=0.60 &
cov<0.85; low_information = cov<0.55.)
"""
import os
import sys
import csv
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import psmplot

PSMV = "/public/groups/carpenterlab/people/emalekos/sORFs/Proteomics/PSM_Verification"
PF = {True: "PASS", False: "FAIL"}


def main():
    ap = argparse.ArgumentParser(description="Write grades.tsv (every rubric metric + PASS/FAIL).")
    ap.add_argument("--qc", default=f"{PSMV}/render_qc.tsv")
    ap.add_argument("--master", default=f"{PSMV}/PSM_master_table.tsv")
    ap.add_argument("--json-dir", default=f"{PSMV}/spectra_macrophage_json")
    ap.add_argument("--out", default=f"{PSMV}/figures_psmplot/grades.tsv")
    args = ap.parse_args()

    qc = list(csv.DictReader(open(args.qc), delimiter="\t"))
    mt = list(csv.DictReader(open(args.master), delimiter="\t"))
    by3 = {(r["peptide"], r["scan"], r["dataset"]): r for r in mt}
    by2 = {}
    for r in mt:
        by2.setdefault((r["peptide"], r["scan"]), r)

    rows, skipped = [], []
    for q in qc:
        pep, z, ds, run, scan = q["peptide"], int(q["charge"]), q["dataset"], q["run"], q["scan"]
        jpath = os.path.join(args.json_dir, f"{run}_{scan}.json")
        if os.path.exists(jpath):                                   # re-downloaded macrophage scan
            spec = psmplot.load_spectrum(jpath, scan)
            tol = ("Da", 0.02) if spec and spec.get("analyzer") == "FTMS" else ("Da", 0.4)
            tmt = False
        else:                                                       # on-disk MGF
            m = by3.get((pep, scan, ds)) or by2.get((pep, scan))
            path = m["spectra_path"] if m else ""
            spec = psmplot.load_spectrum(path, scan) if path and os.path.exists(path) else None
            tol = None
            tmt = psmplot.is_tmt(ds)
        if spec is None or not len(spec["mz"]):
            skipped.append((q["target"], pep))
            continue
        mods = psmplot.base_mods(pep, tmt=tmt)
        mods, nox, ppm = psmplot.infer_oxidations(pep, mods, spec["pmz"], spec["pz"] or z)
        g = psmplot.grade(spec, pep, mods, z, nox, ppm, tmt=tmt, tol=tol)
        rows.append(dict(
            target=q["target"], peptide=pep, dataset=ds, charge=z, grade=g["grade"],
            ms2pip_r=q.get("ms2pip_r", ""), PEP=q.get("pep", ""), unique_to_target=q.get("unique_to_target", ""),
            cov=g["cov"], cov_pass=PF[g["cov_ok"]],
            max_gap=g["gap"], gap_pass=PF[g["gap_ok"]],
            b2=g["b2"], n_y=g["ny"], pep_len=g["L"], series_pass=PF[g["series_ok"]],
            ppm=g["ppm"], ppm_pass=PF[g["ppm_ok"]],
            n_ox=g["nox"], nox_pass=PF[g["nox_ok"]],
            base_peak_assign=g["base_assign"], base_pass=PF[g["base_ok"]],
            max_unexplained=g["max_unexp"], unexp_pass=PF[g["unexp_ok"]],
            caution_trigger=g["caution_trigger"], lowinfo_trigger=g["lowinfo_trigger"],
            n_b=g["nb"], tol_Da=g["tol"]))

    cols = ["target", "peptide", "dataset", "charge", "grade", "ms2pip_r", "PEP", "unique_to_target",
            "cov", "cov_pass", "max_gap", "gap_pass", "b2", "n_y", "pep_len", "series_pass",
            "ppm", "ppm_pass", "n_ox", "nox_pass", "base_peak_assign", "base_pass",
            "max_unexplained", "unexp_pass", "caution_trigger", "lowinfo_trigger", "n_b", "tol_Da"]
    gord = {"excellent": 0, "good": 1, "caution": 2, "low_information": 3}
    rows.sort(key=lambda r: (gord.get(r["grade"], 9), r["target"], r["peptide"]))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    from collections import Counter
    print(f"wrote {len(rows)} graded PSMs -> {args.out}")
    print("grades:", dict(Counter(r["grade"] for r in rows)))
    if skipped:
        print("skipped (no spectrum):", skipped)


if __name__ == "__main__":
    main()
