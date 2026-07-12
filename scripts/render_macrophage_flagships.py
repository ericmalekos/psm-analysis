#!/usr/bin/env python3
"""Render the re-downloaded macrophage-lysate flagships (the render_qc.tsv rows whose iProX
peak lists were purged) through psmplot, in full and minimal styles.

Spectra come from ThermoRawFileParser `query` PROXI JSONs in --json-dir, named
``<run>_<scan>.json`` (as produced by dl_flagship_macs.sh). For each PSM we set the
fragment tolerance from the analyzer (FTMS -> 0.02 Da, ion trap -> 0.4 Da) and pick the
better-correlating MS2PIP model between HCD2019 and CID, matching the original macrophage
figures. PSM metadata (score/q/PEP/unique/ppm) is read from render_qc.tsv.
"""
import os
import sys
import csv
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import psmplot

PSMV = "/public/groups/carpenterlab/people/emalekos/sORFs/Proteomics/PSM_Verification"


def main():
    ap = argparse.ArgumentParser(description="Render re-downloaded macrophage flagships (full + minimal).")
    ap.add_argument("--qc", default=f"{PSMV}/render_qc.tsv")
    ap.add_argument("--json-dir", default=f"{PSMV}/spectra_macrophage_json",
                    help="dir of <run>_<scan>.json PROXI files (the preserved re-downloaded scans)")
    ap.add_argument("--outdir", default=f"{PSMV}/figures_psmplot")
    ap.add_argument("--styles", default="full,minimal")
    ap.add_argument("--badge", action="store_true",
                    help="also render a 'badged' series (full style + rubric grade badge)")
    args = ap.parse_args()

    qc = list(csv.DictReader(open(args.qc), delimiter="\t"))
    styles = [s.strip() for s in args.styles.split(",") if s.strip()]
    if args.badge and "badged" not in styles:
        styles.append("badged")
    theme_of = {"full": None, "minimal": psmplot.MINIMAL, "badged": None}

    done, missing = [], []
    for q in qc:
        jpath = os.path.join(args.json_dir, f"{q['run']}_{q['scan']}.json")
        if not os.path.exists(jpath):
            continue                                        # not one of the re-downloaded macrophage rows
        pep, z, ds = q["peptide"], int(q["charge"]), q["dataset"]
        spec = psmplot.load_spectrum(jpath, q["scan"])
        if spec is None or not len(spec["mz"]):
            missing.append(q)
            continue
        mods = psmplot.base_mods(pep, tmt=False)
        mods, nox, ppm = psmplot.infer_oxidations(pep, mods, spec["pmz"], spec["pz"] or z)
        tolval = 0.02 if spec.get("analyzer") == "FTMS" else 0.4
        sc = psmplot.clean(spec)
        base = sc["it"].max() if len(sc["it"]) else 1.0
        best = None                                          # pick the better-correlating model
        for model in ("HCD2019", "CID"):
            pm, pit, _ = psmplot.predict_ms2pip(pep, mods, z, model)
            if pm is None or not len(pit):
                continue
            rr = psmplot.mirror_corr(sc["mz"], sc["it"] / base, pm, pit / pit.max(), tol=max(0.05, tolval))
            if best is None or (rr == rr and rr > best[1]):
                best = (model, rr)
        model = best[0] if best else "HCD2019"

        share = "unique to target" if q.get("unique_to_target") == "True" else f"shared w/ {q.get('shared_with','')}"
        try:
            scores = f"score {float(q['score']):.2f}, q={float(q['qvalue']):.1e}, PEP={float(q['pep']):.1e}"
        except Exception:
            scores = ""
        oxs = f" · {nox} Ox-Met" if nox else ""
        short_share = "unique" if q.get("unique_to_target") == "True" else "shared"
        title = f"{q['target']}   {pep}   {z}+   ({q['run']} scan {q['scan']})"
        subtitle = (f"{ds}/{q.get('db','')} · {scores} · {short_share} · "
                    f"macrophage ({spec.get('analyzer','?')}) · {ppm:+.1f} ppm{oxs}")
        g = psmplot.grade(spec, pep, mods, z, nox, ppm, tmt=False, tol=("Da", tolval))
        gr = g["grade"]
        stem = f"rank{q['rank']}_{pep[:24]}_{ds}".replace("/", "_")
        for style in styles:
            out = os.path.join(args.outdir, style, q["target"], stem + ".png")
            psmplot.plot_psm(spec, pep, mods, z, model, title, out, subtitle=subtitle,
                             tol=("Da", tolval), theme=theme_of[style],
                             badge=(gr if style == "badged" else None))
        q["_grade"] = gr
        done.append(q)
        r_disp = f"{best[1]:.2f}" if best else "n/a"
        print(f"  [ok] {q['target']:22s} {pep[:26]:26s} {spec.get('analyzer','?'):5s} model={model:8s} r={r_disp} grade={gr}")

    print(f"\nrendered {len(done)} macrophage flagships in styles [{', '.join(styles)}] -> {args.outdir}")
    from collections import Counter
    print("grades:", dict(Counter(q["_grade"] for q in done)))
    exc = [q for q in done if q["_grade"] == "excellent"]
    print(f"excellent ({len(exc)}):" + ("" if exc else " none"))
    for q in exc:
        print(f"   {q['target']:22s} {q['peptide']}")
    if missing:
        print(f"{len(missing)} had a JSON but failed to load:")
        for q in missing:
            print(f"   {q['target']} {q['peptide']} scan {q['scan']}")


if __name__ == "__main__":
    main()
