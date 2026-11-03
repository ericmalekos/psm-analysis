#!/usr/bin/env python3
"""Render a single PSM to an annotated + MS2PIP-mirror figure.

Example:
    python scripts/plot_psm.py \\
        --peptide LMDLSHTNSVAELQETPGK --charge 3 \\
        --spectrum /path/run.mgf --scan 30662 \\
        --dataset iMACProteomics --out fig.png \\
        --title "Gvin-ps7  LMDLSHTNSVAELQETPGK  3+"

If --model is omitted it is derived from --dataset (see psmplot/datasets.py).
"""
import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import psmplot


def main():
    ap = argparse.ArgumentParser(description="Render one PSM (annotated observed + MS2PIP mirror).")
    ap.add_argument("--peptide", required=True, help="bare amino-acid sequence (no mods)")
    ap.add_argument("--charge", type=int, required=True)
    ap.add_argument("--spectrum", required=True, help="MGF or mzML(.gz) file")
    ap.add_argument("--scan", required=True, help="scan number")
    ap.add_argument("--out", default="psm.png", help="output PNG path (a .pdf is written alongside)")
    ap.add_argument("--dataset", default="", help="dataset id -> MS2PIP model + TMT (datasets.py)")
    ap.add_argument("--model", default="", help="MS2PIP model override (e.g. HCD2019, CID, TMT, Immuno-HCD)")
    ap.add_argument("--run", default="", help="run name (only needed to disambiguate PXD008733)")
    ap.add_argument("--tmt", nargs="?", const="tmt", default="", choices=["tmt", "tmtpro"],
                    help="force TMT labelling (N-term + K): --tmt or --tmt tmt for TMT 6/10/11-plex, "
                         "--tmt tmtpro for TMTpro 16/18-plex; reporter, tag and complement ions are removed")
    ap.add_argument("--title", default=None)
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--minimal", action="store_true",
                    help="strip title/subtitle/legend/side-labels; keep ion labels, QC line, axes")
    ap.add_argument("--badge", action="store_true", help="compute + draw the rubric grade badge")
    args = ap.parse_args()

    model = args.model or (psmplot.model_for(args.dataset, args.run) if args.dataset else "HCD2019")
    tmt = args.tmt or (psmplot.is_tmt(args.dataset) if args.dataset else False)
    theme = psmplot.MINIMAL if args.minimal else None

    qc = psmplot.render_psm(
        peptide=args.peptide, charge=args.charge, spectrum_path=args.spectrum, scan=args.scan,
        out_png=args.out, model=model, tmt=tmt, title=args.title, subtitle=args.subtitle,
        theme=theme, badge=args.badge)

    print(f"wrote {qc['out_png']}")
    print(f"  model={model} tmt={tmt}  {'HIGH' if qc['highres'] else 'low'}-res (tol {qc['tol']} Da)"
          + (f"  grade={qc['grade']}" if qc.get('grade') else ""))
    print(f"  b{qc['nb']}/y{qc['ny']}  longest run b{qc['bmax']}/y{qc['ymax']}  "
          f"fragTIC={qc['covfrac']*100:.0f}%  MS2PIP r={qc['ms2pip_r']:.2f}  "
          f"precursor {qc['ppm']:+.1f} ppm  {qc['n_ox']} Ox-Met")


if __name__ == "__main__":
    main()
