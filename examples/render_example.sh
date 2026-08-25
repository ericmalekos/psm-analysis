#!/usr/bin/env bash
# Regenerate the README display image (docs/example_psm.png): the SYNIR peptide
# SPGSESQGWGYGLQPK from the secretome (PXD001905, run 20140910_PP_JB_120min_neat_L3,
# scan 17235) — a complete b13/y13 ladder, MS2PIP CID r=0.86. Also a smoke test.
set -euo pipefail
PY=/public/groups/carpenterlab/people/emalekos/bin/ms2rescore/ms2rescore/bin/python3
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(dirname "$HERE")"
DATA=/public/groups/carpenterlab/people/emalekos/sORFs/Proteomics

"$PY" "$REPO/scripts/plot_psm.py" \
    --peptide SPGSESQGWGYGLQPK --charge 2 \
    --spectrum "$DATA/Secretome/PXD001905/20140910_PP_JB_120min_neat_L3.mgf" \
    --scan 17235 --dataset PXD001905 \
    --out "$REPO/docs/example_psm.png" \
    --title "SYNIR   SPGSESQGWGYGLQPK   2+   (20140910_PP_JB_120min_neat_L3 scan 17235)"
