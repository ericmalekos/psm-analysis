#!/usr/bin/env bash
# Render one known-good PSM (Gvin-ps7, iMAC proteomics) as a smoke test / design sample.
set -euo pipefail
PY=/public/groups/carpenterlab/people/emalekos/bin/ms2rescore/ms2rescore/bin/python3
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(dirname "$HERE")"
DATA=/public/groups/carpenterlab/people/emalekos/sORFs/Proteomics

"$PY" "$REPO/scripts/plot_psm.py" \
    --peptide LMDLSHTNSVAELQETPGK --charge 3 \
    --spectrum "$DATA/iMACProteomics/20220110_EMalekos_Carpenter_12000_1_WT.mgf" \
    --scan 30662 --dataset iMACProteomics \
    --out "$HERE/out/gvin_ps7_LMDLSHTNSVAELQETPGK.png" \
    --title "Gvin-ps7   LMDLSHTNSVAELQETPGK   3+"
