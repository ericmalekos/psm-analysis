# structure_evidence

Paint a predicted protein fold by **PSM support**: every residue is coloured by the strongest
peptide-spectrum match that covers it — **hue = spectral grade** (excellent / good / caution /
low-information), **shade = genome uniqueness** (dark = unique to the locus, light = the peptide also
occurs at other genomic loci). Grey = the fold is predicted but no confident PSM covers it.

Built for the noncanonical ORFs SEMR/Brip1os and SYNIR/Env (reviewer-requested manual PSM verification).
The structure source is decoupled from the evidence: residue numbering is validated against the ORF
sequence, so an **AlphaFold3** model drops in with no code change.

## Scripts

| script | output | notes |
|---|---|---|
| `render_pymol.py` | **cartoon** PNG (+ optional surface, 2-D SS-track SVG) | hardened CLI; needs PyMOL |
| `render_evidence_static.py` | Cα worm + putty | vector SVG + PNG |
| `build_evidence_map.py` | interactive 3Dmol viewer | self-contained HTML |

All three share one category logic and take the same `LOCUS=PATH` model override.

## Inputs (tracked)

- `data/combined_peptide_table.tsv` — best-of-all-PSMs grade × genome-uniqueness per peptide. A **snapshot**;
  regenerate from `../../PSM_Verification/genome_uniqueness/combined_table.py`, or point `--evidence` at the live file.
- `models/SEMR_esmfold.pdb` — ESMFold preview model for SEMR (SYNIR awaits its AlphaFold3 model).
- `models/3Dmol-min.js` — vendored viewer for the interactive HTML.

Generated figures / HTML / per-residue TSVs are git-ignored.

## Run

Use a Python with Biopython + matplotlib (this project: the ms2rescore env); PyMOL (open-source) is needed
for the cartoon/surface/SS-track.

```bash
render_pymol.py                                    # SEMR ESMFold cartoon (registry default)
render_pymol.py --model SEMR=af3/semr.cif \
                --model SYNIR=af3/synir.cif \
                --representations cartoon,surface,ss-track --outdir figures/
render_pymol.py --help                             # every path & option is a flag
```

## Colour scheme

grade × uniqueness, dark = unique / light = multiple loci:

| grade | unique | multiple loci |
|---|---|---|
| excellent | `#124e8c` | `#7fb0da` |
| good | `#1c7d47` | `#8fca9f` |
| caution | `#c0561a` | `#eaa06a` |
| low-information | `#6a4b9c` | `#b4a3d6` |

no PSM support = `#9aa4a0`. Override with `render_pymol.py --palette palette.json`.
