# psm-analysis

Annotate peptide–spectrum matches (PSMs) and mirror them against MS2PIP predictions, for
manual PSM verification. Each figure has three panels: a b/y **fragment-coverage map**, the
**annotated observed** MS2, and an **observed-vs-MS2PIP mirror**. Optional rubric **grade
badge** (excellent / good / caution / low-information). Import package: `psmplot`.

<img src="docs/example_psm.png" alt="example PSM figure" width="56.25%">


## Install

Python ≥ 3.8 with **numpy, matplotlib, pyteomics, ms2pip** (ms2pip does the fragment-intensity
prediction: MS2PIP 3 when installed, else MS2PIP 4, which gets every modification as a mass shift;
`render_table.py` predicts all rows in one MS2PIP 4 call per model, `--processes N`;
`$PSMPLOT_MS2PIP_MODEL_DIR` overrides MS2PIP 4's model directory). No install needed — the scripts add the repo to `sys.path`; `pip install -e .`
also provides a `plot-psm` command.

## Generate a figure

**Input:** peptide sequence, charge, a spectrum file (**MGF / mzML / PROXI-JSON**), scan number.

```bash
python scripts/plot_psm.py \
    --peptide LMDLSHTNSVAELQETPGK --charge 3 \
    --spectrum run.mgf --scan 30662 --out fig.png     # optional: --minimal  --badge
```

**Output:** `fig.png` (+ `fig.pdf`) and a QC line (b/y counts, coverage, precursor ppm, MS2PIP
r, grade). Batch from a TSV (columns `peptide charge scan spectra_path [dataset] [run] [target]`):

```bash
python scripts/render_table.py --table psms.tsv --outdir out
```

## TMT / TMTpro data

`--tmt` (TMT 6/10/11-plex, 229.1629 Da) or `--tmt tmtpro` (TMTpro 16/18-plex, 304.2071 Da) puts the label on the
N-term and every K; in `render_table.py` the optional `tmt` column does the same per row, and `render_psm`,
`plot_psm` and `grade` take `tmt=` (`True`, `"tmt"` or `"tmtpro"`). For a labelled peptide the ions that carry no
sequence information are removed from the spectrum before annotation, plotting and grading, as the intact precursor
is: the reporter ions (126 to 131 for TMT, 126 to 135 for TMTpro), the intact tag ion (230.170 / 305.214) and the
complement ions (precursor less one reporter ion and CO, at charge z-1; in high-resolution spectra only), each with
its 13C isotope, and the tag + H2O ion. A peak that is also one of the peptide's b/y ions is kept. Otherwise a
reporter base peak counts as an unexplained peak in the rubric and drives most TMT spectra to *caution*. The QC line
notes the removal.
`python tests/test_tmt.py` checks this.

## PSM source pipeline

PSMs identified with **MS-GF+**, rescored with **MS²Rescore v2 / Percolator**; high-confidence
= **q < 0.01 AND PEP < 0.1**. Modifications: fixed Carbamidomethyl-C, variable Oxidation-M
(inferred from the precursor mass), TMT-6plex where applicable; MS2PIP model chosen per
acquisition (HCD2019 / CID / TMT / Immuno-HCD).

Restyle everything from `psmplot/theme.py` (one `Theme` dataclass — palette, size, fonts, badge colors).
