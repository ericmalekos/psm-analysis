"""Project-specific mapping from a dataset (and run name) to its MS2PIP model and label chemistry.

This is the one module that encodes *this study's* datasets. If you reuse psmplot for
another project, this is what you change. The MS2PIP model choice must match how each
dataset was acquired (HCD vs CID vs TMT vs immunopeptidome-HCD).
"""

# dataset id -> MS2PIP model, with a couple of run-name-dependent special cases below.
_DATASET_MODEL = {
    "PXD029155": "TMT",                 # TMT 6-plex whole-proteome
    "PXD001905": "CID",                 # secretome, low-res CID
    "iMACProteomics": "CID",
}


def model_for(dataset, run=""):
    """Return the MS2PIP model name for a PSM from ``dataset`` (``run`` disambiguates PXD008733)."""
    if dataset == "PXD029155":
        return "TMT"
    if dataset == "PXD008733_price":
        return "Immuno-HCD" if "HCD" in run else "CID"
    if dataset in ("PXD001905", "iMACProteomics"):
        return "CID"
    return "HCD2019"


def is_tmt(dataset):
    """True if PSMs from ``dataset`` carry TMT labels (N-term + K)."""
    return dataset == "PXD029155"
