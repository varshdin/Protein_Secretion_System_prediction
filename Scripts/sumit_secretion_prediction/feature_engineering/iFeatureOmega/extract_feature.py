"""
Feature extraction for bacterial secretion system classification.

Extracts AAC (amino acid composition) or DPC (dipeptide composition) from a
protein FASTA file. The N-/C-terminal window blocks will be added here later.

Runs in the ifeature_env conda environment (Python 3.7).

Usage:
    python extract_feature.py --block aac --fasta train.fasta
    python extract_feature.py --block dpc --fasta train.fasta

conda create -n ifeature_env python=3.7
conda activate ifeature_env
pip install iFeatureOmegaCLI


All paths are resolved against the repository root, so the script works
from any working directory and on any machine that clones the repo.
"""

import argparse
from pathlib import Path

import numpy as np
from iFeatureOmegaCLI import iProtein


# This file lives at:
#   <repo>/Scripts/sumit_secretion_prediction/feature_engineering/iFeatureOmega/
# so the repository root is four levels up. Anchoring to __file__ rather than
# the current working directory keeps every path below portable.
PROJECT_ROOT = Path(__file__).resolve().parents[4]

FASTA_DIR = (PROJECT_ROOT / "Results" / "sumit_secretion_prediction"
             / "phase0_dataset_prep" / "Train_Test")

FEATURE_DIR = (PROJECT_ROOT / "Results" / "sumit_secretion_prediction"
               / "phase1_feature_engineering" / "features")

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def check_project_root():
    """Fail loudly if PROJECT_ROOT is wrong.

    The parents[4] index depends on where this file sits in the tree. If the
    script is ever moved, every path below silently points somewhere else.
    """
    if not (PROJECT_ROOT / "Results").is_dir():
        raise RuntimeError(
            "repository root looks wrong: %s\n"
            "expected a 'Results' directory there - has this script moved?"
            % PROJECT_ROOT
        )


def read_fasta(path):
    """Read a FASTA file into an ordered {protein_id: sequence} dict."""
    seqs = {}
    pid = None
    with open(str(path)) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                pid = line[1:].strip()
                if pid in seqs:
                    raise ValueError("duplicate id: " + pid)
                seqs[pid] = []
            else:
                seqs[pid].append(line.upper())
    return {pid: "".join(parts) for pid, parts in seqs.items()}


def check_sequences(seqs):
    """Stop if anything other than the 20 standard amino acids is present.

    Composition descriptors divide by the sequence length, so a non-standard
    residue inflates the denominator without filling any column. The row then
    sums to less than 1 and nothing downstream notices.
    """
    observed = set("".join(seqs.values()))
    unexpected = sorted(observed - set(ALPHABET))
    if unexpected:
        bad = [pid for pid, s in seqs.items() if set(s) - set(ALPHABET)]
        raise ValueError(
            "non-standard residues %s in %d sequence(s): %s"
            % (unexpected, len(bad), bad[:10])
        )
    print("input : %d sequences, alphabet OK" % len(seqs))


def run_descriptor(fasta_path, seqs, descriptor, columns):
    """Call iFeatureOmega and return a tidy DataFrame: Protein_id + columns."""
    prot = iProtein(str(fasta_path))
    prot.get_descriptor(descriptor)
    enc = prot.encodings

    # get_descriptor() leaves encodings as None on a bad descriptor name
    # instead of raising, so check here rather than failing later.
    if enc is None:
        raise RuntimeError("no encodings returned for descriptor %r" % descriptor)

    # iFeatureOmega returns rows in FASTA order, indexed by protein id.
    # Realign by id so row order can never matter.
    enc.index.name = "Protein_id"
    enc = enc.reindex(list(seqs.keys()))
    if enc.isna().any().any():
        raise RuntimeError("id mismatch between FASTA and package output")

    return enc[columns].reset_index()


def extract_aac(fasta_path, seqs):
    """Extract the 20-dimensional amino acid composition."""
    columns = ["AAC_" + a for a in ALPHABET]
    return run_descriptor(fasta_path, seqs, "AAC", columns)


def extract_dpc(fasta_path, seqs):
    """Extract the 400-dimensional dipeptide composition."""
    columns = ["DPC_" + a + b for a in ALPHABET for b in ALPHABET]
    return run_descriptor(fasta_path, seqs, "DPC type 1", columns)


BLOCKS = {
    "aac": {"extract": extract_aac, "dim": 20, "prefix": "AAC"},
    "dpc": {"extract": extract_dpc, "dim": 400, "prefix": "DPC"},
}


def check_features(df, n_expected, n_features):
    """Sanity checks before the matrix is written."""
    feat = df.drop(columns="Protein_id")

    assert df.shape == (n_expected, n_features + 1), df.shape
    assert df["Protein_id"].is_unique
    assert df.notna().all().all()
    assert (feat >= 0).all().all()
    assert np.allclose(feat.values.sum(axis=1), 1.0, atol=1e-9)

    print("checks: %d rows x %d features, all passed" % (df.shape[0], n_features))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--block", required=True, choices=sorted(BLOCKS),
        help="feature block to extract"
    )
    parser.add_argument(
        "--fasta", required=True,
        help="FASTA file name inside the Train_Test folder, or a full path"
    )
    parser.add_argument(
        "--out", default=None,
        help="output directory (default: the repo's phase1 features folder)"
    )
    args = parser.parse_args()

    check_project_root()
    block = BLOCKS[args.block]

    # A bare name resolves inside the repo; a full path is used as given, so
    # an unseen set of proteins can be passed from anywhere.
    fasta_path = Path(args.fasta)
    if not fasta_path.is_absolute():
        candidate = FASTA_DIR / fasta_path
        if candidate.exists():
            fasta_path = candidate
    fasta_path = fasta_path.resolve()

    if not fasta_path.exists():
        raise FileNotFoundError("no such FASTA file: %s" % fasta_path)

    out_dir = Path(args.out).resolve() if args.out else FEATURE_DIR / block["prefix"]
    out_dir.mkdir(parents=True, exist_ok=True)

    print("block : %s" % args.block)
    print("fasta : %s" % fasta_path)

    seqs = read_fasta(fasta_path)
    check_sequences(seqs)

    features = block["extract"](fasta_path, seqs)
    check_features(features, len(seqs), block["dim"])

    out_path = out_dir / ("%s_%s.csv" % (block["prefix"], fasta_path.stem))
    features.to_csv(str(out_path), index=False, float_format="%.10f")
    print("wrote : %s" % out_path)


if __name__ == "__main__":
    main()