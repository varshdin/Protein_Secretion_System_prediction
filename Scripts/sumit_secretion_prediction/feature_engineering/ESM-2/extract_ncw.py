"""
N-/C-terminal window feature extraction using ESM-2.

Cuts the first w_n and last w_c residues of each protein, concatenates them,
and embeds the result with a pre-trained ESM-2 model. The window width was
selected as 50/50 by stratified 5-fold cross-validation on the training set
(SVC, macro-F1); see the window-selection experiment.

Needs torch and transformers, which require Python 3.10+. This is why the
block lives in its own script: iFeatureOmega (aac, dpc) requires Python 3.7
and the two cannot share an environment.

    conda activate esm_env
    python extract_ncw.py --fasta train.fasta
    python extract_ncw.py --fasta test.fasta

Output contract, shared with extract_feature.py so the blocks compose:
    column 0      Protein_id
    columns 1..n  NCW_0 .. NCW_{dim-1}
    row order     same as the input FASTA
    path          features/N_C_window/N_C_{w_n}_window_esm2_{fasta_stem}.csv


conda create -n esm_env python=3.11
conda activate esm_env
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install transformers pandas

All paths resolve against the repository root, so the script runs from any
working directory on any machine that clones the repo.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModel


PROJECT_ROOT = Path(__file__).resolve().parents[4]

FASTA_DIR = (PROJECT_ROOT / "Results" / "sumit_secretion_prediction"
             / "phase0_dataset_prep" / "Train_Test")

FEATURE_DIR = (PROJECT_ROOT / "Results" / "sumit_secretion_prediction"
               / "phase1_feature_engineering" / "features")

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

PREFIX = "NCW"
DEFAULT_MODEL = "facebook/esm2_t6_8M_UR50D"
DEFAULT_WN = 50
DEFAULT_WC = 50
BATCH = 32


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

    ESM-2 maps unknown residues to <unk>, which carries no biological meaning
    but still contributes a vector to the pooled average. Rejecting them here
    keeps this block consistent with the composition blocks, which reject them
    for a different reason.
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


def nc_window(seq, w_n, w_c):
    """First w_n residues plus the last w_c residues, concatenated.

    Sequences shorter than w_n + w_c share residues between the two halves.
    That is reported but allowed: every window is still exactly w_n + w_c long,
    and the uniform length is what lets a batch run without padding.
    """
    return seq[:w_n] + seq[-w_c:]


def embed(windows, tokenizer, model, width):
    """Embed equal-length windows. No padding, so no attention mask needed."""
    out = []
    for start in range(0, len(windows), BATCH):
        batch = tokenizer(windows[start:start + BATCH], return_tensors="pt")
        n_tokens = batch["input_ids"].shape[1]
        if n_tokens != width + 2:
            raise RuntimeError(
                "expected %d tokens (%d residues + cls + eos), got %d"
                % (width + 2, width, n_tokens)
            )
        with torch.no_grad():
            hidden = model(**batch).last_hidden_state
        out.append(hidden[:, 1:-1, :].mean(1))
    return torch.cat(out).numpy().astype(np.float32)


def extract(seqs, model_name, w_n, w_c):
    """Return (DataFrame, metadata) for the N-/C-window embedding."""
    width = w_n + w_c
    pids = list(seqs.keys())
    lengths = np.array([len(seqs[p]) for p in pids])

    too_short = lengths < max(w_n, w_c)
    if too_short.any():
        bad = [pids[i] for i in np.flatnonzero(too_short)]
        raise ValueError(
            "%d sequence(s) shorter than %d residues: %s"
            % (len(bad), max(w_n, w_c), bad[:10])
        )

    overlap = lengths < width
    print("window: N=%d C=%d -> %d residues" % (w_n, w_c, width))
    print("        %d of %d shorter than %d, so their halves share residues "
          "(%.2f%%)" % (overlap.sum(), len(pids), width, 100.0 * overlap.mean()))

    print("model : %s" % model_name)
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name, add_pooling_layer=False)
    model.eval()              
    dim = model.config.hidden_size
    print("        %d dimensions, %d layers"
          % (dim, model.config.num_hidden_layers))

    windows = [nc_window(seqs[p], w_n, w_c) for p in pids]
    observed = set(len(w) for w in windows)
    if observed != {width}:
        raise RuntimeError("window lengths not uniform: %s" % sorted(observed))

    matrix = embed(windows, tokenizer, model, width)

    alone = embed(windows[:1], tokenizer, model, width)[0]
    together = embed(windows[:BATCH], tokenizer, model, width)[0]
    drift = float(np.abs(alone - together).max())
    if drift > 1e-4:
        raise RuntimeError("embedding depends on batch composition (%.2e)" % drift)
    print("checks: batch invariance %.2e, OK" % drift)

    df = pd.DataFrame(matrix, columns=["%s_%d" % (PREFIX, i) for i in range(dim)])
    df.insert(0, "Protein_id", pids)

    meta = {
        "block": PREFIX,
        "model": model_name,
        "dims": dim,
        "w_n": w_n,
        "w_c": w_c,
        "pooling": "mean over residues, cls and eos excluded",
        "n_sequences": len(pids),
        "n_overlap": int(overlap.sum()),
        "pct_overlap": round(100.0 * float(overlap.mean()), 2),
        "length_min": int(lengths.min()),
        "length_max": int(lengths.max()),
    }
    return df, meta


def check_features(df, n_expected, n_features):
    """Sanity checks before the matrix is written.

    Embedding dimensions are unbounded reals, so the non-negative and
    sum-to-one checks that apply to composition blocks do not apply here. A
    constant column is checked instead: a dead dimension means the forward
    pass produced nothing for it.
    """
    feat = df.drop(columns="Protein_id")

    assert df.shape == (n_expected, n_features + 1), df.shape
    assert df["Protein_id"].is_unique
    assert df.notna().all().all()
    assert np.isfinite(feat.values).all()
    assert feat.values.std(axis=0).min() > 0, "a feature column is constant"

    print("checks: %d rows x %d features, all passed" % (df.shape[0], n_features))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta", required=True,
        help="FASTA file name inside the Train_Test folder, or a full path"
    )
    parser.add_argument(
        "--out", default=None,
        help="output directory (default: the repo's phase1 features/N_C_window folder)"
    )
    parser.add_argument(
        "--model", default=DEFAULT_MODEL,
        help="Hugging Face model id or local path (default: %s)" % DEFAULT_MODEL
    )
    parser.add_argument(
        "--window", default="%d,%d" % (DEFAULT_WN, DEFAULT_WC),
        help="N,C window widths in residues (default: %d,%d)"
             % (DEFAULT_WN, DEFAULT_WC)
    )
    args = parser.parse_args()

    check_project_root()

    try:
        w_n, w_c = (int(v) for v in args.window.split(","))
    except ValueError:
        raise SystemExit("--window expects two integers, e.g. 50,50")
    if w_n < 1 or w_c < 1:
        raise SystemExit("--window widths must be positive")

    fasta_path = Path(args.fasta)
    if not fasta_path.is_absolute():
        candidate = FASTA_DIR / fasta_path
        if candidate.exists():
            fasta_path = candidate
    fasta_path = fasta_path.resolve()

    if not fasta_path.exists():
        raise FileNotFoundError("no such FASTA file: %s" % fasta_path)

    out_dir = Path(args.out).resolve() if args.out else FEATURE_DIR / "N_C_window"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("block : ncw")
    print("fasta : %s" % fasta_path)

    seqs = read_fasta(fasta_path)
    check_sequences(seqs)

    features, meta = extract(seqs, args.model, w_n, w_c)
    check_features(features, len(seqs), meta["dims"])

    stem = "N_C_%d_window_esm2_%s" % (w_n, fasta_path.stem)

    out_path = out_dir / (stem + ".csv")
    features.to_csv(str(out_path), index=False, float_format="%.10f")
    print("wrote : %s" % out_path)

    meta_path = out_dir / (stem + ".meta.json")
    with open(str(meta_path), "w") as handle:
        json.dump(meta, handle, indent=2)
    print("wrote : %s" % meta_path)


if __name__ == "__main__":
    main()