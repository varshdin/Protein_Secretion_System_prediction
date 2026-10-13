#!/usr/bin/env python3
"""
Cluster-based train/test split for the 9-class secretion system dataset.

WHAT THIS SCRIPT DOES
---------------------
Takes the MMseqs2 clustering result and the combined FASTA, and produces a
train/test split in which:

  1. every homology cluster lies WHOLLY in train or WHOLLY in test
     -> no test sequence has a homologue above 30% identity in the training set

  2. every class receives a test set of approximately 20% of its sequences
     -> the split stays balanced despite clusters that span several classes

WHY CLUSTER-BASED SPLITTING
---------------------------
A random split places homologous sequences on both sides. A model can then score
well on the test set by recognising sequences similar to ones it memorised,
rather than by learning what determines the label. Clustering at 30% identity
(the "twilight zone", Rost 1999) groups everything that may still be homologous,
and each group is then assigned as one unit.

HOW THE BALANCE IS ACHIEVED
---------------------------
A cluster may contain sequences of several classes, so assigning it affects
several class budgets at once. The assignment is therefore:

    for every class c:   target[c] = round(TEST_FRACTION * total[c])

    process clusters in ascending order of size
        need[c] = sequences of class c in this cluster
        if  test_count[c] + need[c] <= target[c]  for EVERY class c present:
                -> test, and charge every one of those classes
        else:   -> train

Ascending size matters: singleton clusters constrain exactly one class, so they
fill the quotas with single-sequence precision. Larger and mixed clusters are
considered last, find the quotas full, and go to train. The test set therefore
consists of mutually independent sequences. Train is simply the complement.

RUNTIME CHECKS
--------------
The script aborts before writing anything if:
  - a protein received no split
  - a protein is missing from the FASTA
  - a cluster ended up spanning both splits
  - a class exceeded its test quota
  - a protein appears in both splits

INPUT
-----
  cluster_result_30.tsv   MMseqs2 createtsv output: representative <TAB> member
  combine_final.fasta     the same sequences

OUTPUT (into TRAIN_TEST_DIR)
----------------------------
  seq_cluster_class_final.csv    every protein with its cluster and class
  Final_Split.csv               the same, plus split and sequence
  cluster_split_assignment.csv  per-cluster audit trail
  train.csv / test.csv
  train.fasta / test.fasta      model input

PATHS
-----
Paths come from secretion_prediction.config, which locates the repository root
by walking up to the folder containing .git. The script therefore runs on any
machine that has the repository cloned, with no edits. Each path can still be
overridden on the command line, which is how the shell pipeline calls it.

USAGE
-----
  python cluster_train_test_split.py                       
  python cluster_train_test_split.py --tsv X --fasta Y --out Z
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

TEST_FRACTION = 0.20


# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #
def resolve_paths(args):
    """
    Decide where to read and write.

    Defaults come from secretion_prediction.config so the script is portable.
    Command-line arguments override them, which lets the shell pipeline point at
    a different clustering run without editing the file.
    """
    if args.tsv and args.fasta and args.out:
        return Path(args.tsv), Path(args.fasta), Path(args.out)

    try:
        from secretion_prediction.config import (
            MMSEQS2_DIR, CLEAN_DATA_DIR, TRAIN_TEST_DIR,
        )
    except ImportError:
        sys.exit(
            "ERROR: cannot import secretion_prediction.config.\n"
            "       Either install the package (pip install -e .) or pass\n"
            "       --tsv, --fasta and --out explicitly."
        )

    tsv = Path(args.tsv) if args.tsv else MMSEQS2_DIR / "30_80_per" / "cluster_result_30.tsv"
    fasta = Path(args.fasta) if args.fasta else CLEAN_DATA_DIR / "combine_final.fasta"
    out = Path(args.out) if args.out else TRAIN_TEST_DIR
    return tsv, fasta, out


def load_fasta(path):
    """
    Read a FASTA file into {protein_id: sequence}.

    Sequence lines are joined, so a file wrapped at a fixed width is read
    correctly. Taking only the line after the header would silently truncate
    every sequence.
    """
    seqs, sid, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if sid:
                seqs[sid] = ''.join(buf)
            sid, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if sid:
        seqs[sid] = ''.join(buf)
    return seqs


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--tsv", help="MMseqs2 cluster TSV (default: from config)")
    ap.add_argument("--fasta", help="combined FASTA (default: from config)")
    ap.add_argument("--out", help="output directory (default: from config)")
    ap.add_argument("--test-fraction", type=float, default=TEST_FRACTION)
    args = ap.parse_args()

    TSV, FASTA, OUT = resolve_paths(args)
    test_fraction = args.test_fraction
    OUT.mkdir(parents=True, exist_ok=True)

    print("TSV  :", TSV)
    print("FASTA:", FASTA)
    print("OUT  :", OUT)
    if not TSV.exists():
        sys.exit(f"ERROR: not found: {TSV}")
    if not FASTA.exists():
        sys.exit(f"ERROR: not found: {FASTA}")

    # ----------------------------------------------------------------- STEP 1
    # Read the MMseqs2 output. Two columns, one row per sequence. The class
    # label is the id prefix: "Non-secreted-885" -> "Non-secreted", so rsplit
    # on the LAST hyphen.
    print("\n=== STEP 1  read cluster TSV ===")
    tsv = pd.read_csv(TSV, sep='\t', header=None, names=['representative', 'member'])
    tsv['label'] = tsv['member'].str.rsplit('-', n=1).str[0]

    print("rows (sequences):", len(tsv))
    print("clusters:", tsv['representative'].nunique())

    # ----------------------------------------------------------------- STEP 2
    # MMseqs2 names each cluster after its representative sequence. Rename to
    # Sec_1, Sec_2, ... for readability. This is a NAME only: it does not change
    # any protein's class, and for a mixed cluster it reflects only the first
    # member seen.
    print("\n=== STEP 2  readable cluster ids ===")
    cluster_id_map, class_counters = {}, {}
    for rep in tsv['representative'].unique():
        cls = tsv.loc[tsv.representative == rep, 'label'].iloc[0]
        class_counters[cls] = class_counters.get(cls, 0) + 1
        cluster_id_map[rep] = f"{cls}_{class_counters[cls]}"

    tsv['cluster_id'] = tsv['representative'].map(cluster_id_map)

    seq_df = (tsv[['label', 'cluster_id', 'member']]
                .rename(columns={'member': 'Protein_id'})
                .sort_values(['label', 'cluster_id'])
                .reset_index(drop=True))
    seq_df.to_csv(OUT / 'seq_cluster_class_final.csv', index=False, sep=';')
    print("per-protein rows:", len(seq_df))

    # ----------------------------------------------------------------- STEP 3
    # Fix the test quota for every class BEFORE any assignment, so no decision
    # can be influenced by what has already been assigned.
    print("\n=== STEP 3  test quota per class ===")
    total = Counter(seq_df['label'])
    target = {c: round(n * test_fraction) for c, n in total.items()}

    quota = pd.DataFrame({'total': pd.Series(total), 'test_quota': pd.Series(target)})
    quota['train_expected'] = quota['total'] - quota['test_quota']
    print(quota.sort_values('total', ascending=False).to_string())

    # ----------------------------------------------------------------- STEP 4
    # For each cluster, count how many sequences of EACH class it holds. A
    # mixed cluster has more than one entry, and every one of them must be
    # checked against its own quota later.
    print("\n=== STEP 4  cluster contents ===")
    grouped = seq_df.groupby('cluster_id')

    clusters = {}
    for cluster_id, group in grouped:
        clusters[cluster_id] = Counter(group['label'])

    sizes = {}
    for cluster_id, label_counts in clusters.items():
        sizes[cluster_id] = sum(label_counts.values())

    assert sum(sizes.values()) == len(seq_df), "sequence count changed"

    mixed = [c for c in clusters if len(clusters[c]) > 1]
    print("clusters:      ", len(clusters))
    print("singletons:    ", sum(1 for c in sizes if sizes[c] == 1))
    print("multi-member:  ", sum(1 for c in sizes if sizes[c] > 1))
    print("mixed-class:   ", len(mixed))
    for c in mixed[:5]:
        print(f"   {c:16s} size={sizes[c]:2d}  {dict(clusters[c])}")

    # ----------------------------------------------------------------- STEP 5
    # Assign each cluster. Smallest first, ties broken by cluster id, so the
    # result is identical on every run without needing a random seed.
    print("\n=== STEP 5  assign clusters ===")

    def get_cluster_size(cid):
        return (sizes[cid], cid)

    order = sorted(clusters, key=get_cluster_size)

    test_count = Counter()
    assign = {}

    for cid in order:
        need = clusters[cid]

        # A cluster may enter the test set only if EVERY class inside it still
        # has room. One full class is enough to send the whole cluster to train.
        can_fit_in_test = True
        for label, count in need.items():
            if test_count[label] + count > target[label]:
                can_fit_in_test = False
                break

        if can_fit_in_test:
            assign[cid] = 'test'
            test_count.update(need)
        else:
            assign[cid] = 'train'

    seq_df['final_split'] = seq_df['cluster_id'].map(assign)

    test_clusters_count = sum(1 for v in assign.values() if v == 'test')
    train_clusters_count = sum(1 for v in assign.values() if v == 'train')
    print("test clusters: ", test_clusters_count)
    print("train clusters:", train_clusters_count)

    # ----------------------------------------------------------------- STEP 6
    print("\n=== STEP 6  attach sequences ===")
    seq_df['Sequence'] = seq_df['Protein_id'].map(load_fasta(FASTA))
    print("sequences attached:", seq_df['Sequence'].notna().sum(), "/", len(seq_df))

    # ----------------------------------------------------------------- STEP 7
    # Each of these would otherwise fail silently and produce a plausible but
    # wrong dataset.
    print("\n=== STEP 7  checks ===")
    assert seq_df['final_split'].notna().all(), "a protein got no split"
    assert seq_df['Sequence'].notna().all(), "a protein is missing from the FASTA"
    assert seq_df.groupby('cluster_id')['final_split'].nunique().max() == 1, \
        "a cluster was split across train and test"

    for c, n in Counter(seq_df.loc[seq_df.final_split == 'test', 'label']).items():
        assert n <= target[c], f"{c} over quota: {n} > {target[c]}"

    train_ids = set(seq_df.loc[seq_df.final_split == 'train', 'Protein_id'])
    test_ids = set(seq_df.loc[seq_df.final_split == 'test', 'Protein_id'])
    assert not (train_ids & test_ids), "a protein is in both splits"

    seq_df.to_csv(OUT / 'Final_Split.csv', index=False, sep=';')
    print("all checks passed")

    # ----------------------------------------------------------------- STEP 8
    print("\n=== STEP 8  write outputs ===")
    out = (seq_df[['Protein_id', 'cluster_id', 'label', 'final_split', 'Sequence']]
             .rename(columns={'cluster_id': 'Cluster_id'}))

    train = out[out.final_split == 'train'].reset_index(drop=True)
    test = out[out.final_split == 'test'].reset_index(drop=True)

    train.to_csv(OUT / 'train.csv', index=False, sep=';')
    test.to_csv(OUT / 'test.csv', index=False, sep=';')

    pd.DataFrame({'cluster_id':   list(assign),
                  'cluster_size': [sizes[c] for c in assign],
                  'n_classes':    [len(clusters[c]) for c in assign],
                  'split':        [assign[c] for c in assign]}) \
      .to_csv(OUT / 'cluster_split_assignment.csv', index=False, sep=';')

    for name, d in (('train.fasta', train), ('test.fasta', test)):
        with open(OUT / name, 'w') as f:
            for _, r in d.iterrows():
                f.write(f">{r['Protein_id']}\n{r['Sequence']}\n")
        print(f"wrote {name}: {len(d)} sequences")

    # ----------------------------------------------------------------- STEP 9
    print("\n=== STEP 9  summary ===")
    summary = seq_df.groupby(['label', 'final_split']).size().unstack(fill_value=0)
    for c in ('train', 'test'):
        if c not in summary:
            summary[c] = 0
    summary['total'] = summary['train'] + summary['test']
    summary['quota'] = summary.index.map(target)
    summary['test%'] = (summary['test'] / summary['total'] * 100).round(1)
    summary = summary[['total', 'train', 'test', 'quota', 'test%']] \
                .sort_values('total', ascending=False)
    print(summary.to_string())
    summary.to_csv(OUT / 'split_summary.csv', sep=';')

    print(f"\ntrain {len(train)}   test {len(test)}   "
          f"({len(test) / len(seq_df) * 100:.1f}% test)")
    print(f"test: {len(test)} sequences in {test['Cluster_id'].nunique()} clusters")


if __name__ == "__main__":
    main()