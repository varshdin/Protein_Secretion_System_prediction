"""
build_train_test_split.py

Purpose
-------
Takes the MMseqs2 clustering output (cluster_result_30.tsv) and the
deduplicated combine_final.fasta, and produces the final leakage-checked
train/test split for all 9 secretion classes.

This script consolidates the transformation logic originally developed
and verified in Formula_perform.ipynb (cells: afb0c1f3, e106a120,
2212cf6d, 62562934, 87713aec, 96a65b62) into one deterministic,
repeatable script with dynamic paths, replacing manual notebook
execution.

Process
-------
1. Read cluster_result_30.tsv (representative -> member pairs).
2. Assign each cluster a readable, per-class sequential ID.
3. Compute cluster sizes and cumulative sequence coverage per class,
   sorting clusters largest-first within each class.
4. Assign train/test at the CLUSTER level: clusters are added to train
   until a class's cumulative coverage exceeds 80%, then remaining
   clusters for that class go to test. This guarantees that no cluster
   (i.e. no group of near-duplicate/homologous sequences) is split
   across train and test, preventing data leakage.
5. Map the cluster-level split decision back onto every individual
   sequence, attach real sequence data from combine_final.fasta, and
   write train.csv / test.csv / train.fasta / test.fasta.

Input
-----
- Results/.../MMsques2/cluster_result_30.tsv
- Results/.../MMsques2/combine_final.fasta

Output
------
- Results/.../MMsques2/seq_cluster_class_final.csv   (per-sequence, with cluster ID)
- Results/.../MMsques2/cluster_summary_train_test.csv (per-cluster, with split decision)
- Results/.../MMsques2/train.csv, test.csv
- Results/.../MMsques2/train.fasta, test.fasta

Usage
-----
    python3 build_train_test_split.py

Called automatically as Step C of Pipeline/run_mmseqs2_full_pipeline.sh;
can also be run standalone if cluster_result_30.tsv and
combine_final.fasta already exist.
"""

import pandas as pd
from secretion_prediction.config import MMSEQS2_DIR

TSV_PATH = MMSEQS2_DIR / "cluster_result_30.tsv"
FASTA_PATH = MMSEQS2_DIR / "combine_final.fasta"
TRAIN_TEST_CUTOFF = 80  # percent


def load_fasta(path):
    seqs = {}
    current_id, current_seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = "".join(current_seq)
    return seqs


def main():
    print("=" * 60)
    print("STEP 1 — Load cluster_result_30.tsv, assign class labels")
    print("=" * 60)
    df = pd.read_csv(TSV_PATH, sep="\t", header=None, names=["representative", "member"])
    df["label"] = df["member"].apply(lambda x: x.rsplit("-", 1)[0])
    print(f"Loaded {len(df)} sequence rows, {df['representative'].nunique()} clusters")

    print("\n" + "=" * 60)
    print("STEP 2 — Assign readable per-class cluster IDs")
    print("=" * 60)
    cluster_id_map = {}
    class_counters = {}
    for rep in df["representative"].unique():
        cls = df[df["representative"] == rep]["label"].iloc[0]
        class_counters.setdefault(cls, 0)
        class_counters[cls] += 1
        cluster_id_map[rep] = f"{cls}_{class_counters[cls]}"
    df["cluster_id"] = df["representative"].map(cluster_id_map)

    seq_cluster_df = df[["label", "cluster_id", "member"]].rename(columns={"member": "Protein_id"})
    seq_cluster_df = seq_cluster_df.sort_values(["label", "cluster_id"]).reset_index(drop=True)
    seq_cluster_path = MMSEQS2_DIR / "seq_cluster_class_final.csv"
    seq_cluster_df.to_csv(seq_cluster_path, index=False, sep=";")
    print(f"Saved: {seq_cluster_path} ({len(seq_cluster_df)} rows)")

    print("\n" + "=" * 60)
    print("STEP 3 — Compute cluster sizes and cumulative coverage")
    print("=" * 60)
    cluster_sizes = seq_cluster_df.groupby(["label", "cluster_id"]).size().reset_index(name="cluster_size")
    cluster_sizes = cluster_sizes.sort_values(["label", "cluster_size"], ascending=[True, False]).reset_index(drop=True)
    cluster_sizes["cumulative_protein"] = cluster_sizes.groupby("label")["cluster_size"].cumsum()
    class_totals = cluster_sizes.groupby("label")["cluster_size"].sum()
    cluster_sizes["class_total_protein"] = cluster_sizes["label"].map(class_totals)
    cluster_sizes["percentage_seq_cover"] = (
        cluster_sizes["cumulative_protein"] / cluster_sizes["class_total_protein"] * 100
    ).round(2)
    print(f"Computed sizes for {len(cluster_sizes)} clusters across {cluster_sizes['label'].nunique()} classes")

    print("\n" + "=" * 60)
    print(f"STEP 4 — Assign train/test at cluster level (cutoff: {TRAIN_TEST_CUTOFF}%)")
    print("=" * 60)
    cluster_sizes["Train/Test_Split"] = cluster_sizes["percentage_seq_cover"].apply(
        lambda x: "train" if x <= TRAIN_TEST_CUTOFF else "test"
    )
    cluster_summary_path = MMSEQS2_DIR / "cluster_summary_train_test.csv"
    cluster_sizes.to_csv(cluster_summary_path, index=False, sep=";")
    print(f"Saved: {cluster_summary_path}")

    print("\n" + "=" * 60)
    print("STEP 5 — Map split decision to sequences, attach real sequence data")
    print("=" * 60)
    split_map = cluster_sizes.set_index("cluster_id")["Train/Test_Split"].to_dict()
    seq_cluster_df["final_split"] = seq_cluster_df["cluster_id"].map(split_map)

    missing = seq_cluster_df["final_split"].isna().sum()
    print(f"Sequences with no split assigned (should be 0): {missing}")
    if missing > 0:
        raise RuntimeError(f"{missing} sequences failed to get a split assignment — check cluster_id consistency")

    seqs = load_fasta(FASTA_PATH)
    seq_cluster_df["Sequence"] = seq_cluster_df["Protein_id"].map(seqs)

    missing_seqs = seq_cluster_df["Sequence"].isna().sum()
    print(f"Sequences with no FASTA match found (should be 0): {missing_seqs}")
    if missing_seqs > 0:
        raise RuntimeError(f"{missing_seqs} Protein_ids not found in {FASTA_PATH.name}")

    final_df = seq_cluster_df[["Protein_id", "cluster_id", "label", "final_split", "Sequence"]]
    final_df = final_df.rename(columns={"cluster_id": "Cluster_id"})

    print("\n" + "=" * 60)
    print("STEP 6 — Split and save train/test CSV + FASTA")
    print("=" * 60)
    train_df = final_df[final_df["final_split"] == "train"].reset_index(drop=True)
    test_df = final_df[final_df["final_split"] == "test"].reset_index(drop=True)

    train_csv_path = MMSEQS2_DIR / "train.csv"
    test_csv_path = MMSEQS2_DIR / "test.csv"
    train_df.to_csv(train_csv_path, index=False, sep=";")
    test_df.to_csv(test_csv_path, index=False, sep=";")

    train_fasta_path = MMSEQS2_DIR / "train.fasta"
    test_fasta_path = MMSEQS2_DIR / "test.fasta"
    with open(train_fasta_path, "w") as f:
        for _, row in train_df.iterrows():
            f.write(f">{row['Protein_id']}\n{row['Sequence']}\n")
    with open(test_fasta_path, "w") as f:
        for _, row in test_df.iterrows():
            f.write(f">{row['Protein_id']}\n{row['Sequence']}\n")

    print(f"Train: {len(train_df)} sequences -> {train_csv_path.name}, {train_fasta_path.name}")
    print(f"Test:  {len(test_df)} sequences -> {test_csv_path.name}, {test_fasta_path.name}")

    print("\n" + "=" * 60)
    print("PER-CLASS SUMMARY")
    print("=" * 60)
    summary = final_df.groupby(["label", "final_split"]).size().unstack(fill_value=0)
    summary["total"] = summary.sum(axis=1)
    summary["train_pct"] = (summary.get("train", 0) / summary["total"] * 100).round(1)
    print(summary)


if __name__ == "__main__":
    main()