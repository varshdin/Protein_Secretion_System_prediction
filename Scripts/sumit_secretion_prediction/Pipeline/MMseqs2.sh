#!/bin/bash
set -e

cd "$(dirname "$0")/../../../"
REPO_ROOT="$(pwd)"

MMSEQS2_DIR="${REPO_ROOT}/Results/sumit_secretion_prediction/phase0_dataset_prep/MMsques2/Re-configure"
TMP_DIR="${MMSEQS2_DIR}/tmp"
SCRIPTS_DIR="${REPO_ROOT}/Scripts/sumit_secretion_prediction/Train_Test_Split"
OUTPUT_D="${REPO_ROOT}/Results/sumit_secretion_prediction/phase0_dataset_prep/Train_Test"

mkdir -p "$MMSEQS2_DIR" "$TMP_DIR"
cd "$MMSEQS2_DIR"

echo "=========================================="
echo "MMSEQS2 CLUSTERING PIPELINE (30% identity, 80% coverage)"
echo "=========================================="

echo ""
echo "--- Step 1: Build MMseqs2 database from combine_final.fasta ---"
mmseqs createdb combine_final.fasta combined_db

echo ""
echo "--- Step 2: Cluster at 30% identity, 80% coverage ---"
mmseqs cluster combined_db combined_clu_30 "$TMP_DIR" --min-seq-id 0.3 -c 0.8

echo ""
echo "--- Step 3: Export cluster assignments to TSV ---"
mmseqs createtsv combined_db combined_db combined_clu_30 cluster_result_30.tsv

echo ""
echo "=========================================="
echo "MMSEQS2 CLUSTERING COMPLETE"
echo "=========================================="
n_clusters=$(cut -f1 cluster_result_30.tsv | sort -u | wc -l)
n_sequences=$(wc -l < cluster_result_30.tsv)
echo "Clusters: $n_clusters"
echo "Sequences: $n_sequences"
echo ""
echo "--- Step 4: Build train/test split ---"
cd "$REPO_ROOT"
python3 "${SCRIPTS_DIR}/train_test_split.py"

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE"
echo "=========================================="
echo "Final output: ${OUTPUT_D}/train.fasta, test.fasta"