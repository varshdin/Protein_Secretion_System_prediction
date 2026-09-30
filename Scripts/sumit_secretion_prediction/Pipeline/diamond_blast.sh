#!/bin/bash
# Steps:
#   1. Build complete secreted reference (Sec, Tat ,T1SS-T6SS combined)
#   2. Build DIAMOND database
#   3. Run DIAMOND BLAST (Non-secreted vs. secreted)
#   4. Flag suspicious sequences (>30% identity AND >70% coverage)
#   5. Remove flagged sequences
#   6. Renumber the cleaned Non-secreted file sequentially
#   7. Copy all other classes (unchanged) into Clean_Data
#   8. Print final summary of all 9 classes
#
# Run this from WSL/Ubuntu (DIAMOND is a Linux tool).
#
# Usage:
#     bash run_diamond_pipeline.sh
# ============================================================

set -e
pwd

cd "$(dirname "$0")/../../../"
REPO_ROOT="$(pwd)"

CD_HIT_DIR="${REPO_ROOT}/Results/sumit_secretion_prediction/phase0_dataset_prep/CD_HIT_60_per"
DIAMOND_DIR="${REPO_ROOT}/Results/sumit_secretion_prediction/phase0_dataset_prep/DIAMOND_BLAST"
CLEAN_DATA_DIR="${REPO_ROOT}/Results/sumit_secretion_prediction/phase0_dataset_prep/Clean_Data"
DIAMOND_SCRIPTS_DIR="${REPO_ROOT}/Scripts/sumit_secretion_prediction/DIAMOND_BLAST"
SECRETED_CLASSES=("T1SS" "T2SS" "T3SS" "T4SS" "T5SS" "T6SS")
REFERENCE_CLASSES=("T1SS" "T2SS" "T3SS" "T4SS" "T5SS" "T6SS" "Sec" "Tat")

mkdir -p "$DIAMOND_DIR" "$CLEAN_DATA_DIR"

echo "=========================================="
echo "DIAMOND BLAST VALIDATION PIPELINE"
echo "(Non-secreted vs. Secreted: T1SS-T6SS + Sec + Tat)"
echo "=========================================="

echo ""
echo "--- Step 1: Building complete reference (secreted + Sec + Tat) ---"
cat "${CD_HIT_DIR}/T1SS_final.fasta" "${CD_HIT_DIR}/T2SS_final.fasta" \
    "${CD_HIT_DIR}/T3SS_final.fasta" "${CD_HIT_DIR}/T4SS_final.fasta" \
    "${CD_HIT_DIR}/T5SS_final.fasta" "${CD_HIT_DIR}/T6SS_final.fasta" \
    "${CD_HIT_DIR}/Sec_final.fasta" "${CD_HIT_DIR}/Tat_final.fasta" \
    > "${DIAMOND_DIR}/all_secreted_class.fasta"
ref_count=$(grep -c "^>" "${DIAMOND_DIR}/all_secreted_class.fasta")
echo "Reference built: $ref_count total sequences (T1SS-T6SS + Sec + Tat)"

echo ""
echo "--- Step 2: Building DIAMOND database ---"
cd "$DIAMOND_DIR"
diamond makedb --in all_secreted_class.fasta -d all_secreted_class

echo ""
echo "--- Step 3: Running DIAMOND BLAST ---"
diamond blastp -d all_secreted_class \
    -q "${CD_HIT_DIR}/Non-secreted_final.fasta" \
    -o non_secreted_vs_secreted.tsv --very-sensitive

echo ""
echo "--- Step 4: Flagging suspicious sequences ---"
python3 "${DIAMOND_SCRIPTS_DIR}/flag_suspicious_data.py" \
    non_secreted_vs_secreted.tsv "${CD_HIT_DIR}/Non-secreted_final.fasta"

echo ""
echo "--- Step 5: Removing flagged sequences ---"
python3 "${DIAMOND_SCRIPTS_DIR}/remove_flagged_sequence.py" \
    "${CD_HIT_DIR}/Non-secreted_final.fasta" flagged_non_secreted_ids.txt \
    "${CLEAN_DATA_DIR}/Non-secreted_final.fasta"

echo ""
echo "--- Step 6: Renumbering cleaned Non-secreted sequences ---"
python3 "${DIAMOND_SCRIPTS_DIR}/renumber.py" \
    "${CLEAN_DATA_DIR}/Non-secreted_final.fasta" "Non-secreted"

before=$(grep -c "^>" "${CD_HIT_DIR}/Non-secreted_final.fasta")
after=$(grep -c "^>" "${CLEAN_DATA_DIR}/Non-secreted_final.fasta")
echo "Non-secreted: $before -> $after (removed $((before - after)) likely mislabeled effectors)"

echo ""
echo "--- Step 7: Copying other classes (unchanged) into Clean_Data ---"
for class in "${REFERENCE_CLASSES[@]}"; do
    cp "${CD_HIT_DIR}/${class}_final.fasta" "${CLEAN_DATA_DIR}/${class}_final.fasta"
    count=$(grep -c "^>" "${CLEAN_DATA_DIR}/${class}_final.fasta")
    echo "$class: copied ($count sequences)"
done

echo ""
echo "=========================================="
echo "DIAMOND PIPELINE COMPLETE"
echo "=========================================="
echo ""
echo "--- Final summary: all 9 classes in Clean_Data ---"
for class in "${REFERENCE_CLASSES[@]}" "Non-secreted"; do
    file="${CLEAN_DATA_DIR}/${class}_final.fasta"
    if [ -f "$file" ]; then
        count=$(grep -c "^>" "$file")
        echo "$class: $count"
    fi
done
echo ""
echo "All final files saved to: $CLEAN_DATA_DIR"