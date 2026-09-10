    #!/bin/bash
# ============================================================
# STAGE 1: CD-HIT (60%) -- Within-class deduplication
# ============================================================
# Runs CD-HIT at 60% identity threshold across every class's
# combined raw FASTA file, automatically -- one consistent
# command applied uniformly, rather than typing individual
# commands per class.
#
# WHY 60%: matches the established threshold used consistently
# across every class in this project (T1SE-T6SE, Non-effector,
# Sec, Tat) -- see combine_train_test_raw.py and prior pipeline
# stages for the same convention.
#
# Run this from WSL/Ubuntu (CD-HIT is a Linux tool).
#
# Usage:
#     bash run_cdhit_60_all_classes.sh
# ============================================================

set -e  # stop immediately if any step fails

# ---- EDIT THESE PATHS to match your real folder structure ----
RAW_DIR="/mnt/c/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Data/sumit_secretion_prediction/Raw/Combined_train_test"
OUTPUT_DIR="/mnt/c/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Results/sumit_secretion_prediction/phase0_dataset_prep/CD_HIT_60%"

# ---- List every class you want processed in this run ----
CLASSES=("T1SS" "T2SS" "T3SS" "T4SS" "T5SS" "T6SS" "Non-secreted" "Sec" "Tat")

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "STAGE 1: CD-HIT (60%) -- all classes"
echo "=========================================="

for class in "${CLASSES[@]}"; do
    input_file="${RAW_DIR}/${class}.fasta"
    output_file="${OUTPUT_DIR}/${class}_final.fasta"

    if [ ! -f "$input_file" ]; then
        echo ""
        echo "SKIPPED: $class -- input file not found at $input_file"
        continue
    fi

    before_count=$(grep -c "^>" "$input_file")

    echo ""
    echo "--- $class ---"
    echo "Input: $input_file ($before_count sequences)"

    cd-hit -i "$input_file" -o "$output_file" -c 0.6 -n 4

    after_count=$(grep -c "^>" "$output_file")
    removed=$((before_count - after_count))

    echo "Result: $before_count -> $after_count (removed $removed near-duplicates)"
done

echo ""
echo "=========================================="
echo "STAGE 1 COMPLETE"
echo "=========================================="
echo "All final files saved to: $OUTPUT_DIR"