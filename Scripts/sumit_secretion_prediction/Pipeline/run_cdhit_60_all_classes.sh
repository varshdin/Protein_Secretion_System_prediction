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
# across every class in this project (T1SE-T6SE, Non-secreted,
# Sec, Tat)     
#
# Run this from WSL/Ubuntu (CD-HIT is a Linux tool).
#
# Usage:
#     bash run_cdhit_60_all_classes.sh
# ============================================================

set -e 

# Paths — adjust ONLY here if folder structure changes
INPUT_DIR="Data/sumit_secretion_prediction/raw/Combined_raw_dataset"
OUTPUT_DIR="Results/sumit_secretion_prediction/phase0_dataset_prep/CD_HIT_60_per"
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "STAGE 1: CD-HIT (60%) -- all classes"
echo "=========================================="

#CLASSES=("T1SS" "T2SS" "T3SS" "T4SS" "T5SS" "T6SS" "Non-secreted" "Sec" "Tat")
CLASSES=("Sec" "Tat")

for class in "${CLASSES[@]}"; do
    echo "Running CD-HIT for class: $class"
    cd-hit -i "$INPUT_DIR/${class}.fasta" \
           -o "$OUTPUT_DIR/${class}_final.fasta" \
           -c 0.6 -n 4
done


echo ""
echo "=========================================="
echo "STAGE 1 COMPLETE"
echo "=========================================="
echo "All final files saved to: $OUTPUT_DIR"