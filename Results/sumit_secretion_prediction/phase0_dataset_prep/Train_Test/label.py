"""
Reads Train.fasta / Test.fasta and writes new FASTA files where each
header is simplified to JUST the class label (e.g. ">T1SS" instead of
">T1SS-1"), dropping the sequence number.

WARNING: This means many sequences will share the EXACT SAME header
(e.g. all 148 T1SS sequences become ">T1SS"). This is valid FASTA
syntax, but means individual sequences can no longer be distinguished
or traced back to their original identity. Consider whether you
actually need this, versus reading the label directly from the
original unique IDs during feature engineering instead.

Usage:
    python3 simplify_fasta_labels.py
"""

import re
from pathlib import Path
from Bio import SeqIO

TRAIN_TEST_DIR = Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Results/sumit_secretion_prediction/phase0_dataset_prep/Train_Test")

TRAIN_FASTA = TRAIN_TEST_DIR / "Train.fasta"
TEST_FASTA = TRAIN_TEST_DIR / "Test.fasta"


def extract_label(sequence_id):
    return re.sub(r"-\d+$", "", sequence_id)


def main():
    for name, fasta_path in [("Train", TRAIN_FASTA), ("Test", TEST_FASTA)]:
        if not fasta_path.exists():
            print(f"ERROR: {name} file not found at {fasta_path}")
            continue

        output_path = TRAIN_TEST_DIR / f"{name}_label_only.fasta"
        count = 0
        with open(output_path, "w") as out:
            for record in SeqIO.parse(fasta_path, "fasta"):
                label = extract_label(record.id)
                out.write(f">{label}\n{str(record.seq)}\n")
                count += 1

        print(f"{name}: {count} sequences written with simplified headers -> {output_path}")

    print("\nWARNING: headers within each class are now identical (e.g. all")
    print("T1SS sequences are labeled '>T1SS') -- individual sequence identity")
    print("has been lost. Keep your original Train.fasta/Test.fasta as backup.")


if __name__ == "__main__":
    main()