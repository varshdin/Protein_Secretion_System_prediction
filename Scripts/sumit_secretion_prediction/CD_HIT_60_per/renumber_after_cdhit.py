#!/usr/bin/env python3
"""
Renumbers every sequence in each class's post-CD-HIT FASTA file with
clean, sequential IDs -- e.g. T1SS-1, T1SS-2, T1SS-3, ... with no gaps,
even though CD-HIT removed some sequences (leaving gaps in the original
numbering).

Usage:
    python3 renumber_after_cdhit.py
"""
 
from pathlib import Path
from Bio import SeqIO
from secretion_prediction.config import CDHIT_DIR

CLASSES = ["T1SS", "T2SS", "T3SS", "T4SS", "T5SS", "T6SS", "Non-secreted", "Sec", "Tat"]
# ============================================================


def main():
    print("=" * 50)
    print("Renumbering sequences after CD-HIT")
    print("=" * 50)

    for class_name in CLASSES:
        input_path = CDHIT_DIR / f"{class_name}_final.fasta"

        if not input_path.exists():
            print(f"\nSKIPPED: {class_name} -- file not found at {input_path}")
            continue

        records = list(SeqIO.parse(input_path, "fasta"))

        # Overwrite the same file with renumbered sequences
        with open(input_path, "w") as out:
            for i, record in enumerate(records, start=1):
                out.write(f">{class_name}-{i}\n{str(record.seq)}\n")

        print(f"\n{class_name}: {len(records)} sequences renumbered "
              f"({class_name}-1 to {class_name}-{len(records)})")

    print("\n" + "=" * 50)
    print("Renumbering complete for all classes")
    print("=" * 50)


if __name__ == "__main__":
    main()