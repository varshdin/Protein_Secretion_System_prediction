"""
Renumbers sequences in a single FASTA file with clean, sequential IDs,
e.g. Non-secreted-1, Non-secreted-2, ... with no gaps -- used after
remove_flagged_sequence.py, which removes flagged sequences but leaves
gaps in the original numbering.

Usage:
    python3 renumber_single_file.py <input_fasta> <class_prefix>

Example:
    python3 renumber_single_file.py Non-secreted_cleaned.fasta Non-secreted
    (overwrites the same file with clean sequential IDs)
"""

import sys
from pathlib import Path
from Bio import SeqIO


def main(input_fasta, class_prefix):
    input_path = Path(input_fasta)

    if not input_path.exists():
        print(f"ERROR: File not found: {input_path}")
        sys.exit(1)

    records = list(SeqIO.parse(input_path, "fasta"))

    with open(input_path, "w") as out:
        for i, record in enumerate(records, start=1):
            out.write(f">{class_prefix}-{i}\n{str(record.seq)}\n")

    print(f"Renumbered {len(records)} sequences: {class_prefix}-1 to {class_prefix}-{len(records)}")
    print(f"File updated in place: {input_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 renumber_single_file.py <input_fasta> <class_prefix>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])