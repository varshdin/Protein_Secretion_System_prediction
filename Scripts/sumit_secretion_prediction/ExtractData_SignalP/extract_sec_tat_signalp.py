"""
Extracts Sec and Tat Gram-negative sequences from SignalP 6.0's own
published training dataset (3-line FASTA format).

File format (one entry = 3 lines):
    >Uniprot_AC|Kingdom|Type|Partition_No
    AMINOACIDSEQUENCE...
    ANNOTATIONSTRING... (region labels, e.g. all "I" for intracellular)

Filtering rule:
    Kingdom = NEGATIVE (Gram-negative bacteria only)
    Type = SP or LIPO      -> goes to Sec.fasta
    Type = Tat or TatLIPO  -> goes to Tat.fasta
    Everything else (NO_SP, PILIN, other Kingdoms) is discarded
    Partition number (0/1/2) is ignored -- all three are kept

Usage:
    python3 extract_sec_tat_signalp.py <signalp_training_filename.fasta>

    Pass ONLY the filename (not a full path) -- the script automatically
    looks for it inside INPUT_DIR (see below). Place your downloaded
    SignalP training file there first.

PORTABILITY NOTE:
    BOTH the input file location AND the output location are computed
    automatically, RELATIVE to this script's own position in the repository --
    neither uses any hardcoded personal path (like "C:/Users/sumit/...").
    This means the script works correctly for anyone who clones the repo,
    on any machine, without editing anything -- as long as the repo's
    folder structure stays the same, and the input FASTA file is placed
    in the expected input folder (see INPUT_DIR below).

"""

import sys
from pathlib import Path


# This script's own location on disk, whoever is running it
THIS_SCRIPT_DIR = Path(__file__).resolve().parent

# How many folders "up" from this script to reach the repo root (comp-ss/)
# Scripts/sumit_secretion_prediction/ExtractData_SignalP/  -> 3 levels up to repo root
REPO_ROOT = THIS_SCRIPT_DIR.parent.parent.parent

# Where the raw input FASTA file is expected to live -- place your downloaded
# SignalP training file here before running the script
INPUT_DIR = REPO_ROOT / "Data" / "sumit_secretion_prediction" / "raw" / "SignalP"

# Output directory, built relative to the repo root -- works identically
# for every collaborator, regardless of their username or drive letter
OUTPUT_DIR = REPO_ROOT / "Data" / "sumit_secretion_prediction" / "raw" / "SignalP"

SEC_TYPES = {"SP"}
TAT_TYPES = {"TAT"}


def parse_3line_fasta(filepath):
    """
    Reads a 3-line-per-entry FASTA file and yields
    (header, sequence, annotation) tuples.
    """
    with open(filepath) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    for i in range(0, len(lines), 3):
        if i + 2 >= len(lines):
            break
        header = lines[i]
        sequence = lines[i + 1]
        annotation = lines[i + 2]
        yield header, sequence, annotation


def main(input_filename):
    input_path = INPUT_DIR / input_filename

    if not input_path.exists():
        print(f"ERROR: Input file not found at expected location:")
        print(f"  {input_path}")
        print(f"\nPlace your SignalP training FASTA file in:")
        print(f"  {INPUT_DIR}")
        print(f"then run again, passing just the filename.")
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Repo root resolved as: {REPO_ROOT}")
    print(f"Reading input from: {input_path}")
    print(f"Output directory: {OUTPUT_DIR}")

    sec_records = []
    tat_records = []
    discarded_counts = {}

    total = 0
    for header, sequence, annotation in parse_3line_fasta(input_path):
        total += 1
        if not header.startswith(">"):
            continue

        parts = header[1:].split("|")
        if len(parts) < 3:
            continue

        uniprot_ac = parts[0]
        kingdom = parts[1].strip().upper()
        seq_type = parts[2].strip().upper()

        if kingdom != "NEGATIVE":
            discarded_counts[f"Kingdom={kingdom}"] = discarded_counts.get(f"Kingdom={kingdom}", 0) + 1
            continue

        if seq_type in SEC_TYPES:
            sec_records.append((uniprot_ac, seq_type, kingdom, sequence))
        elif seq_type in TAT_TYPES:
            tat_records.append((uniprot_ac, seq_type, kingdom, sequence))
        else:
            discarded_counts[f"Type={seq_type}"] = discarded_counts.get(f"Type={seq_type}", 0) + 1

    sec_path = OUTPUT_DIR / "Sec_raw.fasta"
    tat_path = OUTPUT_DIR / "Tat_raw.fasta"
    metadata_path = OUTPUT_DIR / "sec_tat_metadata_log.tsv"

    with open(sec_path, "w") as out, open(metadata_path, "w") as meta:
        meta.write("new_id\toriginal_uniprot_ac\tsubtype\tkingdom\n")
        for i, (uniprot_ac, seq_type, kingdom, sequence) in enumerate(sec_records, start=1):
            new_id = f"Sec-{i}"
            out.write(f">{new_id}\n{sequence}\n")
            meta.write(f"{new_id}\t{uniprot_ac}\t{seq_type}\t{kingdom}\n")

    with open(tat_path, "w") as out, open(metadata_path, "a") as meta:
        for i, (uniprot_ac, seq_type, kingdom, sequence) in enumerate(tat_records, start=1):
            new_id = f"Tat-{i}"
            out.write(f">{new_id}\n{sequence}\n")
            meta.write(f"{new_id}\t{uniprot_ac}\t{seq_type}\t{kingdom}\n")

    print(f"\nTotal entries processed: {total}")
    print(f"\nSec (Gram-negative, SP + LIPO): {len(sec_records)} sequences")
    print(f"  -> saved to: {sec_path}")
    print(f"\nTat (Gram-negative, Tat + TatLIPO): {len(tat_records)} sequences")
    print(f"  -> saved to: {tat_path}")
    print(f"\nMetadata log (original IDs + subtype, for traceability): {metadata_path}")
    print(f"\nDiscarded breakdown:")
    for reason, count in sorted(discarded_counts.items(), key=lambda x: -x[1]):
        print(f"  {reason}: {count}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 extract_sec_tat_signalp.py <signalp_training_filename.fasta>")
        print(f"(File must already be placed inside: {INPUT_DIR})")
        sys.exit(1)
    main(sys.argv[1])