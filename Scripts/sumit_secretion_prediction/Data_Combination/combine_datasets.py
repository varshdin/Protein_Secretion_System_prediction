"""
combine_datasets.py

Combines DeepSecE's dataset with PLM-Effector's dataset,
per secretion type, with continuous numbering to avoid
ID conflicts.

HOW TO USE:
1. Ensure the CONFIG section below matches your actual file paths.
2. Run:  python combine_datasets.py
3. Output: one combined FASTA file per type, in the "output/" folder
   (created next to this script, unless you change OUTPUT_DIR).

No external libraries required (pure Python, no BioPython needed).
"""

import os
import re
from secretion_prediction.config import RAW_DEEPSECE, RAW_PLM_EFFECTOR, RAW_DIR


# DeepSecE's single combined file (all types mixed together)
DEEPSECE_FILE = RAW_DEEPSECE / "Train-2918.fasta"

# Standalone Non-effector file (if you want to treat it separately from PLM-Effector)
NON_EFFECTOR_FILE = RAW_DIR / "Non-effector.fasta"

# Folder containing PLM-Effector's per-type FASTA files
PLM_EFFECTOR_DIR = RAW_PLM_EFFECTOR

# Mapping from our internal type label → exact PLM-Effector filename
# Edit this if your filenames differ (e.g. T1SE_train.fasta, T1SE.fasta, etc.)
PLM_FILES = {
    "T1SE": "T1SE_train.fasta",
    "T2SE": "T2SE_train.fasta",
    "T3SE": "T3SE_train.fasta",
    "T4SE": "T4SE_train.fasta",
    "T6SE": "T6SE_train.fasta",
    "Non-effector": None,
}

TYPE_LABELS = ["T1SE", "T2SE", "T3SE", "T4SE", "T6SE", "Non-effector"]

# Output folder (relative to this script's directory)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = RAW_DIR / "Combined_raw_dataset"

# ============================================================
# CORE FUNCTIONS
# ============================================================

def read_fasta(filepath):
    """
    Reads a FASTA file and returns a list of (header, sequence)
    tuples. Handles both single-line and multi-line (wrapped)
    sequences automatically.
    """
    records = []
    header = None
    seq_lines = []

    if not os.path.exists(filepath):
        print(f"  WARNING: File not found: {filepath}")
        return records

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                # save the previous record before starting a new one
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:]  # remove ">"
                seq_lines = []
            else:
                seq_lines.append(line)
        # save the last record in the file
        if header is not None:
            records.append((header, "".join(seq_lines)))

    return records


def extract_type_label(header, known_labels):
    """
    Given a header string, figures out which secretion type
    it belongs to, by checking which known label it STARTS with.
    Handles both underscore and hyphen separators automatically.

    Example headers this can handle:
      "T1SE-1"      -> "T1SE"
      "T1SE_20"     -> "T1SE"
      "Non-effector-1" -> "Non-effector"
    """
    # Sort labels by length, longest first, so "Non-effector"
    # is checked before any shorter accidental partial match
    for label in sorted(known_labels, key=len, reverse=True):
        if header.startswith(label):
            return label
    return None  # unrecognized label - will be flagged


def find_max_number(records, label):
    """
    Given a list of (header, sequence) records already known
    to belong to ONE type (label), finds the highest existing
    number used in their headers.

    Handles headers like "T1SE-128", "T1SE_128", etc.
    Returns 0 if no valid number is found (safe starting point).
    """
    max_num = 0
    for header, _ in records:
        # extract all digits at the very end of the header
        match = re.search(r"(\d+)$", header)
        if match:
            num = int(match.group(1))
            if num > max_num:
                max_num = num
    return max_num


def write_fasta(filepath, records, line_width=60):
    """
    Writes a list of (header, sequence) tuples to a FASTA file,
    wrapping sequences at line_width characters (matching
    DeepSecE's style).
    """
    with open(filepath, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")


# ============================================================
# MAIN LOGIC
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("STEP 1 — Reading DeepSecE's combined file")
    print("=" * 60)
    deepsece_all = read_fasta(DEEPSECE_FILE)
    print(f"  Total sequences read from DeepSecE: {len(deepsece_all)}")

    # Split DeepSecE's records by type
    deepsece_by_type = {label: [] for label in TYPE_LABELS}
    unrecognized = 0
    for header, seq in deepsece_all:
        label = extract_type_label(header, TYPE_LABELS)
        if label:
            deepsece_by_type[label].append((header, seq))
        else:
            unrecognized += 1

    if unrecognized:
        print(f"  WARNING: {unrecognized} DeepSecE sequences had "
              f"unrecognized labels and were skipped.")

    for label in TYPE_LABELS:
        print(f"  {label}: {len(deepsece_by_type[label])} sequences")

    print()
    print("=" * 60)
    print("STEP 2-5 — Processing each type: read PLM, renumber, combine")
    print("=" * 60)

    summary_rows = []

    for label in TYPE_LABELS:
        print(f"\n--- Processing {label} ---")

        deepsece_records = deepsece_by_type[label]
        deepsece_count = len(deepsece_records)

        # Find the highest existing number in DeepSecE's set for this type
        max_num = find_max_number(deepsece_records, label)
        print(f"  DeepSecE has {deepsece_count} sequences, "
              f"highest number = {max_num}")

        # Determine PLM-Effector file path for this label
        plm_filename = PLM_FILES.get(label)

        # Special handling for Non-effector:
        # - If label is "Non-effector" and NON_EFFECTOR_FILE exists, use that.
        # - Else, if a PLM filename is specified, use that.
        # - Else, skip PLM for this label.
        plm_filepath = None
        if label == "Non-effector" and os.path.exists(NON_EFFECTOR_FILE):
            plm_filepath = NON_EFFECTOR_FILE
            print(f"  Using dedicated Non-effector file: {NON_EFFECTOR_FILE}")
        elif plm_filename:
            plm_filepath = os.path.join(PLM_EFFECTOR_DIR, plm_filename)
            print(f"  Using PLM-Effector file: {plm_filepath}")
        else:
            print(f"  No PLM-Effector file configured for {label}; "
                  f"using only DeepSecE data.")

        plm_records_raw = read_fasta(plm_filepath) if plm_filepath and os.path.exists(plm_filepath) else []
        if plm_filepath:
            print(f"  PLM-Effector file '{plm_filepath}' has "
                  f"{len(plm_records_raw)} sequences")
        else:
            print(f"  No PLM-Effector file read for {label}.")

        # Renumber PLM's sequences to continue after DeepSecE's max
        plm_renumbered = []
        next_num = max_num + 1
        for _, seq in plm_records_raw:
            new_header = f"{label}-{next_num}"
            plm_renumbered.append((new_header, seq))
            next_num += 1

        # Combine DeepSecE's original + PLM's renumbered sequences
        combined = deepsece_records + plm_renumbered
        combined_count = len(combined)

        # Write the combined file for this type
        safe_label = label.replace(" ", "_").replace("-", "_")
        out_path = os.path.join(OUTPUT_DIR, f"{safe_label}_combined.fasta")
        write_fasta(out_path, combined)

        print(f"  -> Combined total: {combined_count} sequences "
              f"({deepsece_count} DeepSecE + {len(plm_renumbered)} PLM)")
        print(f"  -> Written to: {out_path}")

        summary_rows.append((label, deepsece_count, len(plm_renumbered), combined_count))

    print()
    print("=" * 60)
    print("SUMMARY TABLE (before CD-HIT deduplication)")
    print("=" * 60)
    print(f"{'Type':<15}{'DeepSecE':<12}{'PLM Added':<12}{'Combined Total':<15}")
    print("-" * 54)
    total_deepsece, total_plm, total_combined = 0, 0, 0
    for label, d, p, c in summary_rows:
        print(f"{label:<15}{d:<12}{p:<12}{c:<15}")
        total_deepsece += d
        total_plm += p
        total_combined += c
    print("-" * 54)
    print(f"{'TOTAL':<15}{total_deepsece:<12}{total_plm:<12}{total_combined:<15}")

    print()
    print("NOTE: These are RAW combined counts, BEFORE CD-HIT")
    print("deduplication. Run CD-HIT next (60% threshold) on each")
    print("file in the 'output/' folder, per type, to get your")
    print("final deduplicated dataset.")
    print()
    print("Example CD-HIT command for one type:")
    first_label = TYPE_LABELS[0]
    safe_first = first_label.replace(" ", "_").replace("-", "_")
    print(f"  cd-hit -i {OUTPUT_DIR}\\{safe_first}_combined.fasta "
          f"-o {OUTPUT_DIR}\\{safe_first}_final.fasta -c 0.6 -n 4")


if __name__ == "__main__":
    main()