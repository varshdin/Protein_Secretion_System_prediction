"""
Applies SignalP 6.0 predictions to the existing Non-effector pool to find
Sec/Tat candidates hiding there. Matching sequences are REMOVED from
Non-effector and MERGED into the existing Sec_raw.fasta / Tat_raw.fasta
files (continuing the sequential numbering, not restarting it).

Expects SignalP's real "Short" output format (confirmed from actual downloaded results):
    # ID  Prediction  OTHER  SP(Sec/SPI)  LIPO(Sec/SPII)  TAT(Tat/SPI)  TATLIPO(Tat/SPII)  PILIN(Sec/SPIII)  CS Position
    Non-effector-5  SP  0.001310  0.808926  0.188770  0.000352  0.000319  0.000296  CS pos: 23-24...

The Prediction column (index 1) already gives the category label directly.
The confidence score for that prediction is read from its matching column:
    OTHER -> column 2, SP -> column 3, LIPO -> column 4,
    TAT -> column 5, TATLIPO -> column 6, PILIN -> column 7

Usage:
    python3 find_sec_tat_in_nonEffector.py <signalp_output.tsv> <confidence_threshold>

Example:
    python3 find_sec_tat_in_nonEffector.py combined_signalp_nonEffector_results.tsv 0.7

"""

import sys
from pathlib import Path
from Bio import SeqIO

# ============================================================
# FILE CONFIGURATION -- check every path here before running
# ============================================================

THIS_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_SCRIPT_DIR.parent.parent.parent   # comp-ss/ root

# ---- INPUT 1: your current Non-effector FASTA file (sequences to screen) ----
NON_EFFECTOR_INPUT = REPO_ROOT / "Data" / "sumit_secretion_prediction" / "raw" / "Non-effector.fasta"

# ---- INPUT 2: SignalP results file -- pass ONLY the filename when running,
# the script automatically looks for it inside SIGNALP_RESULTS_DIR below
SIGNALP_RESULTS_DIR = REPO_ROOT / "Results" / "sumit_secretion_prediction" / \
    "phase0_dataset_prep" / "SignalP_non_effector"

# ---- OUTPUT 1 & existing Sec file (new hits get APPENDED here) ----
SEC_OUTPUT = REPO_ROOT / "Data" / "sumit_secretion_prediction" / "raw" / "SignalP" / "Sec_raw.fasta"

# ---- OUTPUT 2 & existing Tat file (new hits get APPENDED here) ----
TAT_OUTPUT = REPO_ROOT / "Data" / "sumit_secretion_prediction" / "raw" / "SignalP" / "Tat_raw.fasta"

# ---- OUTPUT 3: final, reduced Non-effector file (sequences NOT reclassified) ----
NON_EFFECTOR_OUTPUT_DIR = REPO_ROOT / "Results" / "sumit_secretion_prediction" / \
    "phase0_dataset_prep" / "SignalP_non_effector"
NON_EFFECTOR_OUTPUT = NON_EFFECTOR_OUTPUT_DIR / "Non-effector_final_v3.fasta"

# ============================================================
# (nothing below this line needs editing)
# ============================================================

SEC_TYPES = {"SP", "LIPO"}
TAT_TYPES = {"TAT", "TATLIPO"}

# Column index (0-based) holding the probability for each prediction type
PROB_COLUMN = {
    "OTHER": 2,
    "SP": 3,
    "LIPO": 4,
    "TAT": 5,
    "TATLIPO": 6,
    "PILIN": 7,
}


def parse_signalp_predictions(signalp_file, threshold):
    """Returns dict: {sequence_id: prediction_type} for predictions >= threshold."""
    results = {}
    with open(signalp_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            seq_id = parts[0]
            prediction = parts[1].strip().upper()

            col_index = PROB_COLUMN.get(prediction)
            if col_index is None:
                continue
            try:
                probability = float(parts[col_index])
            except (ValueError, IndexError):
                continue

            if probability >= threshold:
                results[seq_id] = prediction
    return results


def count_existing_sequences(fasta_path):
    """Returns how many sequences already exist in a FASTA file (for renumbering)."""
    if not fasta_path.exists():
        return 0
    return sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))


def main(signalp_filename, threshold):
    if not NON_EFFECTOR_INPUT.exists():
        print(f"ERROR: Non-effector file not found at: {NON_EFFECTOR_INPUT}")
        sys.exit(1)

    signalp_file = SIGNALP_RESULTS_DIR / signalp_filename
    if not signalp_file.exists():
        print(f"ERROR: SignalP results file not found at:")
        print(f"  {signalp_file}")
        print(f"\nPlace your combined SignalP results file in:")
        print(f"  {SIGNALP_RESULTS_DIR}")
        print(f"then run again, passing just the filename.")
        sys.exit(1)

    NON_EFFECTOR_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    predictions = parse_signalp_predictions(signalp_file, threshold)
    print(f"Predictions loaded (>= {threshold} confidence): {len(predictions)}")

    sec_new, tat_new, remaining_non_effector = [], [], []

    for record in SeqIO.parse(NON_EFFECTOR_INPUT, "fasta"):
        prediction = predictions.get(record.id)
        if prediction in SEC_TYPES:
            sec_new.append(record)
        elif prediction in TAT_TYPES:
            tat_new.append(record)
        else:
            remaining_non_effector.append(record)

    # Continue numbering from existing Sec/Tat counts, don't restart at 1
    sec_start = count_existing_sequences(SEC_OUTPUT) + 1
    tat_start = count_existing_sequences(TAT_OUTPUT) + 1

    with open(SEC_OUTPUT, "a") as out:
        for i, record in enumerate(sec_new, start=sec_start):
            out.write(f">Sec-{i}\n{str(record.seq)}\n")

    with open(TAT_OUTPUT, "a") as out:
        for i, record in enumerate(tat_new, start=tat_start):
            out.write(f">Tat-{i}\n{str(record.seq)}\n")

    with open(NON_EFFECTOR_OUTPUT, "w") as out:
        for i, record in enumerate(remaining_non_effector, start=1):
            out.write(f">Non-effector-{i}\n{str(record.seq)}\n")

    original_total = len(sec_new) + len(tat_new) + len(remaining_non_effector)

    print(f"\nOriginal Non-effector count (input file): {original_total}")
    print(f"Reclassified as Sec: {len(sec_new)} (appended to {SEC_OUTPUT}, IDs Sec-{sec_start} to Sec-{sec_start + len(sec_new) - 1})")
    print(f"Reclassified as Tat: {len(tat_new)} (appended to {TAT_OUTPUT}, IDs Tat-{tat_start} to Tat-{tat_start + len(tat_new) - 1})")
    print(f"Remaining Non-effector: {len(remaining_non_effector)} -> saved to {NON_EFFECTOR_OUTPUT}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 find_sec_tat_in_nonEffector.py <signalp_results_filename.txt> <confidence_threshold>")
        print(f"(File must already be placed inside: {SIGNALP_RESULTS_DIR})")
        sys.exit(1)
    main(sys.argv[1], float(sys.argv[2]))