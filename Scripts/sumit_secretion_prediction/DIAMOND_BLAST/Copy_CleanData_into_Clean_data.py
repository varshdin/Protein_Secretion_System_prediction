"""
Consolidates ALL 9 classes into the Clean_Data folder, so it becomes
the single, unified source of truth for the next pipeline stage
(train/test split) -- regardless of which classes needed DIAMOND
BLAST validation and which didn't.

LOGIC:
    - Non-secreted: already validated via DIAMOND BLAST, already sitting
      in Clean_Data (Non-secreted_cleaned.fasta) -- left untouched.
    - T1SS, T2SS, T3SS, T4SS, T5SS, T6SS, Sec, Tat: these classes rest
      on experimentally validated source databases and were NOT run
      through DIAMOND validation (per the established scope decision).
      Their CD-HIT-60%-deduplicated files are COPIED as-is from
      CD_HIT_60% into Clean_Data, so every class's FINAL file lives
      in one consistent folder going forward.

Usage:
    python3 consolidate_clean_data.py

    No arguments needed -- paths configured below.
"""

import shutil
from pathlib import Path

# ============================================================
# FILE CONFIGURATION -- edit these to match your real setup
# ============================================================
CD_HIT_DIR = Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Results/sumit_secretion_prediction/phase0_dataset_prep/CD_HIT_60%")
CLEAN_DATA_DIR = Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Results/sumit_secretion_prediction/phase0_dataset_prep/Clean_Data")

# Classes that were NOT run through DIAMOND validation -- copied as-is
# from CD_HIT_60% (their CD-HIT-deduplicated file IS their final file)
CLASSES_TO_COPY_DIRECTLY = ["T1SS", "T2SS", "T3SS", "T4SS", "T5SS", "T6SS", "Sec", "Tat"]

# Non-secreted is handled separately -- it already went through DIAMOND
# validation and should already exist in Clean_Data as Non-secreted_cleaned.fasta
# ============================================================


def main():
    CLEAN_DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Consolidating all classes into Clean_Data")
    print("=" * 60)

    for class_name in CLASSES_TO_COPY_DIRECTLY:
        source = CD_HIT_DIR / f"{class_name}_final.fasta"
        destination = CLEAN_DATA_DIR / f"{class_name}_final.fasta"

        if not source.exists():
            print(f"\nSKIPPED: {class_name} -- source not found at {source}")
            continue

        shutil.copy2(source, destination)

        # Count sequences for confirmation
        with open(destination) as f:
            count = sum(1 for line in f if line.startswith(">"))

        print(f"{class_name}: copied ({count} sequences) -> {destination}")

    # Confirm Non-secreted is present (already handled by DIAMOND validation)
    non_secreted_path = CLEAN_DATA_DIR / "Non-secreted_cleaned.fasta"
    if non_secreted_path.exists():
        with open(non_secreted_path) as f:
            count = sum(1 for line in f if line.startswith(">"))
        print(f"\nNon-secreted: already present ({count} sequences) -> {non_secreted_path}")
    else:
        print(f"\nWARNING: Non-secreted_cleaned.fasta not found in Clean_Data -- "
              f"run the DIAMOND validation + removal steps first")

    print("\n" + "=" * 60)
    print("Consolidation complete")
    print("=" * 60)
    print(f"\nAll 9 classes' FINAL files now live in: {CLEAN_DATA_DIR}")
    print("This folder is now the single source for the next pipeline stage.")


if __name__ == "__main__":
    main()