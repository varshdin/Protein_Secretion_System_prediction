"""
Removes flagged sequences from a FASTA file, based on the flagged IDs
produced by flag_suspicious_negatives.py


Purpose
-------
Removes sequences flagged as suspicious (likely mislabeled secreted)
from the Non-secreted class, based on the ID list produced by
flag_suspicious_data.py.

Process
-------
1. Reads the flagged sequence IDs into a lookup set.
2. Iterates through the original Non-secreted FASTA file, writing
   every sequence NOT in the flagged set to a new, cleaned output file.
3. Cross-checks that the number of sequences actually removed matches
   the number of IDs that were flagged, warning if there's a mismatch
   (commonly caused by ID formatting inconsistencies, e.g. stray
   whitespace or carriage-return characters).
4. Writes an audit log (removed_sequences_log.txt) recording exactly
   which sequence IDs were removed, for provenance and reproducibility.

Input
-----
1. The original Non-secreted FASTA file.
2. The flagged ID list (flagged_non_secreted_ids.txt) from Step 4.

Output
------
- A cleaned Non-secreted FASTA file, with suspicious sequences excluded.
- removed_sequences_log.txt: full record of what was removed and why.


Usage:
    python3 remove_flagged_sequences.py <original_fasta> <flagged_ids_txt> <output_cleaned_fasta>

Example:
    python3 remove_flagged_sequences.py Non-effector_blast.fasta flagged_non_effector_ids.txt Non-effector_cleaned.fasta
"""

import sys
from Bio import SeqIO

def main(original_fasta, flagged_ids_file, output_fasta):
    # Load flagged IDs into a set for fast lookup
    with open(flagged_ids_file) as f:
        flagged_ids = set(line.strip() for line in f if line.strip())

    print(f"Number of flagged IDs to remove: {len(flagged_ids)}")

    total_count = 0
    kept_count = 0
    removed_count = 0
    removed_list = []

    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(original_fasta, "fasta"):
            total_count += 1
            if record.id in flagged_ids:
                removed_count += 1
                removed_list.append(record.id)
            else:
                SeqIO.write(record, out_handle, "fasta")
                kept_count += 1

    print(f"\nOriginal total sequences: {total_count}")
    print(f"Removed (flagged as suspicious): {removed_count}")
    print(f"Kept (clean): {kept_count}")
    print(f"\nCleaned file saved as: {output_fasta}")

    if removed_count != len(flagged_ids):
        print(f"\nWARNING: {len(flagged_ids)} IDs were flagged, but only {removed_count} were found and removed.")
        print("Check for ID mismatches (e.g. trailing whitespace, \\r characters) between the two files.")

    # Save a log of exactly what was removed, for your thesis records
    with open("removed_sequences_log.txt", "w") as log:
        log.write(f"Total original sequences: {total_count}\n")
        log.write(f"Total removed: {removed_count}\n")
        log.write(f"Total kept: {kept_count}\n\n")
        log.write("Removed sequence IDs:\n")
        for rid in removed_list:
            log.write(rid + "\n")
    print("Removal log saved as: removed_sequences_log.txt")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 remove_flagged_sequences.py <original_fasta> <flagged_ids_txt> <output_cleaned_fasta>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3]) 