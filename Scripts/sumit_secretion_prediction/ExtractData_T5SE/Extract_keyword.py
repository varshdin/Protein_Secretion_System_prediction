"""
For every protein sequence, extracts the TITLE words only -- not the
amino acid sequence, not the ID, not the OS=/OX=/GN=/PE=/SV= metadata.

Example header:
    >sp|Q8VSL2|SEPA_SHIFL Serine protease SepA autotransporter OS=Shigella flexneri OX=623 GN=sepA PE=1 SV=1

    sp|Q8VSL2|SEPA_SHIFL          <- ID section, IGNORED
    Serine protease SepA autotransporter   <- TITLE, this is what we extract
    OS=Shigella flexneri OX=623 GN=sepA PE=1 SV=1   <- metadata, IGNORED

From the title, individual WORDS are extracted:
    Serine, protease, SepA, autotransporter

This is done for EVERY sequence in the file, then every UNIQUE word found
across the WHOLE dataset is listed once, with how many times it appeared --
so each word can be reviewed and classified (KEEP / DISCARD) with real
evidence, one at a time, before building the final filter.

Usage:
    python3 collect_unique_keywords.py <input_filename>

    Pass ONLY the filename -- the script looks for it inside T5SE_DIR
    (see FILE CONFIGURATION below).
"""

import sys
import re
from pathlib import Path
from collections import Counter
from Bio import SeqIO

# ============================================================
# FILE CONFIGURATION -- check this path before running
# ============================================================
THIS_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_SCRIPT_DIR.parent.parent.parent   # comp-ss/ root

T5SE_DIR = REPO_ROOT / "Data" / "sumit_secretion_prediction" / "raw" / "T5SE"
# ============================================================


def extract_title(header):
    """
    Given a full UniProt header, returns ONLY the title portion --
    everything between the "sp|ID|NAME " section and the " OS=" metadata.

    Example:
        Input:  "sp|Q8VSL2|SEPA_SHIFL Serine protease SepA autotransporter OS=Shigella flexneri OX=623 GN=sepA PE=1 SV=1"
        Output: "Serine protease SepA autotransporter"
    """
    # Step 1: cut off everything from " OS=" onward (organism + metadata)
    before_os = re.split(r"\sOS=", header)[0]

    # Step 2: remove the leading "sp|ID|NAME" section (first whitespace-
    # separated token, which contains the pipe characters)
    parts = before_os.split(" ", 1)
    if len(parts) < 2:
        return ""
    title = parts[1]
    return title.strip()


def extract_words(title):
    """
    Splits a title string into individual words.
    Example: "Serine protease SepA autotransporter"
           -> ["Serine", "protease", "SepA", "autotransporter"]
    """
    return re.findall(r"[A-Za-z][A-Za-z0-9\-]*", title)


def main(input_filename):
    input_fasta = T5SE_DIR / input_filename

    if not input_fasta.exists():
        print(f"ERROR: Input file not found at:")
        print(f"  {input_fasta}")
        print(f"\nPlace your FASTA file in:")
        print(f"  {T5SE_DIR}")
        print(f"then run again, passing just the filename.")
        sys.exit(1)

    word_counts = Counter()
    total_sequences = 0

    for record in SeqIO.parse(input_fasta, "fasta"):
        total_sequences += 1
        title = extract_title(record.description)
        if not title:
            continue
        words = extract_words(title)
        for word in words:
            word_counts[word] += 1

    output_path = T5SE_DIR / "unique_title_words.tsv"
    with open(output_path, "w") as f:
        f.write("count\tword\n")
        for word, count in word_counts.most_common():
            f.write(f"{count}\t{word}\n")

    print(f"Total protein sequences processed: {total_sequences}")
    print(f"Total UNIQUE words found in titles: {len(word_counts)}")
    print(f"\nSaved to: {output_path}")
    print(f"\nFirst 30 words (most frequent first):")
    for word, count in word_counts.most_common(30):
        print(f"  {count:>6}  {word}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 collect_unique_keywords.py <input_filename>")
        print(f"(File must already be placed inside: {T5SE_DIR})")
        sys.exit(1)
    main(sys.argv[1])