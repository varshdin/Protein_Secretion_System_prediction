"""
Filters T5SS candidate sequences by checking each entry's TITLE words
against a KEEP list and a DISCARD list -- both built directly from a
word-frequency analysis of this project's own real, downloaded data
(collect_unique_keywords.py), reviewed and classified word by word.

======================================================================
EVIDENCE FOR EACH DECISION
======================================================================

KEEP -- confirmed genuine T5SS terminology:

  autotransporter        -> core, defining term (Desvaux et al. 2004;
                             Henderson et al. 2004; Abby et al. 2016)
  beta-barrel / barrel /
  beta                    -> confirmed translocator structure:
                             "inserts into the OM as a beta-barrel"
  outer                   -> T5SS translocator is an outer-membrane
                             protein (confirmed throughout mechanism)
  yada                    -> confirmed real, named T5cSS marker protein,
                             directly found in this project's own data
  trimeric                -> confirmed T5cSS-specific marker (field's
                             standard name for this sub-type)
  adhesin                 -> confirmed passenger domain function
                             (Henderson et al. 2004, explicit category)
  serine / protease /
  peptidase                -> confirmed passenger domain function
                             (Henderson et al. 2004: "protease, peptidase");
                             matches real named examples in this data
                             (EspP, SepA, Pet, EatA)
  esterase                -> confirmed passenger domain function
                             (Henderson et al. 2004: "lipase, esterase")
  hemagglutinin           -> confirmed real T5bSS marker, found
                             repeatedly in this project's FhaB/FhaC data
  patatin                 -> confirmed T5dSS marker (PlpD/FplA family)
  tps / tpsa / tpsb        -> confirmed T5bSS terminology (two-partner
                             secretion: TpsA=passenger, TpsB=translocator)

DISCARD -- confirmed real false-positive terms, found directly in this
project's own gene-name search results:

  tonb                    -> confirmed collision found in "FhaB" search
                             (unrelated iron-transport protein)
  sigma                   -> confirmed collision with gene "sigA"
                             (unrelated RNA polymerase sigma factor)
  sulfate                 -> confirmed collision with gene "sat"
                             (unrelated sulfate adenylyltransferase)
  adenylyltransferase     -> same collision as above
  quorum-quenching        -> confirmed collision with gene "aidA" in
                             Acinetobacter baumannii (unrelated enzyme)
  glycosylase             -> unrelated DNA-repair enzyme class
  ribosomal                -> unrelated translation machinery

AMBIGUOUS -- flagged for manual review, not auto-included or excluded,
since evidence is insufficient to decide by keyword alone:

  fragment                -> confirmed to mark INCOMPLETE sequences
                             (needs separate handling, not blind inclusion)
  uncharacterized          -> no functional description at all; relies
                             purely on the Pfam match as sole evidence
  putative                -> explicitly means "not confirmed, only
                             suspected" -- a hedge word, not a
                             keep/discard signal by itself
  porin                    -> no direct, confirmed evidence found linking
                             standalone porin entries to T5SS; typically
                             a distinct outer-membrane protein family;
                             needs individual verification before deciding

TOO GENERIC -- ignored entirely, provide zero discriminating value:

  protein, domain-containing, domain, family

======================================================================

Usage:
    python3 filter_t5se_by_protein_name.py <input_filename> <output_filename>

    Pass ONLY filenames -- both are automatically looked for / saved
    inside T5SE_DIR (see FILE CONFIGURATION below).
"""

import sys
import re
from pathlib import Path
from Bio import SeqIO

# ============================================================
# FILE CONFIGURATION -- check this path before running
# ============================================================
THIS_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_SCRIPT_DIR.parent.parent.parent   # comp-ss/ root

T5SE_DIR = REPO_ROOT / "Data" / "sumit_secretion_prediction" / "raw" / "T5SE"
# ============================================================

KEEP_WORDS = [
    "autotransporter",
    "beta-barrel",
    "barrel",
    "beta",
    "outer",
    "yada",
    "trimeric",
    "adhesin",
    "serine",
    "protease",
    "peptidase",
    "esterase",
    "hemagglutinin",
    "patatin",
    "tps",
    "tpsa",
    "tpsb",
]

DISCARD_WORDS = [
    "tonb",
    "sigma",
    "sulfate",
    "adenylyltransferase",
    "quorum-quenching",
    "glycosylase",
    "ribosomal",
]

AMBIGUOUS_WORDS = [
    "fragment",
    "uncharacterized",
    "putative",
    "porin",
]

MIN_LENGTH = 100


def extract_title(header):
    before_os = re.split(r"\sOS=", header)[0]
    parts = before_os.split(" ", 1)
    return parts[1].strip() if len(parts) > 1 else ""


def main(input_filename, output_filename):
    T5SE_DIR.mkdir(parents=True, exist_ok=True)

    input_fasta = T5SE_DIR / input_filename
    output_fasta = T5SE_DIR / output_filename

    if not input_fasta.exists():
        print(f"ERROR: Input file not found at:")
        print(f"  {input_fasta}")
        print(f"\nPlace your FASTA file in:")
        print(f"  {T5SE_DIR}")
        print(f"then run again, passing just the filename.")
        sys.exit(1)

    kept, discarded, ambiguous = [], [], []

    for record in SeqIO.parse(input_fasta, "fasta"):
        title = extract_title(record.description).lower()
        seq_length = len(record.seq)

        has_keep = any(w in title for w in KEEP_WORDS)
        has_discard = any(w in title for w in DISCARD_WORDS)
        has_ambiguous = any(w in title for w in AMBIGUOUS_WORDS)
        long_enough = seq_length >= MIN_LENGTH

        if has_discard:
            discarded.append((record, "matched confirmed false-positive word"))
        elif not long_enough:
            discarded.append((record, f"too short ({seq_length} aa)"))
        elif has_ambiguous:
            ambiguous.append((record, "matched an ambiguous word, needs manual review"))
        elif has_keep:
            kept.append(record)
        else:
            discarded.append((record, "no genuine T5SS terminology found"))

    SeqIO.write(kept, output_fasta, "fasta")

    ambiguous_path = T5SE_DIR / "T5SS_ambiguous_needs_review.fasta"
    SeqIO.write([r for r, reason in ambiguous], ambiguous_path, "fasta")

    print(f"Total processed: {len(kept) + len(discarded) + len(ambiguous)}")
    print(f"\nKEPT (confirmed genuine T5SS): {len(kept)}")
    print(f"  -> saved to: {output_fasta}")
    print(f"\nAMBIGUOUS (needs manual review): {len(ambiguous)}")
    print(f"  -> saved to: {ambiguous_path}")
    print(f"\nDISCARDED: {len(discarded)}")
    for record, reason in discarded[:20]:
        print(f"  - {record.id}: {reason}")
    if len(discarded) > 20:
        print(f"  ... and {len(discarded) - 20} more")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 filter_t5se_by_protein_name.py <input_filename> <output_filename>")
        print(f"(Input file must already be placed inside: {T5SE_DIR})")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])