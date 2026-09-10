"""
Flag a non-effector sequence as "suspicious" (possibly mislabeled) if
its match to a known effector has:
    - percent identity > 30%
    - AND alignment coverage > 70% of the query (non-effector) sequence length

Usage:
    python3 flag_suspicious_negatives.py <blast_tsv_file> <non_secreted_fasta_file>
"""

import sys
from Bio import SeqIO

def get_sequence_lengths(fasta_file):
    lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths[record.id] = len(record.seq)
    return lengths

def main(tsv_file, fasta_file):
    query_lengths = get_sequence_lengths(fasta_file)

    flagged = []
    all_rows = []

    with open(tsv_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            query_id = parts[0]
            subject_id = parts[1]
            pct_identity = float(parts[2])
            aln_length = int(parts[3])
            evalue = float(parts[10])
            bitscore = float(parts[11])

            query_len = query_lengths.get(query_id)
            if query_len is None:
                continue

            coverage = (aln_length / query_len) * 100

            row = {
                "query": query_id,
                "subject": subject_id,
                "pct_identity": pct_identity,
                "aln_length": aln_length,
                "query_len": query_len,
                "coverage": round(coverage, 1),
                "evalue": evalue,
                "bitscore": bitscore,
            }
            all_rows.append(row)

            if pct_identity > 30 and coverage > 70:
                flagged.append(row)

    print(f"Total alignment rows in file: {len(all_rows)}")
    print(f"Unique non-seceted queries with matches: {len(set(r['query'] for r in all_rows))}")
    print(f"\n=== FLAGGED (identity > 30% AND coverage > 70%) ===")
    print(f"Number of suspicious non-secreted sequences: {len(set(r['query'] for r in flagged))}\n")

    flagged_sorted = sorted(flagged, key=lambda r: r["pct_identity"], reverse=True)
    print(f"{'Query':<20}{'Subject':<15}{'%Identity':<12}{'Coverage%':<12}{'E-value':<15}")
    for r in flagged_sorted:
        print(f"{r['query']:<20}{r['subject']:<15}{r['pct_identity']:<12}{r['coverage']:<12}{r['evalue']:<15}")

    flagged_ids = sorted(set(r["query"] for r in flagged))
    with open("flagged_non_secreted_ids.txt", "w") as out:
        for qid in flagged_ids:
            out.write(qid + "\n")
    print(f"\nFlagged sequence IDs saved to: flagged_non_secreted_ids.txt")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 flag_suspicious_negatives.py <blast_tsv_file> <non_secreted_fasta_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])