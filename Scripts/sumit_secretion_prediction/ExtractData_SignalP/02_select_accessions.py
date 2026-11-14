from collections import Counter
from pathlib import Path
import csv

SRC = Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Data/sumit_secretion_prediction/raw/SignalP/complete_set_unpartitioned.fasta")   

OUT = Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Data/sumit_secretion_prediction/raw/SignalP/sec_tat_accessions.tsv")

KINGDOM  = "NEGATIVE"
CLASS_OF = {"SP": "Sec", "TAT": "Tat"}

lines = [l.rstrip("\n") for l in SRC.open() if l.strip()]
assert len(lines) % 3 == 0, f"not 3-line FASTA: {len(lines)} lines"

rows = []
for i in range(0, len(lines), 3):
    header, seq, ann = lines[i], lines[i+1], lines[i+2]
    assert len(seq) == len(ann), header
    acc, kingdom, typ = header[1:].split("|")[:3]
    if kingdom != KINGDOM or typ not in CLASS_OF:
        continue
    rows.append({"accession": acc, "class": CLASS_OF[typ],
                 "subtype": typ, "nterm": seq})

print(f"selected: {len(rows)}   (expected 2469)")
for (c, s), n in sorted(Counter((r["class"], r["subtype"]) for r in rows).items()):
    print(f"  {n:5d}  {c:4s} {s}")

dups = [a for a, n in Counter(r["accession"] for r in rows).items() if n > 1]
print(f"duplicate accessions: {len(dups)}  {dups[:5]}")

with OUT.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["accession", "class", "subtype", "nterm"],
                       delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"wrote {OUT}")