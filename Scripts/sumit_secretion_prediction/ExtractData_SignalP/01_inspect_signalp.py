from collections import Counter
from pathlib import Path

PATH = Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Data/sumit_secretion_prediction/raw/SignalP/complete_set_unpartitioned.fasta")   # adjust to your config.py path

records = []
with PATH.open() as fh:
    lines = [l.rstrip("\n") for l in fh if l.strip()]

# 3-line FASTA: header / sequence / annotation
assert len(lines) % 3 == 0, f"not a clean 3-line file: {len(lines)} lines"

for i in range(0, len(lines), 3):
    header, seq, ann = lines[i], lines[i+1], lines[i+2]
    assert header.startswith(">"), f"bad header at line {i}: {header[:50]}"
    assert len(seq) == len(ann), f"seq/annotation length mismatch: {header}"
    fields = header[1:].split("|")
    records.append({
        "accession": fields[0],
        "kingdom":   fields[1],
        "type":      fields[2],
        "partition": fields[3] if len(fields) > 3 else None,
        "seq":       seq,
    })

print(f"total records: {len(records)}\n")
print("--- first 2 headers verbatim ---")
for r in records[:2]:
    print(">" + "|".join([r["accession"], r["kingdom"], r["type"], str(r["partition"])]))

print("\n--- kingdom values ---")
for k, n in Counter(r["kingdom"] for r in records).most_common():
    print(f"{n:6d}  {k}")

print("\n--- type values ---")
for t, n in Counter(r["type"] for r in records).most_common():
    print(f"{n:6d}  {t}")

print("\n--- type x kingdom ---")
for (k, t), n in sorted(Counter((r["kingdom"], r["type"]) for r in records).items()):
    print(f"{n:6d}  {k:12s} {t}")

print("\n--- sequence lengths ---")
lens = [len(r["seq"]) for r in records]
print(f"min {min(lens)}  max {max(lens)}  unique values: {len(set(lens))}")