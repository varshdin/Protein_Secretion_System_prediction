import csv
from collections import Counter
from pathlib import Path

BASE  = Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss")
CACHE = BASE / "Data/sumit_secretion_prediction/rawSignalP/uniprot_cache"
TSV   = Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Data/sumit_secretion_prediction/raw/SignalP/sec_tat_accessions.tsv")
OUT   = BASE / "Results/sumit_secretion_prediction/phase0_dataset_prep/SignalP_Repair"
OUT.mkdir(parents=True, exist_ok=True)

def parse_fasta(text):
    """UniProt .fasta: one header line, sequence wrapped over many lines."""
    lines = text.strip().split("\n")
    return lines[0], "".join(lines[1:])

rows = list(csv.DictReader(TSV.open(), delimiter="\t"))

kept, dropped = [], []
for r in rows:
    path = CACHE / f"{r['accession']}.fasta"
    if not path.exists():
        dropped.append((r["accession"], r["class"], r["subtype"], "not_fetched")); continue

    header, full = parse_fasta(path.read_text())
    nterm = r["nterm"]

    if "(Fragment)" in header:
        dropped.append((r["accession"], r["class"], r["subtype"], "fragment")); continue

    if full[:len(nterm)] != nterm:
        dropped.append((r["accession"], r["class"], r["subtype"], "prefix_mismatch")); continue

    if len(full) == len(nterm):
        dropped.append((r["accession"], r["class"], r["subtype"], "not_longer")); continue

    kept.append({**r, "full": full, "length": len(full), "uniprot_header": header})

print(f"input {len(rows)}   kept {len(kept)}   dropped {len(dropped)}\n")
print("--- drop reasons ---")
for reason, n in Counter(d[3] for d in dropped).most_common():
    print(f"{n:5d}  {reason}")

print("\n--- kept by class / subtype ---")
for (c, s), n in sorted(Counter((k["class"], k["subtype"]) for k in kept).items()):
    print(f"{n:5d}  {c:4s} {s}")

print("\n--- new length distribution ---")
import statistics as st
for cls in ("Sec", "Tat"):
    L = sorted(k["length"] for k in kept if k["class"] == cls)
    if L:
        print(f"{cls}: n={len(L)} min={L[0]} median={st.median(L)} max={L[-1]}")

# outputs
with (OUT / "dropped.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "class", "subtype", "reason"]); w.writerows(dropped)

counter = Counter()
with (OUT / "mapping.tsv").open("w", newline="") as mf:
    mw = csv.writer(mf, delimiter="\t")
    mw.writerow(["new_id", "accession", "class", "subtype", "length", "uniprot_header"])
    for cls in ("Sec", "Tat"):
        with (OUT / f"{cls}_full.fasta").open("w") as ff:
            for k in kept:
                if k["class"] != cls:
                    continue
                counter[cls] += 1
                new_id = f"{cls}-{counter[cls]}"
                ff.write(f">{new_id}\n{k['full']}\n")
                mw.writerow([new_id, k["accession"], k["class"], k["subtype"],
                             k["length"], k["uniprot_header"]])
print(f"\nwrote {OUT}")