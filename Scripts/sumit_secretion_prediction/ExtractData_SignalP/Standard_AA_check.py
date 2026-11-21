from collections import Counter
from pathlib import Path

OUT =Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Data/sumit_secretion_prediction/raw/SignalP")
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")

def read_fasta(p):
    txt = p.read_text().strip().split("\n")
    return dict(zip([l[1:] for l in txt[0::2]], txt[1::2]))

flagged = []
for cls in ("Sec", "Tat"):
    seqs = read_fasta(OUT / f"{cls}_full.fasta")
    print(f"\n=== {cls}: {len(seqs)} sequences ===")

    alphabet = Counter(c for s in seqs.values() for c in s)
    odd = {c: n for c, n in alphabet.items() if c not in STANDARD}
    print("non-standard characters found:", odd if odd else "none")

    for sid, s in seqs.items():
        bad = sorted(set(s) - STANDARD)
        if bad:
            flagged.append((sid, cls, "".join(bad),
                            sum(s.count(b) for b in bad), len(s)))

print(f"\n--- sequences containing non-standard residues: {len(flagged)} ---")
for sid, cls, chars, n, L in flagged:
    print(f"{sid:12s} {cls:4s} chars={chars} count={n} length={L}")