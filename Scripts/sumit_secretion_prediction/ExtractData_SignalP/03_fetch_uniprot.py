import csv, time
from pathlib import Path
import requests

CACHE = Path("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Data/sumit_secretion_prediction/rawSignalP/uniprot_cache"); CACHE.mkdir(parents=True, exist_ok=True)
URL   = "https://rest.uniprot.org/uniprotkb/{}.fasta"

rows = list(csv.DictReader(open("C:/Users/sumit/OneDrive/Documenten/Thesis/comp-ss/Data/sumit_secretion_prediction/raw/SignalP/sec_tat_accessions.tsv"), delimiter="\t"))
session = requests.Session()
session.headers["User-Agent"] = "thesis-dataset-build (sumit, Frankfurt UAS)"

failed = []
for i, r in enumerate(rows, 1):
    acc  = r["accession"]
    dest = CACHE / f"{acc}.fasta"
    if dest.exists():
        continue                       # restartable: already fetched
    try:
        resp = session.get(URL.format(acc), timeout=30)
    except Exception as e:
        failed.append((acc, repr(e))); continue

    if resp.status_code == 200 and resp.text.startswith(">"):
        dest.write_text(resp.text)
    else:
        failed.append((acc, resp.status_code))
    time.sleep(0.2)

    if i % 200 == 0:
        print(f"{i}/{len(rows)}  cached={len(list(CACHE.glob('*.fasta')))}  failed={len(failed)}")

print(f"\ndone. cached={len(list(CACHE.glob('*.fasta')))}  failed={len(failed)}")
with open("fetch_failures.tsv", "w") as fh:
    for acc, why in failed:
        fh.write(f"{acc}\t{why}\n")