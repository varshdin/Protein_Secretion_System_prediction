# Predicting bacterial secretion systems — pipeline

Nine-class classification of Gram-negative
bacterial proteins: `T1SS`–`T6SS`, `Sec`, `Tat`, `Non-secreted`.

Dataset after phase 0: **2967 train / 744 test**, homology-reduced, no overlap
between splits.

---

## Requirements

Phases 0.1–0.3 use Linux-only bioinformatics tools and run in **WSL**.
Phase 1 uses conda environments on the **Windows** side and runs in Git Bash
or Anaconda Prompt. That split is why there is no single master script.

### WSL

| tool | used by |
|---|---|
| `cd-hit` | phase 0.1 |
| `diamond` | phase 0.2 |
| `mmseqs` | phase 0.3 |

### Windows — two conda environments

They cannot be merged: iFeatureOmega requires Python 3.7, transformers
requires 3.10+.

```bash
conda create -n ifeature_env python=3.7 -y
conda activate ifeature_env
pip install iFeatureOmegaCLI
```

```bash
conda create -n esm_env python=3.11 -y
conda activate esm_env
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install transformers pandas numpy
```

Versions used for the results in this repository:

```
torch         2.15.0+cpu
transformers  5.23.0
pandas        3.1.1
```

Phase 1 verified byte-identical across torch 2.14.1 → 2.15.0 and
transformers 5.21.0 → 5.23.0, so nearby versions are safe.

The ESM-2 block downloads `facebook/esm2_t6_8M_UR50D` from Hugging Face on
first use and caches it. That run needs network access.

---

## Run order

Every script resolves paths from its own location, so the working directory
does not matter.

### Phase 0.1 — within-class deduplication, CD-HIT at 60% identity

```bash
# WSL
bash Scripts/sumit_secretion_prediction/Pipeline/run_cdhit_60_all_classes.sh
```

All nine classes by default. To redo a subset, pass them as arguments —
do not edit the list in the script:

```bash
bash Scripts/.../run_cdhit_60_all_classes.sh Sec Tat
```

Output: `Results/.../phase0_dataset_prep/CD_HIT_60_per/<class>_final.fasta`

### Phase 0.2 — remove mislabelled negatives, DIAMOND

```bash
# WSL
bash Scripts/sumit_secretion_prediction/Pipeline/diamond_blast.sh
```

Blasts `Non-secreted` against all eight secreted classes and removes anything
above **30% identity AND 70% coverage** as a likely effector.

Output: `Results/.../phase0_dataset_prep/Clean_Data/<class>_final.fasta`
(nine classes; the eight secreted ones are copied unchanged)

### MANUAL STEP — build the combined reference

> **Not yet scripted.** `MMseqs2.sh` reads
> `Results/.../MMsques2/Re-configure/combine_final.fasta`, and nothing
> creates it. It must be the concatenation of all nine
> `Clean_Data/*_final.fasta` files, rebuilt whenever phase 0.2 is re-run.
>
> **If this file is stale, MMseqs2 clusters the uncleaned data, the split is
> wrong, and every number downstream is wrong with no symptom anywhere.**
> Expected header count: 3711.

### Phase 0.3 — cluster at 30% identity and split

```bash
# WSL
bash Scripts/sumit_secretion_prediction/Pipeline/MMseqs2.sh
```

MMseqs2 at `--min-seq-id 0.3 -c 0.8`, then `train_test_split.py` splits by
cluster so no cluster spans both sides.

Output: `Results/.../phase0_dataset_prep/Train_Test/{train,test}.fasta`

### Phase 1 — feature engineering

```bash
# Git Bash or Anaconda Prompt, NOT WSL
python Scripts/sumit_secretion_prediction/Pipeline/feature_engineering.py
```

Runs three blocks over both splits, ~195 s total:

| block | features | script | environment |
|---|---|---|---|
| `aac` | 20 | `feature_engineering/iFeatureOmega/extract_feature.py` | `ifeature_env` |
| `dpc` | 400 | `feature_engineering/iFeatureOmega/extract_feature.py` | `ifeature_env` |
| `ncw` | 320 | `feature_engineering/ESM-2/extract_ncw.py` | `esm_env` |

`ncw` is the ESM-2 embedding of the first 50 and last 50 residues
concatenated, mean-pooled over residues with `<cls>`/`<eos>` excluded. The
50/50 width was selected by stratified 5-fold CV on the training set only.

Useful flags:

```bash
--dry-run                    print the commands, run nothing
--blocks ncw                 one block
--splits test                one split
--out-root C:/tmp/fe_test    write elsewhere; cannot touch Results/
```

Output: `Results/.../phase1_feature_engineering/features/{AAC,DPC,N_C_window}/`
plus a `.meta.json` per ESM-2 file recording model, window and dimensions.

The run ends with a table of every file written, and warns if the blocks
disagree on the id column name or the row count.

---

## Verifying phase 1

Row counts agreeing does not prove row *order* agrees, and a reordered file
would pair one protein's composition with another's embedding — with no
symptom. Run this after a full phase 1 run:

```python
import pandas as pd
from pathlib import Path

ROOT = Path(r"<path to this repository>")
F = ROOT / "Results/sumit_secretion_prediction/phase1_feature_engineering/features"

FILES = {"AAC": "AAC/AAC_%s.csv",
         "DPC": "DPC/DPC_%s.csv",
         "NCW": "N_C_window/N_C_50_window_esm2_%s.csv"}
EXPECTED = {"train": 2967, "test": 744}

ids = {}
for split in ["train", "test"]:
    for name, pat in FILES.items():
        d = pd.read_csv(F / (pat % split))
        ids[(name, split)] = d["Protein_id"]
    n = {k[0]: len(v) for k, v in ids.items() if k[1] == split}
    assert len(set(n.values())) == 1, n
    assert list(n.values())[0] == EXPECTED[split]
    base = ids[("AAC", split)]
    for name in FILES:
        assert ids[(name, split)].equals(base), "%s/%s order differs" % (name, split)

assert not set(ids[("AAC", "train")]) & set(ids[("AAC", "test")]), "LEAKAGE"
print("phase 1 output consistent")
```

---

## Notes

**Feature CSVs are not in git** — about 35 MB per run, rebuilt in 195 s. The
`.meta.json` sidecars are tracked; they hold the provenance.

**FASTA files are marked `-text` in `.gitattributes`.** They have CRLF line
endings and every result here was produced from exactly those bytes, so git
must not normalise them. The Python readers call `.strip()`, so no carriage
return reaches any feature.

**Window conditions other than 50/50** are regenerable by setting `WINDOW` in
the CONFIG block of `feature_engineering.py` and using `--out-root`.