# Secretion Protein & Secretion Channel Prediction

**Thesis project:**  Machine and Deep Learning frameworks for predicting 
bacterial secretion system effectors and channels from protein sequence data.

Part of the **HYDROCOW** project.

---

## 1. Author & Supervision

| Role | Name |
|---|---|
| Student | Sumit — M.Sc. High Integrity Systems, Frankfurt UAS |
| First Supervisor | Prof. Dr. Axel Blokesch |
| External Supervisor | Dr. Karan Kumar (Postdoc, iAMB, RWTH Aachen) |

---

## 2. Project Goal

Build and evaluate machine learning and deep learning models to classify 
bacterial protein sequences into secretion effector types (T1SE–T6SE) 
and secretion system/channel types (T1SS–T6SS), using classical ML 
(Random Forest, XGBoost) and deep learning (ESM2 protein language model 
embeddings).

---

## 3. Prerequisites (before pip install)

This pipeline requires CD-HIT, DIAMOND, MMseqs2, iFeatureomega ==> inux-native 
bioinformatics tools, not installable via pip.

- **Linux/Mac:** install directly via conda:
```bash
  conda install -c bioconda cd-hit diamond mmseqs2 iFeatureomega
```
- **Windows:** install WSL first (one-time, in an admin PowerShell):
```powershell
  wsl --install
```
  Then clone and set up the repo from inside the WSL Ubuntu terminal 
  (not Windows Git Bash), and install the tools there the same way.

**Repository convention:** all git operations, editing, and README work 
happen in Windows. All pipeline scripts (CD-HIT, DIAMOND, MMseqs2 steps) 
must be run from WSL Ubuntu, since these tools require Linux.

---

## 4. Folder Structure
```
Data/sumit_secretion_prediction/
├── raw/
│   ├── deepsecE/                    # DeepSecE train/test FASTA
│   ├── PLM-Effector/                # T1SE–T6SE train/test FASTA
│   ├── SignalP/                     # Sec/Tat raw + training FASTA
│   ├── T5SE/                        # T5SE FASTA
│   └── Combined_raw_dataset/        # all sources merged, per class

Results/sumit_secretion_prediction/
├── phase0_dataset_prep/
│   ├── CD_HIT_60_per/               # CD-HIT dedup output (60% identity)
│   ├── DIAMOND_BLAST/               # DIAMOND-cleaned sequences
│   ├── Clean_Data/                  # final cleaned dataset (CD-HIT + DIAMOND)
│   ├── MMsques2/                    # separate clustering comparison (see Section 9)
│   └── Train_Test/                  # PRIMARY final train/test split
├── phase1_feature_engineering/
│   ├── features/
│   └── processed/
├── phase2_classical_ml/           
├── phase3_esm_embeddings/           
└── final_evaluation/                

Scripts/sumit_secretion_prediction/
├── Data_Combination/
├── ExtractData_SignalP/
├── ExtractData_T5SE/
├── CD_HIT_60_per/
├── DIAMOND_BLAST/
├── feature_engineering/
│   ├── iFeatureOmega/
│   ├── POSSUM/
│   └── combine/
├── Train_Test_Split/
└── Pipeline/

src/secretion_prediction/            # reusable Python package (config.py)
Models/sumit_secretion_prediction/
Figures/sumit_secretion_prediction/
tests/test_secretion_prediction/
docs/sumit_secretion_prediction/
```

---

## 5. Pipeline — Phase 0 (Dataset Preparation)

Run in this order:

1. **`Data_Combination/combine_datasets.py`**
   Merges DeepSecE, PLM-Effector, SignalP, and T5SE raw sources into 
   `Data/.../raw/Combined_raw_dataset/`, organized per secretion class.

2. **`ExtractData_SignalP/`**
   - `extract_sec_tat_signalp.py` — extracts Sec/Tat signal peptide sequences
   - `find_sec_tat_from_non_effector.py` — identifies Sec/Tat sequences 
     among non-effector data

3. **`ExtractData_T5SE/`**
   - `Extract_keyword.py`, `filter_t5se_by_protein_name.py` — extract 
     and filter T5SE sequences by annotation keywords

4. **`CD_HIT_60_per/`** (via `Pipeline/run_cdhit_60_all_classes.sh`)
   Deduplicates sequences within each class at 60% identity threshold. 
   Output → `Results/.../CD_HIT_60_per/`. Post-processed with 
   `renumber_after_cdhit.py`.

5. **`DIAMOND_BLAST/`** (via `Pipeline/diamond_blast.sh`)
   - `flag_suspicious_data.py` — flags non-secreted sequences with high 
     homology to secreted classes (very-sensitive DIAMOND BLAST mode)
   - `remove_flagged_sequence.py` — removes flagged sequences
   - `renumber.py` — renumbers/relabels sequences post-cleaning
   - `Copy_CleanData_into_Clean_data.py` — copies the final cleaned 
     output into `Results/.../Clean_Data/`

6. **`Train_Test_Split/Formula_perform.ipynb`**
   Produces the final, leakage-checked train/test split → 
   `Results/.../Train_Test/train.fasta` and `test.fasta`.

**Output of Phase 0:** deduplicated, cross-contamination-checked dataset, 
split into train/test FASTA files, ready for feature engineering.

---

## 6. Pipeline — Phase 1 (Feature Engineering) — in progress

- **iFeatureOmega** (`feature_engineering/iFeatureOmega/feature_extraction.ipynb`) 
  — extracts sequence-based features (e.g. AAC, DPC, CTD)
- **POSSUM** (`feature_engineering/POSSUM/`) — extracts PSSM-based 
  evolutionary features
- **combine/** — merges feature sets from both tools into a single 
  feature table per class, stored in `Results/.../phase1_feature_engineering/processed/`

---

## 7. Environment Setup

```bash
python -m venv venv
source venv/Scripts/activate      # Windows Git Bash
pip install -r requirements.txt
```

Editable install of reusable code (`src/secretion_prediction`):
```bash
pip install -e .
```

Bioinformatics tools used outside pip (install separately):
- CD-HIT — [version : 4.8.1]
- DIAMOND BLAST — [version: 2.2.5]
- MMseqs2 — [version : 18.8.cc5c]
- iFeatureOmega-CLI — [fill: source/version, if re-added later]

---

## 8. Current Status

- [x] Phase 0 — Dataset combination, cleaning (CD-HIT + DIAMOND), 
      leakage-checked train/test split — **COMPLETE**
- [ ] Phase 1 — Feature engineering (iFeatureOmega) — in progress
- [ ] Phase 2 — Classical ML baseline (Random Forest, XGBoost, SVM)
- [ ] Phase 3 — ESM2 embeddings  (RF, XGBoost, SVM)
- [ ] Final evaluation on independent test sets

---

## 9. Note on MMseqs2

A separate clustering pipeline using MMseqs2 
(`Results/.../phase0_dataset_prep/MMsques2/`) is maintained alongside the 
primary CD-HIT + DIAMOND pipeline. It clusters sequences by homology at 
multiple identity thresholds (30%, 40%, 50%, all at 80% coverage) for 
comparison purposes, and produces its own train/test split 
(`MMsques2/Re-configure/train.csv`, `test.csv`).

**The primary train/test split used for downstream modeling is 
`Train_Test/train.fasta` / `test.fasta`**, not the MMseqs2 output, 
unless stated otherwise in a specific analysis.

---

## 10. Conventions

- Raw data (`Data/.../raw/`) is treated as read-only; all derived 
  outputs go under `Results/`.
- FASTA files (`.fasta`) are tracked with `-text` in `.gitattributes` 
  to preserve original LF line endings and avoid CRLF corruption on 
  Windows.
- Folder/file names avoid special characters (e.g. `%`) for shell and 
  cross-platform compatibility (`CD_HIT_60_per`, not `CD_HIT_60%`).
- Large or regeneratable binary files (`.dmnd` databases, model 
  checkpoints, embeddings) are excluded via `.gitignore` rather than 
  committed, unless small enough to track directly (e.g. MMseqs2's 
  `combined_db`, ~1.5MB).
- All path references in code use `src/secretion_prediction/config.py`, 
  which resolves paths relative to the repository root (found via the 
  `.git` marker) — no hardcoded absolute paths.