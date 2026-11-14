# Feature Engineering — iFeatureOmega

## Purpose
Extracts AAC, DPC, N-terminal window (30 residues), and C-terminal 
window (30 residues) features from protein sequences in train.csv
and test.csv

## Why a separate environment
iFeatureOmegaCLI requires Python 3.7 (per its own GitHub README). 
The main project environment runs Python 3.14, incompatible. 
This tool therefore runs in an isolated conda environment, `ifeature_env`.

## Version 
Python = 3.7.1
iFeatureOmegaCLI= 1.0.2

## Setup (one-time)
```bash
conda env create -f environment_ifeature_env.yml
conda activate ifeature_env
```

## Usage
```bash
conda activate ifeature_env
python extract_features.py
conda deactivate
```

## Output
Writes 8 CSV files to `Results/sumit_secretion_prediction/phase1_feature_engineering/features/`:
train/test × AAC, DPC, Nterm, Cterm

## Source
- Repo: https://github.com/Superzchen/iFeatureOmega-CLI
- Citation: Chen et al., iFeatureOmega: an integrative platform for
engineering, visualization and analysis of features
from molecular sequences, structural and ligand data
sets

## Known issues encountered
- `rdkit` import required at package load time even though ligand 
  descriptors (iLigand class) are unused — undocumented dependency, 
  install via `conda install -c rdkit rdkit` if ModuleNotFoundError occurs.
- 3 sequences (T1SS-97, T3SS-361, T6SS-231) contained non-standard 
  character 'X' — removed before extraction.