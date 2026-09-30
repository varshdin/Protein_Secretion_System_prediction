from pathlib import Path

def find_repo_root(marker=".git"):
    """
    Walk UP from this file's location until a folder containing 
    `.git` is found — that folder is the repo root.
    Works no matter how deep this file, or any script importing it, 
    is located, and works on ANY machine (Windows, Linux, Mac).
    """
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / marker).exists():
            return parent
    raise RuntimeError(f"Repo root not found — no '{marker}' folder in any parent directory")

# ---- Core root ----
PROJECT_ROOT = find_repo_root()

# ---- Data ----
DATA_DIR = PROJECT_ROOT / "Data" / "sumit_secretion_prediction"
RAW_DIR = DATA_DIR / "raw"
RAW_DEEPSECE = RAW_DIR / "deepsecE"
RAW_PLM_EFFECTOR = RAW_DIR / "PLM-Effector"
RAW_SIGNALP = RAW_DIR / "SignalP"
RAW_T5SE = RAW_DIR / "T5SE"
RAW_COMBINED = RAW_DIR / "Combined_raw_dataset"

# ---- Results ----
RESULTS_DIR = PROJECT_ROOT / "Results" / "sumit_secretion_prediction"
PHASE0_DIR = RESULTS_DIR / "phase0_dataset_prep"
CDHIT_DIR = PHASE0_DIR / "CD_HIT_60_per"
DIAMOND_DIR = PHASE0_DIR / "DIAMOND_BLAST"
CLEAN_DATA_DIR = PHASE0_DIR / "Clean_Data"
MMSEQS2_DIR = PHASE0_DIR / "MMsques2" / "Re-configure"
TRAIN_TEST_DIR = PHASE0_DIR / "Train_Test"

PHASE1_DIR = RESULTS_DIR / "phase1_feature_engineering"
FEATURES_DIR = PHASE1_DIR / "features"
PROCESSED_DIR = PHASE1_DIR / "processed"

PHASE2_DIR = RESULTS_DIR / "phase2_classical_ml"
PHASE3_DIR = RESULTS_DIR / "phase3_esm_embeddings"
FINAL_EVAL_DIR = RESULTS_DIR / "final_evaluation"

# ---- Scripts ----
SCRIPTS_DIR = PROJECT_ROOT / "Scripts" / "sumit_secretion_prediction"

# ---- Models ----
MODELS_DIR = PROJECT_ROOT / "Models" / "sumit_secretion_prediction"
MODELS_CLASSICAL = MODELS_DIR / "classical"
MODELS_DL = MODELS_DIR / "deep_learning"

# ---- Figures ----
FIGURES_DIR = PROJECT_ROOT / "Figures" / "sumit_secretion_prediction"