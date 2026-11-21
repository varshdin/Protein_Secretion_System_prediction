r"""
Feature engineering pipeline - phase 1.

Runs every feature block over the train and test splits, in order:

    1. aac   20-dim amino acid composition     extract_feature.py  (iFeatureOmega)
    2. dpc  400-dim dipeptide composition      extract_feature.py  (iFeatureOmega)
    3. ncw   ESM-2 N-/C-terminal window 50/50  extract_ncw.py      (transformers)

iFeatureOmega needs Python 3.7 and transformers needs Python 3.10+, so the two
cannot share an environment. That incompatibility is the only reason this file
exists. Each block names the conda environment it needs and is launched with
`conda run`, which switches environment for one command and switches back. A
block whose env is None runs with the current interpreter.

Run it from any shell where conda is on PATH - Git Bash and Anaconda Prompt
both work on Windows. Not from WSL: the conda environments are on the Windows
side. Phase 0 (CD-HIT, DIAMOND, MMseqs2) is the opposite case and runs in WSL,
because those are Linux-only tools. That is why the two phases are separate
commands rather than one chained script.

    python feature_engineering.py                    # all blocks, train + test
    python feature_engineering.py --dry-run          # print the commands only
    python feature_engineering.py --blocks ncw
    python feature_engineering.py --splits test


Environments a fresh clone needs:

    conda create -n ifeature_env python=3.7
    conda activate ifeature_env
    pip install iFeatureOmegaCLI

    conda create -n esm_env python=3.11
    conda activate esm_env
    pip install torch --index-url https://download.pytorch.org/whl/cpu
    pip install transformers pandas

The ncw block downloads facebook/esm2_t6_8M_UR50D from Hugging Face on first
use and caches it, so that run needs network access.

All paths resolve against the repository root, so this runs from any working
directory on any machine that clones the repo.
"""

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path


# --------------------------------------------------------------------------
# CONFIG - the only part you should need to edit.
# Set an env to None to run that block with the current interpreter.
# --------------------------------------------------------------------------
ENV_IFEATURE = "ifeature_env"      # Python 3.7, has iFeatureOmegaCLI
ENV_ESM = "esm_env"
WINDOW = "50,50"                   # selected by CV on the training set
# --------------------------------------------------------------------------


SCRIPT_DIR = Path(__file__).resolve().parent

# This file belongs at:
#   <repo>/Scripts/sumit_secretion_prediction/Pipeline/
# beside the other phase pipelines, so the repository root is three levels up.
# Anchoring to __file__ rather than the working directory is what makes every
# path below portable.
_PARENTS = Path(__file__).resolve().parents
if len(_PARENTS) <= 3:
    raise SystemExit(
        "this file sits at %s, only %d level(s) below the filesystem root.\n"
        "  It belongs at <repo>/Scripts/sumit_secretion_prediction/Pipeline/,\n"
        "  three levels below the repository root."
        % (SCRIPT_DIR, len(_PARENTS) - 1)
    )
PROJECT_ROOT = _PARENTS[3]

# The block scripts live in feature_engineering/, not beside this pipeline.
# Resolved from PROJECT_ROOT so moving this file between Pipeline/ and
# feature_engineering/ does not change where the blocks are looked up.
BLOCK_DIR = (PROJECT_ROOT / "Scripts" / "sumit_secretion_prediction"
             / "feature_engineering")

FASTA_DIR = (PROJECT_ROOT / "Results" / "sumit_secretion_prediction"
             / "phase0_dataset_prep" / "Train_Test")

DEFAULT_FEATURE_DIR = (PROJECT_ROOT / "Results" / "sumit_secretion_prediction"
                       / "phase1_feature_engineering" / "features")

# Each block names the environment it needs, the script that implements it,
# and the arguments that select the block. Output paths are deliberately NOT
# listed: each script reports what it wrote, which is the truth rather than a
# prediction this file would have to keep in sync.
BLOCKS = {
    "aac": {
        "env": ENV_IFEATURE,
        "script": BLOCK_DIR / "iFeatureOmega" / "extract_feature.py",
        "args": ["--block", "aac"],
        "subdir": "AAC",
    },
    "dpc": {
        "env": ENV_IFEATURE,
        "script": BLOCK_DIR / "iFeatureOmega" / "extract_feature.py",
        "args": ["--block", "dpc"],
        "subdir": "DPC",
    },
    "ncw": {
        "env": ENV_ESM,
        "script": BLOCK_DIR / "ESM-2" / "extract_ncw.py",
        "args": ["--window", WINDOW],
        "subdir": "N_C_window",
    },
}

# Explicit, so the log reads in the intended order rather than alphabetically.
BLOCK_ORDER = ["aac", "dpc", "ncw"]
SPLITS = ["train", "test"]


def check_project_root():
    """Stop if PROJECT_ROOT is wrong rather than writing to the wrong place."""
    if not (PROJECT_ROOT / "Results").is_dir():
        raise SystemExit(
            "repository root looks wrong: %s\n"
            "  expected a 'Results' directory there. Has this script been\n"
            "  moved or copied out of Scripts/sumit_secretion_prediction/\n"
            "  Pipeline/?" % PROJECT_ROOT
        )


def summarise(feature_dir, since):
    """Report every CSV written during this run.

    Files are found by modification time rather than by expected name, so this
    reports what is actually on disk instead of what the pipeline predicted -
    and it ignores older files in the same folders, such as the rejected
    window conditions kept from the selection experiment.

    No pandas: this file may run in an interpreter that does not have it.
    """
    if not feature_dir.is_dir():
        print("no features directory at %s" % feature_dir)
        return

    written = sorted(p for p in feature_dir.rglob("*.csv")
                     if p.stat().st_mtime >= since)
    if not written:
        print("no CSV files were written")
        return

    print("%-46s %7s %7s  %s" % ("file", "rows", "cols", "id column"))
    print("-" * 72)

    id_columns = {}
    for path in written:
        with open(str(path)) as handle:
            header = handle.readline().rstrip("\n").split(",")
            n_rows = sum(1 for _ in handle)
        id_columns.setdefault(header[0], []).append(path.name)
        print("%-46s %7d %7d  %s"
              % (path.relative_to(feature_dir), n_rows, len(header) - 1,
                 header[0]))

    # Every block must name its id column identically or the phase 2 merge
    # fails. This has gone wrong once already.
    if len(id_columns) > 1:
        print("\nWARNING: the id column is not named the same in every file:")
        for name, files in sorted(id_columns.items()):
            print("  %-14s %s" % (name, ", ".join(files)))
        print("  a phase 2 merge on the id column will fail - fix this first.")

    # Row counts must agree within a split, or the blocks saw different input.
    by_split = {}
    for path in written:
        for split in SPLITS:
            if split in path.stem:
                with open(str(path)) as handle:
                    handle.readline()
                    by_split.setdefault(split, {})[path.name] = sum(1 for _ in handle)
    for split, counts in sorted(by_split.items()):
        if len(set(counts.values())) > 1:
            print("\nWARNING: %s row counts disagree: %s" % (split, counts))


def build_command(name, split, fasta_dir, feature_dir):
    """The command that runs one block on one split.

    --fasta is passed as a full path and --out as an explicit directory, so the
    block scripts never have to work out where anything lives. That is what
    lets --out-root redirect the output of a trial run.
    """
    block = BLOCKS[name]
    inner = [str(block["script"]),
             "--fasta", str(fasta_dir / ("%s.fasta" % split)),
             "--out", str(feature_dir / block["subdir"])] + block["args"]

    if block["env"]:
        # `conda run` activates the environment for this one command only, so
        # the parent shell is never touched and no `conda activate` is needed.
        # --no-capture-output lets the block's own progress reach the terminal.
        return ["conda", "run", "--no-capture-output",
                "-n", block["env"], "python"] + inner
    return [sys.executable] + inner


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--blocks", nargs="+", default=BLOCK_ORDER,
                        choices=BLOCK_ORDER, help="blocks to run")
    parser.add_argument("--splits", nargs="+", default=SPLITS,
                        choices=SPLITS, help="splits to run")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the commands without running anything")
    parser.add_argument("--out-root", default=None, metavar="DIR",
                        help="write features under DIR instead of the repo, for"
                             " a trial run that cannot touch existing output")
    args = parser.parse_args()

    check_project_root()
    blocks = [b for b in BLOCK_ORDER if b in args.blocks]

    # Input always comes from the real repository; only the destination moves.
    # A trial run therefore exercises the real data without being able to
    # overwrite anything under Results/.
    feature_dir = (Path(args.out_root).resolve() if args.out_root
                   else DEFAULT_FEATURE_DIR)

    print("repo   : %s" % PROJECT_ROOT)
    print("fasta  : %s" % FASTA_DIR)
    print("output : %s%s" % (feature_dir,
                             "   (TRIAL RUN)" if args.out_root else ""))
    print("blocks : %s" % ", ".join(blocks))
    print("splits : %s\n" % ", ".join(args.splits))

    # Fail before launching anything rather than halfway through. SystemExit
    # rather than a raise, so a wrong path reads as a message and not as a
    # traceback that looks like a bug in the pipeline itself.
    for split in args.splits:
        path = FASTA_DIR / ("%s.fasta" % split)
        if not path.exists():
            raise SystemExit(
                "missing input FASTA: %s\n"
                "  phase 0 should have written %s.fasta to the Train_Test\n"
                "  folder. Run the MMseqs2 split stage first." % (path, split)
            )
    for name in blocks:
        script = BLOCKS[name]["script"]
        if not script.exists():
            raise SystemExit(
                "missing script for block '%s': %s\n"
                "  check the folder name in the BLOCKS section." % (name, script)
            )

    # If any block needs a conda environment, conda itself must be callable.
    # Without this the failure arrives as a bare FileNotFoundError from
    # subprocess, which says nothing useful.
    # Skipped for --dry-run, whose whole purpose is to show the commands on a
    # machine that may not be set up to run them yet.
    if (not args.dry_run
            and any(BLOCKS[n]["env"] for n in blocks)
            and shutil.which("conda") is None):
        raise SystemExit(
            "a block needs a conda environment but 'conda' is not on PATH.\n"
            "  open an Anaconda Prompt, or set that block's env to None in\n"
            "  the CONFIG section to use the current interpreter."
        )

    results = []
    failed = False
    started = time.time()

    for name in blocks:
        for split in args.splits:
            label = "%s / %s" % (name, split)
            cmd = build_command(name, split, FASTA_DIR, feature_dir)
            where = BLOCKS[name]["env"] or "current interpreter"

            if args.dry_run:
                print("would run  %-14s [%s]" % (label, where))
                print("           %s\n" % " ".join(cmd))
                results.append((label, "dry-run", 0.0))
                continue

            print("=" * 72)
            print("running %s   [%s]" % (label, where))
            print("=" * 72)

            start = time.time()
            try:
                returncode = subprocess.run(cmd).returncode
            except (FileNotFoundError, OSError) as exc:
                print("\nFAILED %s: could not launch %r\n       %s"
                      % (label, cmd[0], exc))
                results.append((label, "FAILED", time.time() - start))
                failed = True
                break

            elapsed = time.time() - start
            if returncode != 0:
                print("\nFAILED %s (exit %d) after %.1fs" % (label, returncode, elapsed))
                results.append((label, "FAILED", elapsed))
                failed = True
                break

            print("\nok %s in %.1fs\n" % (label, elapsed))
            results.append((label, "ok", elapsed))

        if failed:
            break

    print("\n" + "=" * 72)
    print("%-16s %-10s %8s" % ("block / split", "status", "seconds"))
    print("-" * 72)
    for label, status, elapsed in results:
        print("%-16s %-10s %8.1f" % (label, status, elapsed))

    if failed:
        # Stopping on the first failure is deliberate. A half-built feature
        # set that looks complete is worse than an obvious stop.
        print("\n1 step failed - nothing after it was run")
        return 1

    print("\nall steps completed")
    if not args.dry_run:
        print("\n" + "=" * 72)
        print("FEATURE MATRICES WRITTEN")
        print("=" * 72)
        summarise(feature_dir, started)
    return 0


if __name__ == "__main__":
    sys.exit(main())