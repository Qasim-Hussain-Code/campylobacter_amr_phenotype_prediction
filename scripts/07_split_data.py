#!/usr/bin/env python3
"""Assign every isolate to a train or test fold, two ways.

random    -- a stratified random split, the textbook default
grouped   -- SNP clusters kept whole, so no lineage spans the split

Splits are written to disk once and read by every model script, which makes
it structurally impossible to compare two models on different splits.

Usage:
    python scripts/07_split_data.py [drug]

Reads  data/processed/labels_<drug>.parquet
Writes data/processed/splits_<drug>.parquet
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GroupShuffleSplit

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"

SEED = 42
TEST_FRACTION = 0.25


def main():
    drug = sys.argv[1] if len(sys.argv) > 1 else "ciprofloxacin"
    slug = drug.replace(" ", "_")

    path = PROCESSED / f"labels_{slug}.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 06_build_features.py {drug!r} first.")

    y = pd.read_parquet(path)
    n = len(y)
    print(f"[drug] {drug}")
    print(f"[read] {n:,} isolates, {int(y['resistant'].sum()):,} resistant")

    idx = np.arange(n)

    # Random split. Stratified so both folds hold the same class ratio;
    # without it a small test set can drift away from the base rate by
    # chance and the accuracy becomes hard to interpret.
    train_i, test_i = train_test_split(
        idx,
        test_size=TEST_FRACTION,
        random_state=SEED,
        stratify=y["resistant"],
    )
    y["split_random"] = "train"
    y.loc[test_i, "split_random"] = "test"

    # Grouped split. Isolates with no cluster assignment are treated as
    # their own singleton group, which is the conservative reading: we
    # cannot show they are related to anything, so we do not assume it.
    groups = y["snp_cluster"].fillna(
        pd.Series([f"__unclustered_{i}" for i in idx], index=y.index)
    )

    gss = GroupShuffleSplit(n_splits=1, test_size=TEST_FRACTION, random_state=SEED)
    train_g, test_g = next(gss.split(idx, y["resistant"], groups))
    y["split_grouped"] = "train"
    y.loc[test_g, "split_grouped"] = "test"

    for name in ("split_random", "split_grouped"):
        tr = y[y[name] == "train"]
        te = y[y[name] == "test"]
        print()
        print(f"[{name.split('_')[1]:>7}] train {len(tr):,}  test {len(te):,}")
        print(f"          train resistant {tr['resistant'].mean():.3f}  "
              f"test resistant {te['resistant'].mean():.3f}")

        shared = set(tr["snp_cluster"].dropna()) & set(te["snp_cluster"].dropna())
        print(f"          clusters spanning both folds: {len(shared):,}")

    PROCESSED.mkdir(parents=True, exist_ok=True)
    out = PROCESSED / f"splits_{slug}.parquet"
    y.to_parquet(out, index=False)
    print(f"\n[done] wrote {out.relative_to(ROOT)}  seed={SEED}")


if __name__ == "__main__":
    main()
