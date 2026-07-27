#!/usr/bin/env python3
"""Build a binary feature matrix from the AMR_genotypes field.

Each distinct gene entry across the cohort becomes one column. An isolate
scores 1 in that column when AMRFinderPlus reported the entry for it.
This is multi-hot encoding: one categorical set turned into many binary
columns, which is the standard way to hand set-valued data to a model.

Usage:
    python scripts/06_build_features.py [drug]

Reads  data/interim/cohort_<drug>.parquet
Writes data/processed/features_<drug>.parquet
       data/processed/labels_<drug>.parquet
"""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"

# Drop features carried by fewer than this many isolates. A column with
# two positives cannot support a reliable estimate of its effect, and it
# gives a flexible model something to memorise.
MIN_COUNT = int(sys.argv[2]) if len(sys.argv) > 2 else 3

# Prefixes for genes with an established mechanistic link to a drug in
# this chapter. Dropping one of these is worth a warning, because the
# rarity threshold is an arbitrary choice and the mechanism is not.
MECHANISTIC = ("gyrA_", "parC_", "tet(", "23S_", "erm(")


def entries(genotypes):
    """Split an AMR_genotypes field into individual gene entries."""
    if not isinstance(genotypes, str) or genotypes in ("", "NULL"):
        return []
    return [e.strip() for e in genotypes.strip('"').split(",") if e.strip()]


def main():
    drug = sys.argv[1] if len(sys.argv) > 1 else "ciprofloxacin"
    slug = drug.replace(" ", "_")

    path = INTERIM / f"cohort_{slug}.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 02_build_cohort.py {drug!r} first.")

    df = pd.read_parquet(path)
    print(f"[drug] {drug}")
    print(f"[read] {len(df):,} isolates")

    gene_sets = df["AMR_genotypes"].apply(entries)

    # Every distinct entry seen anywhere in the cohort.
    vocabulary = sorted({g for gs in gene_sets for g in gs})
    print(f"[feat] {len(vocabulary):,} distinct gene entries before filtering")

    # Build the matrix one column at a time.
    X = pd.DataFrame(
        {g: gene_sets.apply(lambda gs, g=g: int(g in gs)) for g in vocabulary},
        index=df.index,
    )

    counts = X.sum()
    keep = counts[counts >= MIN_COUNT].index
    dropped = len(vocabulary) - len(keep)
    X = X[keep]
    print(f"[feat] dropped {dropped:,} entries seen in fewer than {MIN_COUNT} isolates")
    print(f"[feat] {X.shape[1]:,} features retained "
          f"(threshold {MIN_COUNT})")

    lost = [g for g in vocabulary
            if g not in keep and g.startswith(MECHANISTIC)]
    if lost:
        print("[warn] dropped features with a known mechanism:")
        for g in lost:
            print(f"         {int(counts[g]):5d}  {g}")

    X.insert(0, "asm_acc", df["asm_acc"].values)

    y = pd.DataFrame({
        "asm_acc": df["asm_acc"].values,
        "resistant": (df[drug] == "R").astype(int).values,
        "snp_cluster": df["snp_cluster"].values,
    })

    assert len(X) == len(y), "feature and label row counts differ"
    assert (X["asm_acc"].values == y["asm_acc"].values).all(), \
        "feature and label rows are not aligned"

    PROCESSED.mkdir(parents=True, exist_ok=True)
    X.to_parquet(PROCESSED / f"features_{slug}.parquet", index=False)
    y.to_parquet(PROCESSED / f"labels_{slug}.parquet", index=False)

    n, p = len(X), X.shape[1] - 1
    density = X[keep].to_numpy().mean()

    print()
    print(f"[done] X: {n:,} rows x {p} features")
    print(f"       y: {int(y['resistant'].sum()):,} resistant, "
          f"{int((1 - y['resistant']).sum()):,} susceptible")
    print(f"       samples per feature: {n / p:.0f}")
    print(f"       matrix density: {density:.3f} "
          f"({density * 100:.1f}% of cells are 1)")
    print()
    print("       most common features:")
    for g, c in counts[keep].sort_values(ascending=False).head(8).items():
        print(f"         {c:5d}  {g}")


if __name__ == "__main__":
    main()
