#!/usr/bin/env python3
"""Describe the Chapter 1 cohort before any modelling.

Reports class balance, SNP cluster structure, and the agreement between
the AMRFinderPlus tet(O) call and the measured tetracycline phenotype.

Reads  data/interim/cohort.parquet
Writes results/metrics/cohort_summary.txt
"""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
METRICS = ROOT / "results" / "metrics"

DRUG = "tetracycline"


def has_teto(genotypes):
    """True when AMRFinderPlus reported a tet(O) family gene.

    The genotype field looks like:
        "blaOXA-193,gyrA_T86I=POINT,tet(O)"

    Acquired genes appear as bare names; point mutations carry =POINT.
    Matching on 'tet(' catches tet(O), tet(O/32/O) and other variants.
    """
    if not isinstance(genotypes, str):
        return False
    return any(g.strip().startswith("tet(") for g in genotypes.strip('"').split(","))


def main():
    path = INTERIM / "cohort.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 02_build_cohort.py first.")

    df = pd.read_parquet(path)
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Cohort: {len(df):,} isolates")
    out()

    out("== Class balance ==")
    counts = df[DRUG].value_counts()
    out(counts.to_string())
    majority = counts.max() / len(df)
    out(f"majority-class baseline accuracy: {majority:.3f}")
    out()

    out("== SNP cluster structure ==")
    clustered = df[df["snp_cluster"].notna()]
    sizes = clustered["snp_cluster"].value_counts()
    out(f"isolates with a cluster: {len(clustered):,}")
    out(f"isolates without:        {len(df) - len(clustered):,}")
    out(f"distinct clusters:       {sizes.nunique():,}")
    out(f"largest cluster:         {sizes.max():,} isolates")
    out(f"median cluster size:     {sizes.median():.0f}")
    out(f"singleton clusters:      {(sizes == 1).sum():,}")
    out("largest five clusters:")
    out(sizes.head(5).to_string())
    out()

    out("== Genotype vs measured phenotype ==")
    df["teto"] = df["AMR_genotypes"].apply(has_teto)
    table = pd.crosstab(df["teto"], df[DRUG])
    out(table.to_string())
    out()

    if "R" in table.columns and True in table.index:
        tp = table.loc[True, "R"]
        fn = table.loc[False, "R"] if False in table.index else 0
        fp = table.loc[True, "S"] if "S" in table.columns else 0
        tn = table.loc[False, "S"] if False in table.index else 0

        out(f"tet(O) present and resistant     (TP): {tp:,}")
        out(f"tet(O) absent  but resistant     (FN): {fn:,}")
        out(f"tet(O) present but susceptible   (FP): {fp:,}")
        out(f"tet(O) absent  and susceptible   (TN): {tn:,}")
        out()
        out(f"sensitivity: {tp / (tp + fn):.3f}")
        out(f"specificity: {tn / (tn + fp):.3f}")
        out(f"accuracy:    {(tp + tn) / len(df):.3f}")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / "cohort_summary.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
