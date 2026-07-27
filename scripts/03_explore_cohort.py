#!/usr/bin/env python3
"""Describe a Chapter 1 cohort before any modelling.

Reports class balance, SNP cluster structure, and the agreement between a
single AMRFinderPlus marker and the measured phenotype for one drug.

Usage:
    python scripts/03_explore_cohort.py [drug]

Reads  data/interim/cohort_<drug>.parquet
Writes results/metrics/cohort_summary_<drug>.txt
"""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
METRICS = ROOT / "results" / "metrics"

MARKERS = {
    "tetracycline": "tet(",
    "ciprofloxacin": "gyrA_T86I",
    "nalidixic acid": "gyrA_T86I",
    "erythromycin": "23S_",
    "azithromycin": "23S_",
}


def has_marker(genotypes, prefix):
    """True when AMR_genotypes contains an entry starting with prefix.

    The field looks like:
        "blaOXA-193,gyrA_T86I=POINT,tet(O)"

    Acquired genes appear as bare names; point mutations carry =POINT.
    Prefix matching catches tet(O), tet(O/32/O), gyrA_T86I=POINT, etc.
    """
    if not isinstance(genotypes, str):
        return False
    entries = genotypes.strip('"').split(",")
    return any(e.strip().startswith(prefix) for e in entries)


def main():
    drug = sys.argv[1] if len(sys.argv) > 1 else "tetracycline"
    slug = drug.replace(" ", "_")

    if drug not in MARKERS:
        sys.exit(f"[fail] no marker defined for '{drug}'. Known: {list(MARKERS)}")
    marker = MARKERS[drug]

    path = INTERIM / f"cohort_{slug}.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 02_build_cohort.py {drug!r} first.")

    df = pd.read_parquet(path)
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Drug:   {drug}")
    out(f"Marker: {marker}")
    out(f"Cohort: {len(df):,} isolates")
    out()

    out("== Class balance ==")
    counts = df[drug].value_counts()
    out(counts.to_string())
    out(f"majority-class baseline accuracy: {counts.max() / len(df):.3f}")
    out()

    out("== SNP cluster structure ==")
    clustered = df[df["snp_cluster"].notna()]
    sizes = clustered["snp_cluster"].value_counts()
    out(f"isolates with a cluster: {len(clustered):,}")
    out(f"isolates without:        {len(df) - len(clustered):,}")
    out(f"distinct clusters:       {len(sizes):,}")
    out(f"largest cluster:         {sizes.max():,} isolates")
    out(f"median cluster size:     {sizes.median():.0f}")
    out(f"singleton clusters:      {(sizes == 1).sum():,}")
    out("largest five clusters:")
    out(sizes.head(5).to_string())
    out()

    out("== Marker vs measured phenotype ==")
    df["marker"] = df["AMR_genotypes"].apply(lambda g: has_marker(g, marker))
    table = pd.crosstab(df["marker"], df[drug])
    out(table.to_string())
    out()

    def cell(row, col):
        if row in table.index and col in table.columns:
            return int(table.loc[row, col])
        return 0

    tp, fn = cell(True, "R"), cell(False, "R")
    fp, tn = cell(True, "S"), cell(False, "S")

    out(f"marker present, resistant    (TP): {tp:,}")
    out(f"marker absent,  resistant    (FN): {fn:,}")
    out(f"marker present, susceptible  (FP): {fp:,}")
    out(f"marker absent,  susceptible  (TN): {tn:,}")
    out()

    total = tp + fn + fp + tn
    if tp + fn:
        out(f"sensitivity: {tp / (tp + fn):.3f}")
    if tn + fp:
        out(f"specificity: {tn / (tn + fp):.3f}")
    if total:
        out(f"accuracy:    {(tp + tn) / total:.3f}")
        out(f"errors:      {fp + fn:,} of {total:,}")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"cohort_summary_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
