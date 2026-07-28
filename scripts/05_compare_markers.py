#!/usr/bin/env python3
"""Compare candidate feature definitions for a single-marker rule.

Feature engineering: the same data and the same model, with the feature
defined to match the mechanism rather than one instance of it.

Usage:
    python scripts/05_compare_markers.py [drug]

Reads  data/interim/cohort_<drug>.parquet
Writes results/metrics/marker_comparison_<drug>.txt
"""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
METRICS = ROOT / "results" / "metrics"

# Candidate feature definitions, narrowest first. Each is a list of
# prefixes; an isolate scores positive if any AMR_genotypes entry
# starts with any prefix in the list.
CANDIDATES = {
    "ciprofloxacin": [
        ("gyrA_T86I only",        ["gyrA_T86I"]),
        ("any change at gyrA 86", ["gyrA_T86"]),
        ("any gyrA change",       ["gyrA_"]),
    ],
    "nalidixic acid": [
        ("gyrA_T86I only",        ["gyrA_T86I"]),
        ("any change at gyrA 86", ["gyrA_T86"]),
        ("any gyrA change",       ["gyrA_"]),
    ],
    "tetracycline": [
        ("tet(O) only",           ["tet(O)"]),
        ("any tet gene",          ["tet("]),
    ],
}


def entries(genotypes):
    """Split an AMR_genotypes field into individual gene entries."""
    if not isinstance(genotypes, str) or genotypes in ("", "NULL"):
        return []
    return [e.strip() for e in genotypes.strip('"').split(",") if e.strip()]


def evaluate(df, drug, prefixes):
    """Score one feature definition. Returns tp, fn, fp, tn."""
    flag = df["entries"].apply(
        lambda es: any(e.startswith(p) for e in es for p in prefixes)
    )
    resistant = df[drug] == "R"
    return (
        int((flag & resistant).sum()),
        int((~flag & resistant).sum()),
        int((flag & ~resistant).sum()),
        int((~flag & ~resistant).sum()),
    )


def main():
    if len(sys.argv) < 2:
        print(f"usage: 05_compare_markers.py <drug>   drugs: {', '.join(CANDIDATES)}",
              file=sys.stderr)
        sys.exit(2)
    drug = sys.argv[1]
    slug = drug.replace(" ", "_")

    if drug not in CANDIDATES:
        sys.exit(f"[fail] no candidates defined for '{drug}'.")

    path = INTERIM / f"cohort_{slug}.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 02_build_cohort.py {drug!r} first.")

    df = pd.read_parquet(path)
    df["entries"] = df["AMR_genotypes"].apply(entries)

    n = len(df)
    n_res = int((df[drug] == "R").sum())
    baseline_err = min(n_res, n - n_res)

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Drug:     {drug}")
    out(f"Cohort:   {n:,} isolates ({n_res:,} resistant)")
    out(f"Baseline: always guess the majority class")
    out(f"          {baseline_err:,} errors, accuracy {1 - baseline_err/n:.4f}")
    out()

    header = f"{'feature definition':<26}{'TP':>6}{'FN':>6}{'FP':>6}{'TN':>7}" \
             f"{'sens':>8}{'spec':>8}{'acc':>8}{'err':>6}"
    out(header)
    out("-" * len(header))

    for name, prefixes in CANDIDATES[drug]:
        tp, fn, fp, tn = evaluate(df, drug, prefixes)
        sens = tp / (tp + fn) if tp + fn else float("nan")
        spec = tn / (tn + fp) if tn + fp else float("nan")
        acc = (tp + tn) / n
        out(f"{name:<26}{tp:>6}{fn:>6}{fp:>6}{tn:>7}"
            f"{sens:>8.3f}{spec:>8.3f}{acc:>8.4f}{fp+fn:>6}")

    out()
    out("Note: these definitions were chosen on mechanistic grounds before")
    out("scoring, not selected by whichever performed best. Choosing a")
    out("feature by its score on the full dataset is a form of overfitting,")
    out("and the reported accuracy would no longer be honest.")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"marker_comparison_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
