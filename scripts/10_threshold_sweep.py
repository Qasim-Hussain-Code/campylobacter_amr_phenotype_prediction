#!/usr/bin/env python3
"""Test whether the Chapter 1 conclusion depends on the rarity threshold.

MIN_COUNT in 06_build_features.py is an arbitrary choice. A conclusion that
moves when it moves is a conclusion about the threshold, not the biology.

The sweep also tracks which features are admitted at each setting. Rare
features are where passengers live: a gene carried by five isolates that
happen to be resistant looks causal and is not.

Usage:
    python scripts/10_threshold_sweep.py [drug]

Reads  data/interim/cohort_<drug>.parquet
Writes results/metrics/threshold_sweep_<drug>.txt
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
METRICS = ROOT / "results" / "metrics"

THRESHOLDS = [1, 3, 5, 10, 20]
N_SPLITS = 5
N_REPEATS = 3

RULE_PREFIX = {
    "ciprofloxacin": "gyrA_T86",
    "nalidixic acid": "gyrA_T86",
    "tetracycline": "tet(",
}

# Prefixes with an accepted mechanistic link to some drug in this chapter.
# Anything outside this set that tracks the rule perfectly is a passenger.
MECHANISTIC = ("gyrA_", "parC_", "tet(", "23S_", "erm(")


def entries(genotypes):
    if not isinstance(genotypes, str) or genotypes in ("", "NULL"):
        return []
    return [e.strip() for e in genotypes.strip('"').split(",") if e.strip()]


def cv_errors(Xv, y, groups, scheme, seed_base):
    """Mean error rate for the rule and for logistic regression."""
    rule, log = [], []
    for seed in range(N_REPEATS):
        if scheme == "random":
            cv = StratifiedKFold(N_SPLITS, shuffle=True, random_state=seed)
            folds = cv.split(Xv["all"], y)
        else:
            cv = StratifiedGroupKFold(N_SPLITS, shuffle=True, random_state=seed)
            folds = cv.split(Xv["all"], y, groups)
        for tr, te in folds:
            rule.append((Xv["rule"][te] != y[te]).mean())
            clf = LogisticRegression(max_iter=2000, random_state=seed)
            clf.fit(Xv["all"][tr], y[tr])
            log.append((clf.predict(Xv["all"][te]) != y[te]).mean())
    return np.array(rule), np.array(log)


def main():
    drug = sys.argv[1] if len(sys.argv) > 1 else "ciprofloxacin"
    slug = drug.replace(" ", "_")

    path = INTERIM / f"cohort_{slug}.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 02_build_cohort.py {drug!r} first.")

    df = pd.read_parquet(path)
    gene_sets = df["AMR_genotypes"].apply(entries)
    y = (df[drug] == "R").to_numpy().astype(int)

    g = df["snp_cluster"].copy()
    solo = g.isna()
    g[solo] = [f"__solo_{i}" for i in np.where(solo)[0]]
    groups = g.to_numpy()

    # Build the full matrix once, then subset columns per threshold.
    vocab = sorted({x for gs in gene_sets for x in gs})
    full = pd.DataFrame(
        {x: gene_sets.apply(lambda gs, x=x: int(x in gs)) for x in vocab}
    )
    counts = full.sum()

    prefix = RULE_PREFIX[drug]
    rule_cols_all = [c for c in vocab if c.startswith(prefix)]

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Drug:      {drug}")
    out(f"Cohort:    {len(df):,} isolates, {int(y.sum()):,} resistant")
    out(f"Vocabulary: {len(vocab)} distinct gene entries")
    out(f"Scheme:    {N_SPLITS}-fold x {N_REPEATS} repeats "
        f"= {N_SPLITS * N_REPEATS} folds per scheme")
    out()

    header = (f"{'min':>4}{'feats':>7}{'rule cols':>11}"
              f"{'rule rand':>11}{'log rand':>10}"
              f"{'rule grp':>10}{'log grp':>10}{'pass':>6}")
    out(header)
    out("-" * len(header))

    passenger_log = {}

    for t in THRESHOLDS:
        keep = [c for c in vocab if counts[c] >= t]
        Xall = full[keep].to_numpy()

        rule_cols = [c for c in keep if c.startswith(prefix)]
        if rule_cols:
            idx = [keep.index(c) for c in rule_cols]
            rule_pred = Xall[:, idx].max(axis=1)
        else:
            rule_pred = np.zeros(len(y), dtype=int)

        Xv = {"all": Xall, "rule": rule_pred}

        r_rand, l_rand = cv_errors(Xv, y, groups, "random", 0)
        r_grp, l_grp = cv_errors(Xv, y, groups, "grouped", 0)

        # Passengers: admitted, no mechanistic prefix, and every isolate
        # carrying the feature also satisfies the rule.
        passengers = []
        for c in keep:
            if c.startswith(MECHANISTIC):
                continue
            col = Xall[:, keep.index(c)].astype(bool)
            if col.sum() and rule_pred[col].mean() == 1.0:
                passengers.append((c, int(col.sum())))
        passenger_log[t] = passengers

        out(f"{t:>4}{len(keep):>7}{len(rule_cols):>11}"
            f"{1 - r_rand.mean():>11.4f}{1 - l_rand.mean():>10.4f}"
            f"{1 - r_grp.mean():>10.4f}{1 - l_grp.mean():>10.4f}"
            f"{len(passengers):>6}")

    out()
    out("rule cols = how many gyrA position-86 features survive the filter")
    out("pass      = admitted features with no mechanistic link that occur")
    out("            only in isolates already flagged by the rule")
    out()

    out("=" * 66)
    out("PASSENGERS BY THRESHOLD")
    out("=" * 66)
    for t in THRESHOLDS:
        ps = passenger_log[t]
        if not ps:
            out(f"min {t:>3}: none")
            continue
        listed = ", ".join(f"{c} (n={n})" for c, n in sorted(ps, key=lambda x: -x[1]))
        out(f"min {t:>3}: {listed}")

    out()
    out("A passenger is not a modelling error. It is a real correlation in")
    out("the data. The model cannot distinguish it from a cause using")
    out("presence and absence alone, which is why the mechanism has to")
    out("come from outside the dataset.")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"threshold_sweep_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()