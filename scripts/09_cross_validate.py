#!/usr/bin/env python3
"""Cross-validate the rule and logistic regression under both schemes.

One train/test split gives one number. When two models differ by a handful
of isolates, that number cannot say whether the difference is real or an
accident of which rows landed in the test fold. Repeated cross-validation
gives a distribution instead of a point.

Also reports how often each heavily weighted feature co-occurs with the
mechanistic rule, which is how a passenger gene is identified.

Usage:
    python scripts/09_cross_validate.py [drug]

Reads  data/processed/features_<drug>.parquet
       data/processed/splits_<drug>.parquet
Writes results/metrics/cross_validation_<drug>.txt
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
METRICS = ROOT / "results" / "metrics"

N_SPLITS = 5
N_REPEATS = 5

RULE_PREFIX = {
    "ciprofloxacin": "gyrA_T86",
    "nalidixic acid": "gyrA_T86",
    "tetracycline": "tet(",
}


def main():
    if len(sys.argv) < 2:
        print(f"usage: 09_cross_validate.py <drug>   drugs: {', '.join(RULE_PREFIX)}",
              file=sys.stderr)
        sys.exit(2)
    drug = sys.argv[1]
    slug = drug.replace(" ", "_")

    fx = PROCESSED / f"features_{slug}.parquet"
    sp = PROCESSED / f"splits_{slug}.parquet"
    for p in (fx, sp):
        if not p.exists():
            sys.exit(f"[fail] missing {p}")

    X = pd.read_parquet(fx)
    S = pd.read_parquet(sp)
    assert (X["asm_acc"].values == S["asm_acc"].values).all(), "rows misaligned"

    names = [c for c in X.columns if c != "asm_acc"]
    Xv = X[names].to_numpy()
    y = S["resistant"].to_numpy()

    # Isolates with no cluster become their own singleton group. We cannot
    # show they are related to anything, so we do not assume they are.
    g = S["snp_cluster"].copy()
    solo = g.isna()
    g[solo] = [f"__solo_{i}" for i in np.where(solo)[0]]
    groups = g.to_numpy()

    prefix = RULE_PREFIX[drug]
    rule_idx = [i for i, c in enumerate(names) if c.startswith(prefix)]
    rule_all = (Xv[:, rule_idx].max(axis=1) if rule_idx
                else np.zeros(len(y), dtype=int))

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Drug:     {drug}")
    out(f"Features: {len(names)}")
    out(f"Rule:     {[names[i] for i in rule_idx]}")
    out(f"Scheme:   {N_SPLITS}-fold, repeated {N_REPEATS} times "
        f"({N_SPLITS * N_REPEATS} folds per scheme)")
    out()

    for scheme in ("random", "grouped"):
        rule_err, log_err = [], []

        for seed in range(N_REPEATS):
            if scheme == "random":
                cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                                     random_state=seed)
                folds = cv.split(Xv, y)
            else:
                cv = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True,
                                          random_state=seed)
                folds = cv.split(Xv, y, groups)

            for tr, te in folds:
                rule_err.append(int((rule_all[te] != y[te]).sum()) / len(te))

                clf = LogisticRegression(max_iter=2000, random_state=seed)
                clf.fit(Xv[tr], y[tr])
                log_err.append(int((clf.predict(Xv[te]) != y[te]).sum()) / len(te))

        rule_err = np.array(rule_err)
        log_err = np.array(log_err)

        out("=" * 66)
        out(f"{scheme.upper()} FOLDS")
        out("=" * 66)
        out(f"{'model':<12}{'mean acc':>10}{'sd':>9}{'worst':>10}{'best':>10}")
        out("-" * 51)
        for name, e in (("rule", rule_err), ("logistic", log_err)):
            out(f"{name:<12}{1 - e.mean():>10.4f}{e.std():>9.4f}"
                f"{1 - e.max():>10.4f}{1 - e.min():>10.4f}")
        out()

        diff = rule_err - log_err
        wins = int((diff > 0).sum())
        losses = int((diff < 0).sum())
        ties = int((diff == 0).sum())
        out(f"logistic better in {wins} folds, worse in {losses}, tied in {ties}")
        out(f"mean accuracy difference: {diff.mean():+.4f} "
            f"(sd {diff.std():.4f})")
        out()

    # Which features travel with the mechanism rather than causing it.
    out("=" * 66)
    out("FEATURE CO-OCCURRENCE WITH THE RULE")
    out("=" * 66)
    clf = LogisticRegression(max_iter=2000, random_state=0).fit(Xv, y)
    coef = pd.Series(clf.coef_[0], index=names)

    out(f"{'feature':<28}{'n':>6}{'coef':>9}"
        f"{'P(rule|feat)':>14}{'P(feat|rule)':>14}")
    out("-" * 71)
    for feat in coef.sort_values(key=abs, ascending=False).head(8).index:
        col = Xv[:, names.index(feat)].astype(bool)
        n = int(col.sum())
        p_rule_given = rule_all[col].mean() if n else np.nan
        p_feat_given = col[rule_all == 1].mean() if rule_all.sum() else np.nan
        mark = "  <- rule" if names.index(feat) in rule_idx else ""
        out(f"{feat:<28}{n:>6}{coef[feat]:>9.2f}"
            f"{p_rule_given:>14.3f}{p_feat_given:>14.3f}{mark}")

    out()
    out("A feature with no mechanistic link to this drug but a high")
    out("P(rule|feature) is a passenger: it travels with the causal")
    out("mutation on a shared lineage or element, and the model cannot")
    out("separate the two from presence and absence alone.")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"cross_validation_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
