#!/usr/bin/env python3
"""Train and evaluate three models on both splits.

  baseline   always predict the majority class. No features, no parameters.
  rule       gyrA_T86 present. One feature, no fitted parameters.
  logistic   logistic regression on all features. 40 fitted parameters.

Each is scored on the training rows it saw and the test rows it did not.
The gap between those two numbers is the overfitting signal.

Usage:
    python scripts/08_train_models.py [drug]

Reads  data/processed/features_<drug>.parquet
       data/processed/splits_<drug>.parquet
Writes results/metrics/model_results_<drug>.txt
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
METRICS = ROOT / "results" / "metrics"

SEED = 42

RULE_PREFIX = {
    "ciprofloxacin": "gyrA_T86",
    "nalidixic acid": "gyrA_T86",
    "tetracycline": "tet(",
}


def scores(y_true, y_pred):
    """Return the four confusion cells plus the usual rates."""
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        "tp": int(tp), "fn": int(fn), "fp": int(fp), "tn": int(tn),
        "sens": tp / (tp + fn) if tp + fn else np.nan,
        "spec": tn / (tn + fp) if tn + fp else np.nan,
        "acc": (tp + tn) / len(y_true),
        "err": int(fp + fn),
    }


def main():
    drug = sys.argv[1] if len(sys.argv) > 1 else "ciprofloxacin"
    slug = drug.replace(" ", "_")

    fx = PROCESSED / f"features_{slug}.parquet"
    sp = PROCESSED / f"splits_{slug}.parquet"
    for p in (fx, sp):
        if not p.exists():
            sys.exit(f"[fail] missing {p}")

    X = pd.read_parquet(fx)
    S = pd.read_parquet(sp)

    assert (X["asm_acc"].values == S["asm_acc"].values).all(), "rows misaligned"

    feature_names = [c for c in X.columns if c != "asm_acc"]
    Xv = X[feature_names].to_numpy()
    y = S["resistant"].to_numpy()

    prefix = RULE_PREFIX[drug]
    rule_cols = [c for c in feature_names if c.startswith(prefix)]

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Drug:     {drug}")
    out(f"Features: {len(feature_names)}")
    out(f"Rule:     any of {rule_cols}")
    out()

    for split_col in ("split_random", "split_grouped"):
        kind = split_col.split("_")[1]
        tr = (S[split_col] == "train").to_numpy()
        te = ~tr

        out("=" * 74)
        out(f"{kind.upper()} SPLIT   train {tr.sum():,}   test {te.sum():,}")
        out(f"test set is {y[te].mean():.1%} resistant, "
            f"so the majority-class floor on test is "
            f"{max(y[te].mean(), 1 - y[te].mean()):.4f}")
        out("=" * 74)

        preds = {}

        # Baseline. Learns one thing from the training set: which class
        # is commoner. Then says that, always.
        majority = int(round(y[tr].mean()))
        preds["baseline"] = np.full(len(y), majority)

        # Rule. Nothing is fitted. The feature was chosen from mechanism
        # before any accuracy was computed.
        preds["rule"] = (Xv[:, [feature_names.index(c) for c in rule_cols]]
                         .max(axis=1) if rule_cols else np.zeros(len(y), int))

        # Logistic regression. Fits one weight per feature plus an
        # intercept, using only the training rows.
        clf = LogisticRegression(max_iter=2000, random_state=SEED)
        clf.fit(Xv[tr], y[tr])
        preds["logistic"] = clf.predict(Xv)

        header = (f"{'model':<12}{'set':<7}{'TP':>6}{'FN':>5}{'FP':>5}{'TN':>6}"
                  f"{'sens':>8}{'spec':>8}{'acc':>9}{'err':>6}")
        out(header)
        out("-" * len(header))

        for name, p in preds.items():
            for label, mask in (("train", tr), ("test", te)):
                s = scores(y[mask], p[mask])
                out(f"{name:<12}{label:<7}{s['tp']:>6}{s['fn']:>5}{s['fp']:>5}"
                    f"{s['tn']:>6}{s['sens']:>8.3f}{s['spec']:>8.3f}"
                    f"{s['acc']:>9.4f}{s['err']:>6}")
            out()

        # What the fitted model actually learned.
        coef = pd.Series(clf.coef_[0], index=feature_names).sort_values(
            key=abs, ascending=False)
        out("largest logistic coefficients:")
        for gene, w in coef.head(6).items():
            out(f"  {w:+8.3f}  {gene}")
        out()

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"model_results_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
