#!/usr/bin/env python3
"""Test feature definitions that account for annotation quality.

AMRFinderPlus attaches two different kinds of suffix to a gene call, and
they must not be confused:

  method flags   POINT, HMM        how the call was made
  quality flags  PARTIAL           only part of the reference matched
                 MISTRANSLATION    internal stop or frameshift detected
                 PARTIAL_END_OF_CONTIG
                                   the match ran off a contig edge, which is
                                   an assembly boundary rather than a
                                   property of the organism

A prefix match on the gene name discards all of it. A truncated gene may
produce no functional protein; a gene split by an assembler is intact.

POST HOC WARNING. The intact-versus-fragment distinction was noticed while
reading this dataset. The mechanism justifies it independently, but any
improvement reported here is measured on the data that suggested it and is
therefore optimistic. Independent confirmation is required.

Usage:
    python scripts/12_test_refined_features.py [drug]

Reads  data/interim/cohort_<drug>.parquet
Writes results/metrics/refined_features_<drug>.txt
"""

from pathlib import Path
from collections import Counter
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
METRICS = ROOT / "results" / "metrics"

N_SPLITS = 5
N_REPEATS = 5
MAX_LISTED = 25

# Qualifiers indicating the coding sequence is compromised. Defined by
# exclusion so that an unfamiliar qualifier defaults to a valid call
# rather than silently removing every isolate that carries it.
BROKEN_QUALS = {"PARTIAL", "MISTRANSLATION"}


def parse(entry):
    """Split 'tet(O)=PARTIAL' into ('tet(O)', 'PARTIAL')."""
    gene, _, qual = entry.partition("=")
    return gene.strip(), qual.strip()


def entries(genotypes):
    if not isinstance(genotypes, str) or genotypes in ("", "NULL"):
        return []
    return [parse(e) for e in genotypes.strip('"').split(",") if e.strip()]


def intact(qual):
    return qual not in BROKEN_QUALS


def make_defs(drug):
    if drug in ("ciprofloxacin", "nalidixic acid"):
        return [
            ("any gyrA 86 change",
             lambda es: any(g.startswith("gyrA_T86") for g, q in es)),
            ("gyrA 86, intact calls",
             lambda es: any(g.startswith("gyrA_T86") and intact(q)
                            for g, q in es)),
            ("gyrA 86, excluding T86A",
             lambda es: any(g.startswith("gyrA_T86")
                            and not g.startswith("gyrA_T86A")
                            for g, q in es)),
            ("gyrA 86, intact, no T86A",
             lambda es: any(g.startswith("gyrA_T86")
                            and not g.startswith("gyrA_T86A")
                            and intact(q) for g, q in es)),
        ]
    if drug == "tetracycline":
        return [
            ("any tet call",
             lambda es: any(g.startswith("tet(") for g, q in es)),
            ("intact tet call",
             lambda es: any(g.startswith("tet(") and intact(q) for g, q in es)),
            ("full-length tet only",
             lambda es: any(g.startswith("tet(") and q == "" for g, q in es)),
        ]
    sys.exit(f"[fail] no definitions for {drug!r}")


def cells(pred, y):
    return (int((pred & (y == 1)).sum()), int((~pred & (y == 1)).sum()),
            int((pred & (y == 0)).sum()), int((~pred & (y == 0)).sum()))


def main():
    if len(sys.argv) < 2:
        print("usage: 12_test_refined_features.py <drug>   drugs: ciprofloxacin, "
              "nalidixic acid, tetracycline", file=sys.stderr)
        sys.exit(2)
    drug = sys.argv[1]
    slug = drug.replace(" ", "_")

    path = INTERIM / f"cohort_{slug}.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 02_build_cohort.py {drug!r} first.")

    df = pd.read_parquet(path)
    parsed = df["AMR_genotypes"].apply(entries)
    y = (df[drug] == "R").to_numpy().astype(int)

    g = df["snp_cluster"].copy()
    solo = g.isna()
    g[solo] = [f"__solo_{i}" for i in range(int(solo.sum()))]
    groups = g.to_numpy()

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Drug:   {drug}")
    out(f"Cohort: {len(df):,} isolates, {int(y.sum()):,} resistant")
    out(f"Folds:  {N_SPLITS} x {N_REPEATS} repeats, grouped by SNP cluster")
    out()

    # Show the qualifier vocabulary before applying any definition, so a
    # misclassified qualifier is visible rather than silent.
    quals = Counter(q for es in parsed for _, q in es)
    out("qualifiers present in AMR_genotypes:")
    for q, n in quals.most_common():
        tag = "broken" if q in BROKEN_QUALS else "valid"
        out(f"  {(q or '(none)'):<24}{n:>7}   treated as {tag}")
    out()

    out("POST HOC. These definitions were shaped by reading this dataset.")
    out("The mechanism supports them independently, but the numbers below")
    out("are measured on the data that suggested them.")
    out()

    defs = make_defs(drug)
    preds = {name: parsed.apply(fn).to_numpy() for name, fn in defs}

    header = (f"{'definition':<28}{'TP':>6}{'FN':>5}{'FP':>5}{'TN':>6}"
              f"{'acc':>9}{'cv sd':>9}")
    out(header)
    out("-" * len(header))

    fold_acc = {}

    for name, p in preds.items():
        tp, fn, fp, tn = cells(p, y)
        acc = (tp + tn) / len(y)

        accs = []
        for seed in range(N_REPEATS):
            cv = StratifiedGroupKFold(N_SPLITS, shuffle=True, random_state=seed)
            for _, te in cv.split(np.zeros(len(y)), y, groups):
                accs.append((p[te] == y[te]).mean())
        accs = np.array(accs)
        fold_acc[name] = accs

        out(f"{name:<28}{tp:>6}{fn:>5}{fp:>5}{tn:>6}"
            f"{acc:>9.4f}{accs.std():>9.4f}")

    out()
    out("A rule has no fitted parameters, so its cross-validated mean equals")
    out("its accuracy on the whole cohort. Only the spread is new: it says")
    out("whether a gap between two rules exceeds fold-to-fold variation.")
    out()

    names = list(preds)
    base = names[0]
    out(f"{'comparison':<46}{'mean diff':>11}{'sd':>9}{'ratio':>8}")
    out("-" * 74)
    for other in names[1:]:
        d = fold_acc[other] - fold_acc[base]
        ratio = abs(d.mean()) / d.std() if d.std() else float("inf")
        out(f"{(other + ' vs ' + base):<46}{d.mean():>+11.4f}"
            f"{d.std():>9.4f}{ratio:>8.1f}")

    out()
    out("ratio = difference relative to fold-to-fold spread. Below about 1")
    out("        the difference is indistinguishable from noise.")
    out()

    out("=" * 74)
    out("WHERE THE DEFINITIONS DISAGREE")
    out("=" * 74)
    for other in names[1:]:
        differ = preds[base] != preds[other]
        n = int(differ.sum())
        if not n:
            out(f"\n{other}: identical to '{base}'")
            continue
        correct = int((preds[other][differ] == y[differ]).sum())
        out(f"\n{other}: {n} isolates change, {correct} now correct, "
            f"{n - correct} now wrong")
        sub = df[differ].head(MAX_LISTED)
        for acc_, geno, call in zip(sub["asm_acc"], sub["AMR_genotypes"],
                                    sub[drug]):
            out(f"  [{call}] {acc_}  {geno}")
        if n > MAX_LISTED:
            out(f"  ... {n - MAX_LISTED} more not shown")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"refined_features_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()