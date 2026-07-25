#!/usr/bin/env python3
"""Test feature definitions that account for annotation quality.

AMRFinderPlus qualifies some calls: =PARTIAL means only part of the
reference matched, =PARTIAL_END_OF_CONTIG means the match ran off the edge
of a contig, =MISTRANSLATION flags an internal stop or frameshift.

A prefix match on the gene name discards all of that. A truncated gene may
produce no functional protein, whereas a gene split by an assembler is
intact. These are different biological situations with the same gene name.

POST HOC WARNING. The intact-versus-fragment distinction was noticed while
reading this dataset. The mechanism justifies it independently, but the
improvement reported here is measured on the data that suggested it and is
therefore optimistic. Independent confirmation is required before claiming
the gain.

Usage:
    python scripts/12_test_refined_features.py [drug]

Reads  data/interim/cohort_<drug>.parquet
Writes results/metrics/refined_features_<drug>.txt
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
METRICS = ROOT / "results" / "metrics"

N_SPLITS = 5
N_REPEATS = 5

# Qualifiers that leave the coding sequence intact. An absent qualifier
# means a full-length match. END_OF_CONTIG is an assembly boundary, not a
# property of the organism.
INTACT_QUALS = {"", "PARTIAL_END_OF_CONTIG"}


def parse(entry):
    """Split 'tet(O)=PARTIAL' into ('tet(O)', 'PARTIAL')."""
    gene, _, qual = entry.partition("=")
    return gene.strip(), qual.strip()


def entries(genotypes):
    if not isinstance(genotypes, str) or genotypes in ("", "NULL"):
        return []
    return [parse(e) for e in genotypes.strip('"').split(",") if e.strip()]


def make_defs(drug):
    """Named predicates over the parsed (gene, qualifier) list."""
    if drug in ("ciprofloxacin", "nalidixic acid"):
        return [
            ("any gyrA 86 change",
             lambda es: any(g.startswith("gyrA_T86") for g, q in es)),
            ("gyrA 86, intact calls",
             lambda es: any(g.startswith("gyrA_T86") and q in INTACT_QUALS
                            for g, q in es)),
            ("gyrA 86, excluding T86A",
             lambda es: any(g.startswith("gyrA_T86") and not g.startswith("gyrA_T86A")
                            for g, q in es)),
        ]
    if drug == "tetracycline":
        return [
            ("any tet call",
             lambda es: any(g.startswith("tet(") for g, q in es)),
            ("intact tet call",
             lambda es: any(g.startswith("tet(") and q in INTACT_QUALS
                            for g, q in es)),
            ("unqualified tet call only",
             lambda es: any(g.startswith("tet(") and q == "" for g, q in es)),
        ]
    sys.exit(f"[fail] no definitions for {drug!r}")


def cells(pred, y):
    tp = int((pred & (y == 1)).sum())
    fn = int((~pred & (y == 1)).sum())
    fp = int((pred & (y == 0)).sum())
    tn = int((~pred & (y == 0)).sum())
    return tp, fn, fp, tn


def main():
    drug = sys.argv[1] if len(sys.argv) > 1 else "tetracycline"
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

    defs = make_defs(drug)
    preds = {name: parsed.apply(fn).to_numpy() for name, fn in defs}

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Drug:   {drug}")
    out(f"Cohort: {len(df):,} isolates, {int(y.sum()):,} resistant")
    out(f"Folds:  {N_SPLITS} x {N_REPEATS} repeats, grouped by SNP cluster")
    out()
    out("POST HOC. These definitions were shaped by reading this dataset.")
    out("The mechanism supports them independently, but the numbers below")
    out("are measured on the data that suggested them.")
    out()

    header = (f"{'definition':<28}{'TP':>6}{'FN':>5}{'FP':>5}{'TN':>6}"
              f"{'acc':>9}{'cv mean':>10}{'cv sd':>9}")
    out(header)
    out("-" * len(header))

    fold_acc = {}

    for name in preds:
        p = preds[name]
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
            f"{acc:>9.4f}{accs.mean():>10.4f}{accs.std():>9.4f}")

    out()
    out("A rule has no fitted parameters, so its cross-validated mean equals")
    out("its accuracy on the whole cohort. The spread is what matters: it")
    out("says whether a difference between two rules is larger than the")
    out("variation between folds.")
    out()

    names = list(preds)
    base = names[0]
    out(f"{'comparison':<44}{'mean diff':>11}{'sd':>9}{'ratio':>8}")
    out("-" * 72)
    for other in names[1:]:
        d = fold_acc[other] - fold_acc[base]
        ratio = abs(d.mean()) / d.std() if d.std() else np.inf
        out(f"{other + ' vs ' + base:<44}{d.mean():>+11.4f}"
            f"{d.std():>9.4f}{ratio:>8.1f}")

    out()
    out("ratio = size of the difference relative to fold-to-fold spread.")
    out("Below about 1 the difference is indistinguishable from noise.")
    out()

    out("=" * 72)
    out("WHERE THE DEFINITIONS DISAGREE")
    out("=" * 72)
    for other in names[1:]:
        differ = preds[base] != preds[other]
        n = int(differ.sum())
        if not n:
            out(f"\n{other}: identical to '{base}'")
            continue
        sub = df[differ]
        correct = int((preds[other][differ] == y[differ]).sum())
        out(f"\n{other}: {n} isolates change, "
            f"{correct} of them now correct")
        for acc_, geno, call in zip(sub["asm_acc"], sub["AMR_genotypes"],
                                    sub[drug]):
            out(f"  [{call}] {acc_}  {geno}")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"refined_features_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()