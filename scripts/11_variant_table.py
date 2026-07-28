#!/usr/bin/env python3
"""Tabulate individual resistance-gene variants against measured phenotype.

Counts alone are misleading in surveillance data. Four isolates from one
outbreak are one observation repeated, not four. This reports how many
distinct SNP clusters carry each variant alongside the raw count, so the
effective sample size is visible rather than assumed.

Usage:
    python scripts/11_variant_table.py [drug] [prefix ...]

Reads  data/interim/cohort_<drug>.parquet
Writes results/metrics/variant_table_<drug>.txt
"""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
METRICS = ROOT / "results" / "metrics"

DEFAULT_PREFIXES = {
    "ciprofloxacin": ["gyrA_", "cmeB", "cmeR", "blaOXA-493"],
    "nalidixic acid": ["gyrA_", "cmeB", "cmeR"],
    "tetracycline": ["tet(", "cmeB", "cmeR"],
}

# Report every isolate for a variant when the resistant fraction falls
# below this, since those are the ones worth reading individually.
DETAIL_BELOW = 0.95


def entries(genotypes):
    if not isinstance(genotypes, str) or genotypes in ("", "NULL"):
        return []
    return [e.strip() for e in genotypes.strip('"').split(",") if e.strip()]


def main():
    if len(sys.argv) < 2:
        print(f"usage: 11_variant_table.py <drug>   drugs: "
              f"{', '.join(DEFAULT_PREFIXES)}", file=sys.stderr)
        sys.exit(2)
    drug = sys.argv[1]
    slug = drug.replace(" ", "_")
    prefixes = sys.argv[2:] or DEFAULT_PREFIXES.get(drug, [])

    path = INTERIM / f"cohort_{slug}.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 02_build_cohort.py {drug!r} first.")

    df = pd.read_parquet(path)
    df["entries"] = df["AMR_genotypes"].apply(entries)
    resistant = df[drug] == "R"

    # Unclustered isolates count as their own group. We cannot show they
    # are related to anything, so we do not assume that they are.
    clusters = df["snp_cluster"].copy()
    solo = clusters.isna()
    clusters[solo] = [f"__solo_{i}" for i in range(int(solo.sum()))]

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Drug:     {drug}")
    out(f"Cohort:   {len(df):,} isolates, {int(resistant.sum()):,} resistant "
        f"({resistant.mean():.1%})")
    out(f"Prefixes: {prefixes}")
    out()

    variants = sorted({
        e for es in df["entries"] for e in es
        if any(e.startswith(p) for p in prefixes)
    })

    if not variants:
        out("No matching variants found.")
        return

    header = (f"{'variant':<30}{'n':>5}{'R':>6}{'S':>5}"
              f"{'% R':>8}{'alone':>7}{'clust':>7}")
    out(header)
    out("-" * len(header))

    detail = []

    for v in variants:
        has = df["entries"].apply(lambda es, v=v: v in es)
        n = int(has.sum())
        n_r = int((has & resistant).sum())
        n_s = n - n_r
        frac = n_r / n

        fam = next(p for p in prefixes if v.startswith(p))
        alone = int(df.loc[has, "entries"].apply(
            lambda es, fam=fam: sum(1 for e in es if e.startswith(fam)) == 1
        ).sum())

        n_clust = clusters[has].nunique()

        out(f"{v:<30}{n:>5}{n_r:>6}{n_s:>5}"
            f"{frac:>7.1%}{alone:>7}{n_clust:>7}")

        if frac < DETAIL_BELOW:
            detail.append((v, has, n, n_clust))

    out()
    out("alone = isolates carrying this variant and no other from the same")
    out("        gene family. A variant never seen alone cannot be told")
    out("        apart from the variants it travels with.")
    out("clust = distinct SNP clusters carrying it. When clust is much")
    out("        smaller than n, the isolates are near-identical and the")
    out("        effective sample size is closer to clust than to n.")

    if detail:
        out()
        out("=" * 72)
        out(f"VARIANTS BELOW {DETAIL_BELOW:.0%} RESISTANT, ISOLATE BY ISOLATE")
        out("=" * 72)
        for v, has, n, n_clust in detail:
            out()
            out(f"{v}   n={n}, {n_clust} distinct cluster"
                f"{'s' if n_clust != 1 else ''}")
            sub = df[has]
            for acc, cl, geno, call in zip(
                sub["asm_acc"], sub["snp_cluster"].fillna("(none)"),
                sub["AMR_genotypes"], sub[drug]
            ):
                out(f"  [{call}] {acc}  {cl}")
                out(f"        {geno}")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"variant_table_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()