#!/usr/bin/env python3
"""Tabulate individual resistance-gene variants against measured phenotype.

The threshold sweep implied that not every substitution at gyrA 86 behaves
the same way. That was inferred from differences in accuracy. This checks
it directly, which is the only honest way to make the claim.

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


def entries(genotypes):
    if not isinstance(genotypes, str) or genotypes in ("", "NULL"):
        return []
    return [e.strip() for e in genotypes.strip('"').split(",") if e.strip()]


def main():
    drug = sys.argv[1] if len(sys.argv) > 1 else "ciprofloxacin"
    slug = drug.replace(" ", "_")
    prefixes = sys.argv[2:] or DEFAULT_PREFIXES.get(drug, [])

    path = INTERIM / f"cohort_{slug}.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 02_build_cohort.py {drug!r} first.")

    df = pd.read_parquet(path)
    df["entries"] = df["AMR_genotypes"].apply(entries)
    resistant = df[drug] == "R"

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
    else:
        header = f"{'variant':<28}{'n':>6}{'R':>6}{'S':>6}{'% R':>8}{'alone':>7}"
        out(header)
        out("-" * len(header))

        for v in variants:
            has = df["entries"].apply(lambda es, v=v: v in es)
            n = int(has.sum())
            n_r = int((has & resistant).sum())
            n_s = n - n_r

            # How many carry this variant and no other from the same prefix
            # family. A variant that never appears alone cannot be assessed
            # independently of the ones it travels with.
            fam = [p for p in prefixes if v.startswith(p)][0]
            alone = int(df.loc[has, "entries"].apply(
                lambda es, v=v, fam=fam:
                sum(1 for e in es if e.startswith(fam)) == 1
            ).sum())

            out(f"{v:<28}{n:>6}{n_r:>6}{n_s:>6}"
                f"{n_r / n:>7.0%}{alone:>7}")

        out()
        out("alone = isolates carrying this variant and no other from the")
        out("        same gene family. A variant never seen alone cannot")
        out("        be told apart from its companions.")

    out()
    out("=" * 60)
    out("ISOLATES CARRYING A VARIANT THAT IS MOSTLY SUSCEPTIBLE")
    out("=" * 60)
    for v in variants:
        has = df["entries"].apply(lambda es, v=v: v in es)
        n = int(has.sum())
        if n == 0 or (has & resistant).sum() / n >= 0.5:
            continue
        out(f"\n{v}  ({n} isolates)")
        sub = df[has]
        for acc, geno, call in zip(sub["asm_acc"], sub["AMR_genotypes"],
                                   sub[drug]):
            out(f"  [{call}] {acc}  {geno}")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"variant_table_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()