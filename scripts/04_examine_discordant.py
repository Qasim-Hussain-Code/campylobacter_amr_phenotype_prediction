#!/usr/bin/env python3
"""Inspect isolates where the single-marker rule disagrees with the plate.

These are the only isolates carrying information the marker cannot explain,
so they set the ceiling on what any model of this feature can achieve.

Usage:
    python scripts/04_examine_discordant.py [drug]

Reads  data/interim/cohort_<drug>.parquet
Writes results/metrics/discordant_<drug>.txt
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

# Gene families worth counting separately when they appear in a
# discordant isolate. Prefix matched, same convention as MARKERS.
RELATED = {
    "tetracycline": ["tet("],
    "ciprofloxacin": ["gyrA_", "parC_", "cmeR_", "cmeB"],
    "nalidixic acid": ["gyrA_", "parC_"],
    "erythromycin": ["23S_", "erm("],
    "azithromycin": ["23S_", "erm("],
}


def entries(genotypes):
    """Split an AMR_genotypes field into its individual gene entries."""
    if not isinstance(genotypes, str) or genotypes in ("", "NULL"):
        return []
    return [e.strip() for e in genotypes.strip('"').split(",") if e.strip()]


def main():
    if len(sys.argv) < 2:
        print(f"usage: 04_examine_discordant.py <drug>   drugs: {', '.join(MARKERS)}",
              file=sys.stderr)
        sys.exit(2)
    drug = sys.argv[1]
    slug = drug.replace(" ", "_")

    if drug not in MARKERS:
        sys.exit(f"[fail] no marker defined for '{drug}'. Known: {list(MARKERS)}")
    marker = MARKERS[drug]

    path = INTERIM / f"cohort_{slug}.parquet"
    if not path.exists():
        sys.exit(f"[fail] missing {path}. Run 02_build_cohort.py {drug!r} first.")

    df = pd.read_parquet(path)
    df["entries"] = df["AMR_genotypes"].apply(entries)
    df["marker"] = df["entries"].apply(
        lambda es: any(e.startswith(marker) for e in es)
    )

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"Drug:   {drug}")
    out(f"Marker: {marker}")
    out(f"Cohort: {len(df):,} isolates")
    out()

    fn = df[(df[drug] == "R") & (~df["marker"])]
    fp = df[(df[drug] == "S") & (df["marker"])]

    out("=" * 60)
    out(f"RESISTANT on the plate, no {marker} detected: {len(fn):,}")
    out("=" * 60)
    out("Some other mechanism, a marker this rule does not cover,")
    out("an assembly that missed the gene, or a plate error.")
    out()

    if len(fn):
        related = RELATED.get(drug, [])
        hits = {}
        for es in fn["entries"]:
            for e in es:
                if any(e.startswith(p) for p in related):
                    hits[e] = hits.get(e, 0) + 1
        if hits:
            out(f"Related genes present in these isolates:")
            for gene, n in sorted(hits.items(), key=lambda kv: -kv[1]):
                out(f"  {n:4d}  {gene}")
        else:
            out("No related resistance gene found in any of them.")
        out()

        n_empty = fn["entries"].apply(len).eq(0).sum()
        out(f"Isolates with no AMR genes called at all: {n_empty:,}")
        out()
        out("Full genotype strings:")
        for acc, g in zip(fn["asm_acc"], fn["AMR_genotypes"].fillna("(none)")):
            out(f"  {acc}  {g}")
        out()

    out("=" * 60)
    out(f"SUSCEPTIBLE on the plate, {marker} detected: {len(fp):,}")
    out("=" * 60)
    out("A disrupted or non-expressed gene, an MIC sitting just below")
    out("the breakpoint, or a plate error.")
    out()

    if len(fp):
        out("Full genotype strings:")
        for acc, g in zip(fp["asm_acc"], fp["AMR_genotypes"].fillna("(none)")):
            out(f"  {acc}  {g}")
        out()

    total_err = len(fn) + len(fp)
    out(f"Total disagreements: {total_err:,} of {len(df):,} "
        f"({total_err / len(df):.1%})")
    out(f"Ceiling for any model using this feature alone: "
        f"{1 - total_err / len(df):.3f}")

    METRICS.mkdir(parents=True, exist_ok=True)
    dest = METRICS / f"discordant_{slug}.txt"
    dest.write_text("\n".join(lines) + "\n")
    print(f"\n[done] wrote {dest.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
