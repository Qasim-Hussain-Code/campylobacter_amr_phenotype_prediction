#!/usr/bin/env python3
"""Build the Chapter 1 cohort from NCBI Pathogen Detection metadata.

Filters Campylobacter isolates to those with a measured tetracycline
phenotype and a linked genome assembly, then attaches SNP cluster
assignments for use in grouped splitting.

Outputs data/interim/cohort.parquet
"""

from pathlib import Path
import sys
import pandas as pd

PDG = "PDG000000003.2859"
ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"

DRUG = "tetracycline"

USE_COLS = [
    "target_acc",
    "asm_acc",
    "biosample_acc",
    "collection_date",
    "geo_loc_name",
    "isolation_source",
    "epi_type",
    "AST_phenotypes",
    "AMR_genotypes",
]


def parse_ast(field, drug):
    """Extract one drug's S/I/R call from an AST_phenotypes string.

    The field looks like:
        "azithromycin=S,ciprofloxacin=S,nalidixic acid=S,tetracycline=R"

    Quotes are part of the stored value. Drug names may contain spaces.
    Returns None when the field is absent or the drug was not tested.
    """
    if not isinstance(field, str) or field in ("", "NULL"):
        return None
    for pair in field.strip('"').split(","):
        if "=" not in pair:
            continue
        name, _, call = pair.partition("=")
        if name.strip() == drug:
            return call.strip()
    return None


def main():
    meta_path = RAW / f"{PDG}.metadata.tsv"
    clust_path = RAW / f"{PDG}.cluster_list.tsv"

    for p in (meta_path, clust_path):
        if not p.exists():
            sys.exit(f"[fail] missing {p}. Run 01_fetch_ncbi_metadata.sh first.")

    print(f"[read] {meta_path.name}")
    meta = pd.read_csv(
        meta_path,
        sep="\t",
        usecols=USE_COLS,
        dtype=str,
        na_values=["NULL"],
        low_memory=False,
    )
    print(f"       {len(meta):,} isolates in release")

    meta[DRUG] = meta["AST_phenotypes"].apply(lambda s: parse_ast(s, DRUG))

    cohort = meta[meta[DRUG].notna()].copy()
    print(f"[filt] {len(cohort):,} with a {DRUG} result")
    print(cohort[DRUG].value_counts().to_string())

    cohort = cohort[cohort["asm_acc"].notna()].copy()
    print(f"[filt] {len(cohort):,} of those with an assembly")

    print(f"[read] {clust_path.name}")
    clusters = pd.read_csv(
        clust_path,
        sep="\t",
        usecols=["PDS_acc", "gencoll_acc"],
        dtype=str,
    ).rename(columns={"PDS_acc": "snp_cluster", "gencoll_acc": "asm_acc"})

    clusters = clusters.drop_duplicates(subset="asm_acc")

    before = len(cohort)
    cohort = cohort.merge(clusters, on="asm_acc", how="left", validate="one_to_one")
    assert len(cohort) == before, "merge changed row count"

    missing = cohort["snp_cluster"].isna().sum()
    print(f"[join] {missing:,} isolates have no SNP cluster")

    assert cohort["asm_acc"].is_unique, "duplicate assembly accessions"
    assert cohort[DRUG].isin(["S", "I", "R"]).all(), "unexpected phenotype code"

    INTERIM.mkdir(parents=True, exist_ok=True)
    out = INTERIM / "cohort.parquet"
    cohort.to_parquet(out, index=False)

    print(f"\n[done] wrote {out.relative_to(ROOT)}")
    print(f"       {len(cohort):,} isolates x {cohort.shape[1]} columns")
    print(f"       {cohort['snp_cluster'].nunique():,} distinct SNP clusters")


if __name__ == "__main__":
    main()
