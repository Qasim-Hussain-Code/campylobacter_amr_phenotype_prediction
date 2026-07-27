# campylobacter_amr_phenotype_prediction

Predicting antimicrobial resistance in *Campylobacter jejuni* from genomic
features, using phenotypes measured on a plate rather than inferred from
sequence.

Chapter 1 of **Machine Learning for Biology**, an open educational series.

---

## The question

Susceptibility testing means growing the organism and exposing it to the
drug. *Campylobacter* is microaerophilic and slow, so that takes days, and
the patient is treated on a guess in the meantime. Sequencing takes a day.

Can the genome answer the question the plate answers?

---

## What this repository finds

**A single mechanistically specified feature predicts ciprofloxacin
resistance at 99.67%, and 84 features do not beat it.**

| model | features | accuracy | errors |
|---|---|---|---|
| majority class | 0 | 0.7932 | 824 |
| any substitution at gyrA 86 | 1 | 0.9967 | 13 |
| logistic regression | 58 | 0.9961 | – |

Cohort: 3,984 isolates, 824 resistant. Accuracies for the rule are exact on
the full cohort; the logistic figure is the mean over 25 grouped
cross-validation folds.

Four further results, each reproducible from the scripts here.

**The result is stable under grouped splitting.** 344 SNP clusters had
isolates on both sides of a random split and none after grouping by
cluster. Accuracy did not move: 0.9965 either way. A model exploiting
lineage would have dropped. This one predicts resistance in genetic
backgrounds it never saw.

**A beta-lactamase received the same weight as a real gyrase mutation.**
In the grouped fit, `gyrA_T86I` scored 8.845, `gyrA_T86V` scored 2.364,
and `blaOXA-493` scored 2.364. `blaOXA-493` encodes an enzyme that cuts
beta-lactam rings; ciprofloxacin has none. All five isolates carrying it
also carry a gyrA 86 substitution, and no isolate in the cohort carries it
alone, so nothing in the data measures its effect independently.

**An arbitrary rarity threshold deleted the best-supported feature.** An
early version filtered out gene entries seen in fewer than 10 isolates.
`gyrA_T86V` appears in 8, across 7 distinct SNP clusters, and is resistant
in all of them. It was removed silently, and no accuracy figure revealed
it. `scripts/10_threshold_sweep.py` now sweeps the threshold from 1 to 20;
accuracy ranges from 0.9954 to 0.9975 and the conclusion does not depend
on the choice.

**Annotation quality carries signal.** For tetracycline, any *tet* call
gives 0.9932. Requiring an intact call, excluding `=PARTIAL` and
`=MISTRANSLATION` but retaining `=PARTIAL_END_OF_CONTIG`, gives 0.9945.
Four isolates carrying only a truncated `tet(O)`, in four separate
clusters, are all susceptible; the one isolate carrying only an
end-of-contig call is resistant. A truncated gene and a gene split by an
assembler share a name and do not share a phenotype.

Full outputs are in `results/metrics/`.

---

## Data

**Source:** NCBI Pathogen Detection, *Campylobacter* organism group,
release **PDG000000003.2859**.

The release is pinned rather than tracking `latest_snps`, so the cohort is
reconstructible. Two files are downloaded: isolate metadata carrying
AST phenotypes and AMRFinderPlus genotypes, and SNP cluster assignments.

**Cohort construction.** Isolates are retained when they have a measured
phenotype for the drug in question and a linked genome assembly. SNP
cluster assignments are joined on assembly accession, never by position.

| drug | isolates | resistant | susceptible | clusters |
|---|---|---|---|---|
| ciprofloxacin | 3,984 | 824 | 3,160 | 1,128 |
| tetracycline | 3,983 | 1,902 | 2,081 | 1,128 |

464 isolates have no cluster assignment and are treated as singleton
groups when splitting.

**Features** are a multi-hot encoding of the `AMR_genotypes` field: one
binary column per distinct gene entry. **Labels** come from the
`AST_phenotypes` field, which records laboratory susceptibility testing.

These are two different columns of the same table, and keeping them apart
is the point. Deriving the label from AMRFinderPlus, as an earlier version
of this analysis did, makes the target a function of the input and
produces near-perfect accuracy that means nothing.

See `data/README.md` for regeneration instructions.

---

## Reproducing

```bash
bash scripts/00_setup_environment.sh
conda activate campy-amr

bash scripts/01_fetch_ncbi_metadata.sh

python scripts/02_build_cohort.py ciprofloxacin
python scripts/03_explore_cohort.py ciprofloxacin
python scripts/04_examine_discordant.py ciprofloxacin
python scripts/05_compare_markers.py ciprofloxacin
python scripts/06_build_features.py ciprofloxacin
python scripts/07_split_data.py ciprofloxacin
python scripts/08_train_models.py ciprofloxacin
python scripts/09_cross_validate.py ciprofloxacin
python scripts/10_threshold_sweep.py ciprofloxacin
python scripts/11_variant_table.py ciprofloxacin
python scripts/12_test_refined_features.py ciprofloxacin
```

Every script takes the drug as an argument and defaults sensibly. Substitute
`tetracycline` to run the second cohort. Random seeds are fixed. Data
downloads are idempotent and write through a `.part` file, so an
interrupted download cannot masquerade as a complete one.

---

## Layout

```
envs/environment.yml          conda specification
scripts/                      numbered, modular, runnable end to end
data/
  raw/                        NCBI downloads (gitignored, regenerable)
  interim/                    cohort tables
  processed/                  feature matrix, labels, fold assignments
results/
  metrics/                    all reported numbers (in version control)
  figures/                    plots (gitignored, regenerable)
posts/                        the LinkedIn series, verbatim, with the
                              file behind every number in each post
```

`results/metrics/` is deliberately tracked. Those files are small, and they
are the record of what was claimed on which day.

---

## Limitations

**One species, two drugs.** Tetracycline and ciprofloxacin resistance in
*C. jejuni* are both dominated by single well-characterised mechanisms.
Nothing here generalises to organisms where resistance is genuinely
multi-genic, and the finding that one feature suffices is a statement
about this biology rather than about AMR prediction in general.

**The label is a threshold.** S and R are a cut across a continuous MIC,
placed by expert committees and revised over time. Modulation that moves
an isolate within a category, such as efflux via CmeABC, is invisible to a
binary label. This chapter predicts which side of a line an isolate falls
on, not how resistant it is.

**Features are bounded by the detection tool.** AMRFinderPlus reports what
its reference database contains. Every gyrA point mutation observed here
is at position 86, which reflects the database as much as the organism.
Mechanisms absent from the catalogue are absent from the feature matrix.

**Three refinements are underpowered.** Excluding `gyrA_T86A`, requiring
intact *tet* calls, and treating `=PARTIAL_END_OF_CONTIG` as intact each
move 2 to 5 isolates out of roughly 4,000. All three are mechanistically
motivated, all three were suggested by reading this dataset, and the
differences are smaller than the fold-to-fold variation. They are
hypotheses requiring independent confirmation, not established
improvements.

**`gyrA_T86A` rests on two independent observations.** Four isolates carry
it and all are susceptible, but three sit in the same SNP cluster with
near-identical genotypes. The effective sample size is closer to two than
to four.

**No external validation.** All figures come from one cohort. Cross
validation estimates variability within these data; it says nothing about
another laboratory, another country, or another decade.

**Cohort composition is uncharacterised.** This is surveillance data,
aggregated from submitting laboratories. Its distribution across
geography, isolation source and collection year has not been analysed, and
any of those could bias the result.

---

## Series

Part of Machine Learning for Biology. Chapter 1 covers supervised binary
classification: features and labels, baselines, confusion matrices,
train and test splits, leakage through population structure, confounding,
and cross validation.

The posts in `posts/` are reproduced verbatim, including anything later
found to be imprecise, with a footer on each mapping every figure to the
file that produced it.

---

## Licence

MIT. See `LICENSE`.