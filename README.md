# campylobacter_amr_phenotype_prediction

Predicting antimicrobial resistance in *Campylobacter jejuni* from genomic features, using phenotypes measured on a plate rather than inferred from sequence.

Chapter 1 of **Machine Learning for Biology**, an open educational series.


## The question

Susceptibility testing means growing the organism and exposing it to the drug. *Campylobacter* is microaerophilic and slow, so that takes days, and the patient is treated on a guess in the meantime. Sequencing takes a day.

Can the genome answer the question the plate answers?


## What this repository finds

**A single mechanistically specified feature predicts ciprofloxacin resistance at 99.67%, and 84 features do not beat it.**

| Model | Features | Accuracy | Errors |
|---|---|---|---|
| Majority class | 0 | 0.7932 | 824 |
| Any substitution at gyrA 86 | 1 | 0.9967 | 13 |
| Logistic regression (25-fold grouped CV) | 58 | 0.9960 | n/a |

The cohort contains 3,984 isolates, 824 resistant and 3,160 susceptible (`results/metrics/cohort_summary_ciprofloxacin.txt`). The majority-class baseline is 824/3,984; the marker accuracy is exact on the full cohort, counting any of the four substitutions seen at gyrA 86, which are T86A, T86I, T86K and T86V (`results/metrics/marker_comparison_ciprofloxacin.txt`). The logistic regression figure is the mean accuracy over 25 grouped cross-validation folds at the default rarity threshold of 3, which admits 58 features (`results/metrics/cross_validation_ciprofloxacin.txt`).

The 84 features of the headline are the full vocabulary, obtained by setting the rarity threshold to 1 so that every distinct gene entry is retained. Fitted on those, logistic regression reaches 0.9960 under random cross-validation and 0.9957 under grouped, against 0.9967 for the rule (`results/metrics/threshold_sweep_ciprofloxacin.txt`).

Four results follow. Each is reproducible from the scripts here, and each number cites the file that produced it.

**The result is stable under grouped splitting.** A stratified random split placed 344 SNP clusters on both sides of the train/test boundary; the grouped split placed zero (`results/metrics/split_summary_ciprofloxacin.txt`). The rule's mean accuracy was 0.9965 under both random and grouped 25-fold cross-validation (`results/metrics/cross_validation_ciprofloxacin.txt`). That figure is 0.9965 rather than the 0.9967 above because the cross-validated rule is built from the three substitutions that survive the default rarity threshold, T86A, T86I and T86V. The single isolate carrying T86K falls below the threshold and is missed, which costs one error. A model exploiting lineage would have dropped. This one predicts resistance in genetic backgrounds it never saw.

**A beta-lactamase received the same weight as a real gyrase mutation.** In the grouped-split logistic regression, `gyrA_T86I` scored +8.845, `gyrA_T86V` scored +2.364, and `blaOXA-493` scored +2.364 (`results/metrics/model_results_ciprofloxacin.txt`, grouped split). `blaOXA-493` encodes an enzyme that cuts beta-lactam rings; ciprofloxacin has none. Every isolate carrying `blaOXA-493` also carries a gyrA 86 substitution: the co-occurrence table gives P(rule|feature) = 1.000 for it across all five (`results/metrics/cross_validation_ciprofloxacin.txt`). The gene never appears without the mutation, so nothing in the data measures its effect independently. (The `alone` column in `results/metrics/variant_table_ciprofloxacin.txt` reads 5 for this variant and answers a different question: whether an isolate carries any other variant from the same gene family. It says nothing about gyrA.)

**An arbitrary rarity threshold deleted the best-supported feature.** An early version filtered out gene entries seen in fewer than 10 isolates. `gyrA_T86V` appears in 8 isolates, across 7 distinct SNP clusters, and is resistant in all of them (`results/metrics/variant_table_ciprofloxacin.txt`). It was removed silently, and no accuracy figure revealed it. `scripts/10_threshold_sweep.py` sweeps the threshold from 1 to 20; the rule's cross-validated accuracy ranges from 0.9954 to 0.9975 across grouped and random schemes (`results/metrics/threshold_sweep_ciprofloxacin.txt`), and the conclusion does not depend on the choice.

**Annotation quality carries signal.** For tetracycline, any *tet* call gives 0.9932 accuracy; requiring an intact call, excluding `=PARTIAL` and `=MISTRANSLATION` but retaining `=PARTIAL_END_OF_CONTIG`, gives 0.9945 (`results/metrics/refined_features_tetracycline.txt`). Four isolates carrying only a truncated `tet(O)`, in four separate groups (three unclustered, one in its own cluster), are all susceptible. The one isolate carrying only an end-of-contig `tet(O)` call is resistant (`results/metrics/variant_table_tetracycline.txt`). A truncated gene and a gene split by an assembler share a name and do not share a phenotype.

Full outputs are in `results/metrics/`.


## Data

**Source.** NCBI Pathogen Detection, *Campylobacter* organism group, release **PDG000000003.2859** (`results/metrics/provenance_ciprofloxacin.txt`).

The release is pinned rather than tracking `latest_snps`, so the cohort is reconstructible. Two files are downloaded: isolate metadata carrying AST phenotypes and AMRFinderPlus genotypes, and SNP cluster assignments. See `data/README.md` for regeneration instructions.

**Cohort construction.** Isolates are retained when they have a measured phenotype for the drug in question and a linked genome assembly. SNP cluster assignments are joined on assembly accession. The filter sequence and column definitions are documented in `data/README.md`.

| Drug | Isolates | Resistant | Susceptible | Clusters |
|---|---|---|---|---|
| Ciprofloxacin | 3,984 | 824 | 3,160 | 1,128 |
| Tetracycline | 3,983 | 1,902 | 2,081 | 1,128 |

Source: `results/metrics/cohort_summary_ciprofloxacin.txt` and `results/metrics/cohort_summary_tetracycline.txt`.

464 isolates have no cluster assignment and are treated as singleton groups when splitting.

**Features** are a multi-hot encoding of the `AMR_genotypes` field: one binary column per distinct gene entry. **Labels** come from the `AST_phenotypes` field, which records laboratory susceptibility testing.

These are two different columns of the same table, and keeping them apart is the point. Deriving the label from AMRFinderPlus, as an earlier version of this analysis did, makes the target a function of the input and produces near-perfect accuracy that means nothing.


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

Every script requires the drug as its first argument. None of them has a default. Run one without an argument and it prints a usage line naming the drugs it accepts, then exits 2. Earlier versions did have defaults, and because they disagreed with each other a bare run of the pipeline built one drug's cohort and modelled another, without failing anywhere. Substitute `tetracycline` to run the second cohort.

Random seeds are fixed. Data downloads are idempotent and write through a `.part` file, so an interrupted download cannot masquerade as a complete one.

The pipeline was rerun end to end for both drugs on 28 July 2026, under Python 3.11.15, pandas 2.2.3, scikit-learn 1.5.2 and numpy 1.26.4. Every file in `results/metrics/` regenerated byte-identically except the two provenance records, which carry the run date and therefore change on every run.


## Layout

```
envs/environment.yml          conda specification
scripts/                      numbered, modular, runnable end to end
data/
  raw/                        NCBI downloads (gitignored, regenerable)
  interim/                    cohort tables (gitignored, regenerable)
  processed/                  feature matrix, labels, fold assignments (gitignored)
results/
  metrics/                    all reported numbers (in version control)
  figures/                    plots (gitignored, regenerable)
posts/                        the LinkedIn series, verbatim, with the
                              file behind every number in each post
notebooks/                    empty, reserved for later chapters
models/                       empty, reserved for later chapters (gitignored)
```

`results/metrics/` is deliberately tracked. Those files are small, and they are the record of what was claimed on which day.

Eight of those files are cited nowhere above. They are the tetracycline outputs of scripts 04, 05, 08, 09 and 10, the tetracycline split summary and provenance record, and the ciprofloxacin discordance listing. Both drugs were run through the same pipeline; the tetracycline results are reported here only as far as the annotation-quality finding, and the rest is left in place so the record is complete rather than selective.


## Limitations

**One species, two drugs.** Tetracycline and ciprofloxacin resistance in *C. jejuni* are both dominated by single well-characterised mechanisms. Nothing here generalises to organisms where resistance is genuinely multi-genic, and the finding that one feature suffices is a statement about this biology rather than about AMR prediction in general.

**The label is a threshold.** S and R are a cut across a continuous MIC, placed by expert committees and revised over time. Modulation that moves an isolate within a category, such as efflux via CmeABC, is invisible to a binary label. This chapter predicts which side of a line an isolate falls on, not how resistant it is.

**Features are bounded by the detection tool.** AMRFinderPlus reports what its reference database contains. Every gyrA point mutation observed here is at position 86, which reflects the database as much as the organism. Mechanisms absent from the catalogue are absent from the feature matrix.

**Three refinements are underpowered.** Excluding `gyrA_T86A` reclassifies 4 isolates. Requiring intact *tet* calls reclassifies 5. Treating `=PARTIAL_END_OF_CONTIG` as broken rather than intact reclassifies 8, two of which become wrong. Each is a handful out of roughly 4,000 (`results/metrics/refined_features_ciprofloxacin.txt`, `results/metrics/refined_features_tetracycline.txt`).

Both files report each difference as a ratio against fold-to-fold spread, where a value below about 1 cannot be told apart from noise. The `gyrA_T86A` exclusion scores 0.6 and the full-length *tet* definition 0.7, both below that line. Requiring intact *tet* calls scores 1.1, which sits at the boundary rather than beyond it. All three are mechanistically motivated and all three were suggested by reading this dataset. They are hypotheses requiring independent confirmation, not established improvements.

**`gyrA_T86A` rests on two independent observations.** Four isolates carry it and all are susceptible, but three sit in the same SNP cluster (PDS000023369.33) with near-identical genotypes (`results/metrics/variant_table_ciprofloxacin.txt`). The effective sample size is closer to two than to four.

**No external validation.** All figures come from one cohort. Cross-validation estimates variability within these data; it says nothing about another laboratory, another country, or another decade.

**Cohort composition is uncharacterised.** This is surveillance data, aggregated from submitting laboratories. Its distribution across geography, isolation source, and collection year has not been analysed, and any of those could bias the result.

**The environment is pinned to minor versions, not patch versions.** `envs/environment.yml` fixes pandas 2.2 and scikit-learn 1.5 but not the patch release, so a fresh install may resolve to different builds than the ones used here. The byte-identical rerun reported above is an observation from one machine on one day, not a guarantee that a reader will reproduce it exactly.


## Series

Part of Machine Learning for Biology. Chapter 1 covers supervised binary classification: features and labels, baselines, confusion matrices, train and test splits, leakage through population structure, confounding, and cross-validation.

The posts in `posts/` are reproduced verbatim, including anything later found to be imprecise, with a footer on each mapping every figure to the file that produced it.


## Licence and citation

MIT. See `LICENSE`.

If you use this repository, cite it together with the release it was built from, since the numbers belong to that snapshot and not to NCBI Pathogen Detection in general.

```
Hussain, Q. (2026). campylobacter_amr_phenotype_prediction: predicting
ciprofloxacin and tetracycline resistance in Campylobacter jejuni from
AMRFinderPlus genotypes and measured AST phenotypes. Chapter 1 of Machine
Learning for Biology. NCBI Pathogen Detection release PDG000000003.2859.
https://github.com/Qasim-Hussain-Code/campylobacter_amr_phenotype_prediction
```

The underlying data are from NCBI Pathogen Detection and carry their own terms. The MIT licence covers the code and the derived metrics files in this repository, not the source records.
