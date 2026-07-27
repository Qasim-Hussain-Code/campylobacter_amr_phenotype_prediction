# Day 5: Data Leakage and Grouped Splitting

Chapter 1: Predicting antibiotic resistance from a genome

Yesterday I said there was a way to catch a model reading correlation in place of causation.

Hide whole families of isolates from it.

To see why that helps, start with how a model gets tested at all.

You hold data back. Train on three quarters of the isolates, then score the model on the quarter it has never seen. If the model performs well on the unseen, it has learned something.

That only works if the held-back quarter is genuinely new to it.

Now look at what these isolates are.

NCBI sorts them into SNP clusters. Two genomes in the same cluster differ by a handful of bases out of 1.6 million. The same organism, sampled twice, often from two patients in one outbreak.

My 3,984 isolates fall into 1,128 clusters.

Split them at random and 344 of those clusters end up with isolates in both halves.

So the model learns from one genome and is then tested on its near-identical twin. It scores well and it has learned nothing. It is recognising a genome it has already seen.

That is data leakage. Information reaching the test set through a back door.

So I split again, keeping every cluster whole this time. No family on both sides.

344 became zero.

Then I repeated it twenty five times, because one quarter of 3,984 isolates can be accidentally kind.

0.9965 with families split. 0.9965 with families kept whole.

The number did not move.

Which is good news.

If the model was reading family resemblance, it would have failed on families it had never seen. It did not fail. The model predicted resistance in families it had never encountered, and that is what learning a mechanism looks like.

But yesterday's beta-lactamase is still sitting in there.

In the random split its weight was 1.469. In the grouped split, 2.364.

Splitting by family did not remove it.

It could not. blaOXA-493 and the gyrA mutation appear together in every isolate that carries either one. Divide the data any way you like and they are still together on both sides.

So there are two problems here and they are not the same problem.

Leakage is about how you divide the data. Dividing it better fixes it.

Confounding is about what is in the data. No division fixes it at all.

Tomorrow: the gene was there. The protein was not.

#MachineLearningForBiology #MachineLearning #ArtificialIntelligence #AntimicrobialResistance #Bioinformatics #ComputationalBiology

---

## Number and Metric Traceability

- **3,984 isolates across 1,128 SNP clusters:** `results/metrics/cohort_summary_ciprofloxacin.txt`.
- **Random split 344 clusters spanning both folds vs Grouped split 0 clusters spanning both folds:** `results/metrics/split_summary_ciprofloxacin.txt`.
- **0.9965 rule accuracy across 25 random vs 25 grouped folds:** `results/metrics/cross_validation_ciprofloxacin.txt` (`Random folds: rule 0.9965`, `Grouped folds: rule 0.9965`).
- **`blaOXA-493` logistic weight (1.469 random split vs 2.364 grouped split):** `results/metrics/model_results_ciprofloxacin.txt`.
