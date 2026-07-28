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

Numbers in this post

```
3,984 isolates                results/metrics/cohort_summary_ciprofloxacin.txt
1,128 clusters                results/metrics/cohort_summary_ciprofloxacin.txt
train on three quarters       results/metrics/split_summary_ciprofloxacin.txt
                              (test_fraction=0.25, seed=42)
344 clusters in both halves   results/metrics/split_summary_ciprofloxacin.txt
344 became zero               results/metrics/split_summary_ciprofloxacin.txt
                              (grouped split: 0 clusters spanning both folds)
repeated twenty five times    results/metrics/cross_validation_ciprofloxacin.txt
0.9965 split, 0.9965 whole    (25 folds per scheme; the cross-validated rule
                              covers T86A, T86I and T86V)
1.469 in the random split     results/metrics/model_results_ciprofloxacin.txt
2.364 in the grouped split    results/metrics/model_results_ciprofloxacin.txt
Release                       PDG000000003.2859
```
