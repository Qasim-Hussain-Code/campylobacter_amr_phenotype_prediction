# Day 4: Confounding and Passenger Genes

Chapter 1: Predicting antibiotic resistance from a genome

Yesterday the model paid a beta-lactamase the same as a gyrase mutation.

So what does a beta-lactamase do?

It cuts open a beta-lactam ring. Penicillin has one. Amoxicillin has one.

Ciprofloxacin has none. There is nothing there for the enzyme to cut.

The gene is blaOXA-493. Five isolates carry it. All five are resistant.

So the model paid it.

Now look at those five again.

Every one also carries a gyrase mutation. All five. No exceptions.

Put yourself where the model sits. You see two columns. They rise and fall together. Where both are present, the isolate is resistant.

Which one is responsible?

You cannot tell. Nothing in the table can tell you.

So the model split the credit. 2.364 each. It had no other option.

That is confounding, and it is the oldest problem in statistics wearing a genomic coat.

Now the part that should worry you.

Suppose a sixth isolate turns up. It carries blaOXA-493 and no gyrase mutation.

The model calls it resistant. The model is wrong.

There is no sixth isolate. Not in these 3,984.

The mistake is already in the model. Nothing in the data will show it to you.

The model did not learn what the gene does. It learned what the gene travels with.

Why they travel together, I cannot tell you from five isolates. Shared ancestry, perhaps. Or nothing at all.

The model cannot tell you either. It only ever saw them arrive as a pair.

Which is the limit of the whole method.

A model sees co-occurrence. That is all it ever sees. Causes and companions look identical from inside a table.

What separates them is knowing what the gene does, and that never comes from the data.

It comes from biology.

It comes from somebody who knows ciprofloxacin has no beta-lactam ring.

If the model is reading company rather than mechanism, there is a way to catch it.

Hide whole families of isolates from it. Then see whether it still works.

Tomorrow.

#MachineLearningForBiology #MachineLearning #ArtificialIntelligence #AntimicrobialResistance #Bioinformatics #ComputationalBiology

---

## Number and Metric Traceability

- **`blaOXA-493` co-occurrence (5 isolates carry it, 5 resistant, 100% carrying gyrA 86 rule):** `results/metrics/cross_validation_ciprofloxacin.txt` (Co-occurrence table: `blaOXA-493 n=5 P(rule|feat)=1.000`) & `results/metrics/variant_table_ciprofloxacin.txt` (`blaOXA-493  5  5  100.0%`).
- **Equal coefficient weight (+2.364):** `results/metrics/model_results_ciprofloxacin.txt` (Grouped split: `+2.364 gyrA_T86V=POINT` and `+2.364 blaOXA-493`).
- **Contingency table counts (Both=5, Mutation alone=818, Gene alone=0, Neither=3,161):** Derived from `results/metrics/cross_validation_ciprofloxacin.txt` and total cohort n=3,984 (Resistant=824, Susceptible=3,160).
