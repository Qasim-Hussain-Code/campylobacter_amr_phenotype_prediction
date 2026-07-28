# Day 3: Features, Labels, and Baselines

Chapter 1: Predicting antibiotic resistance from a genome

Yesterday I said position 257 is not a switch.

Today I put it to the test.

I took 3,984 Campylobacter genomes from NCBI. Every one had a ciprofloxacin result measured on a plate.

Then I asked one question of each. Is there a substitution at gyrA position 86.

Nothing else. No pump, no context, no second gene.

It was right for all but thirteen of them. 99.67%.

The pump is real. Two isolates with the same mutation do survive different doses of the drug.

I modelled none of it. I asked one question and it was right 3,971 times.

The reason has nothing to do with the bacterium.

It is about what the laboratory measured.

Not the word resistant. A number.

The organism goes into a row of tubes. Each holds twice as much ciprofloxacin as the one before it.

The lowest concentration that stops it growing is the result.

That number is continuous. Some stop at a low dose, some need far more.

Then somebody draws a line across it, susceptible below it, resistant above it.

A good pump buys the organism another tube or two. 

The pump did something. The label recorded nothing.

I did not solve yesterday's problem. I answered a newer one.

Before that 99.67% impresses anybody, including me.

Try to guess susceptible every time. Read no DNA at all. Most Campylobacter is susceptible, so you get 79.32%.

That is the floor. Find the floor before you believe any number.

One question about one position beats it by twenty points.

I gave a second model all 84 resistance genes and let it weight them freely.

It did not beat the one question.

Which is the lesson, and it arrives before any algorithm.

What a model looks at are its features. What it predicts is its label.

And the label is a decision somebody made, usually before you arrived, often for convenience.

Choose the label and you have fixed the ceiling. Every model after that, however clever, works underneath it.

Everything here arrives with its answers attached. That is supervised learning.

The 84-gene model gave position 86 a weight of 8.9.

It gave a beta-lactamase 2.364.

It gave a real gyrase mutation 2.364.

The same number, to three decimals, for a gene that does nothing to fluoroquinolones.

Why? Tomorrow.

#MachineLearningForBiology #MachineLearning #ArtificialIntelligence #AntimicrobialResistance #Bioinformatics #ComputationalBiology

---

Numbers in this post

```
3,984 genomes                 results/metrics/cohort_summary_ciprofloxacin.txt
99.67%, all but thirteen      results/metrics/marker_comparison_ciprofloxacin.txt
                              (row: any change at gyrA 86, acc 0.9967, err 13)
79.32% guessing susceptible   results/metrics/marker_comparison_ciprofloxacin.txt
                              (printed as 0.7932; cohort_summary rounds to 0.793)
all 84 resistance genes       results/metrics/threshold_sweep_ciprofloxacin.txt
                              (vocabulary, rarity threshold 1)
weight of 8.9                 results/metrics/model_results_ciprofloxacin.txt
2.364, twice                  (grouped split; 8.9 is printed as +8.845)
Release                       PDG000000003.2859

Arithmetic from the above, printed in no file:

right 3,971 times             3,984 isolates minus 13 errors
beats it by twenty points     99.67 minus 79.32
```
