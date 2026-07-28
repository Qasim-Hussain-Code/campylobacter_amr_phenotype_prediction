# Day 2: The Biological Target

Chapter 1: Predicting antibiotic resistance from a genome

Campylobacter has 1.6 million bases of DNA. Change one of them, a single C to a T, and ciprofloxacin stops working.

One letter out of 1.6 million.

That is not a figure of speech.

Fluoroquinolones kill bacteria by jamming the enzymes that untangle DNA during replication.

E. coli has two of them, DNA gyrase and topoisomerase IV. Mutate one and the drug still has the other, which is why resistance in E. coli usually needs hits in both.

Campylobacter jejuni has no topoisomerase IV. Repeated attempts to find parC have failed. The genome does not carry it.

One target. No backup.

A single base change at position 257 of gyrA, turning threonine 86 into isoleucine, raises the ciprofloxacin MIC 128-fold. 

But, appearing is one thing. Lasting is another.

Most human Campylobacter infection comes from poultry. The question is not whether the mutation survives inside a patient. It is whether it survives inside a chicken, with no drug anywhere near it.

It does. Resistant and susceptible strains both colonised chickens equally well with no antibiotic present. Nothing pushes the mutation back out. Once it appears, it stays.

So the answer is in the genome. Known gene, known position.

Meanwhile the clinic still does it the slow way. Grow the organism, expose it to the drug, wait.

Campylobacter is microaerophilic and slow. That is days, and the patient is treated on a guess in the meantime.

Sequencing takes a day. Read position 257 and you are done.

But you are not.

Two isolates. Both carry the change at 257. One dies at a concentration of ciprofloxacin that the other survives.

Why?

Campylobacter has a pump, CmeABC, that pushes the drug back out of the cell before it ever reaches the gyrase. Disable the pump and the resistance collapses, mutation still intact.

Position 257 is not a switch. It is a contribution, and how much it contributes depends on what else the genome is carrying.

Which is the question this chapter is built on.

If resistance is written in the sequence, why can nobody simply read it off?

Tomorrow: the second half of that problem, which is not in the biology. It is in the machine learning.

#MachineLearningForBiology #MachineLearning #ArtificialIntelligence #AntimicrobialResistance #Bioinformatics #ComputationalBiology

---

Numbers in this post

```
gyrA T86I as the marker       results/metrics/cohort_summary_ciprofloxacin.txt
                              (Marker: gyrA_T86I)
Release                       PDG000000003.2859
                              results/metrics/provenance_ciprofloxacin.txt

The mechanism in this post is published biology, not an output of this
pipeline:

1.6 million bases             C. jejuni NCTC 11168 reference genome
position 257 of gyrA          codon 86, fluoroquinolone literature
128-fold rise in MIC          fluoroquinolone literature
no topoisomerase IV, no parC  C. jejuni genome literature
CmeABC efflux                 C. jejuni efflux literature
```
