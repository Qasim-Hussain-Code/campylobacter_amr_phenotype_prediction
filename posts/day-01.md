# Day 1: Provenance and Forensic Bioinformatics

Cancer patients at Duke University were assigned to chemotherapy regimens chosen by a machine learning model.

The model was wrong. Proving it took two statisticians thousands of hours, because the code and the processed data were never released.

The method, published in Nature Medicine in 2006, read the gene expression profile of a patient's tumour and predicted which drug that tumour would respond to. Oncology wants this badly. Three clinical trials were built on it.

What actually went wrong is the reason I am starting this series.

It was not the algorithm.

Keith Baggerly and Kevin Coombes, biostatisticians at MD Anderson, tried to reuse the method, could not reproduce it, and so reverse engineered the analysis from the published figures. They called the practice forensic bioinformatics. What they found were sample labels reversed between the responder and non-responder groups. Gene lists shifted by one position against their identifiers. 

Off-by-one errors. Mislabelled columns. Things that never appear in a methods section.

The trials were suspended.

The lesson generalises far beyond this case. In biology, most machine learning failures are not modelling failures. They are data failures wearing a model's clothes. An accuracy figure means nothing until you know the provenance of every label that produced it, and checking that requires understanding the assay, not the architecture.

Which brings me to why I am writing this.

On Kaggle I rank 68th of 10,761 for building biological datasets (a Grandmaster). I have open-sourced around sixty computational biology analyses. Metagenomic co-assemblies, pangenomes, single nucleus atlases, cross tissue TWAS, spatial transcriptomics. I can take raw reads to a merged Anvi'o profile without opening the documentation.

But, almost none of it is machine learning.

It is the half of the problem that Baggerly and Coombes were doing. Provenance, labels, matrix orientation. The other half is missing; building models. Both halves are required, and most people I meet are missing one of them.

So I am spending the next year(s) to help everyone close that gap. Foundations first, then classical methods on omics data, then deep learning, sequence models, generative and agentic systems. Code, notebooks, mistakes and dead ends, all open-sourced, all public.

Tell me in the comments which half you are missing.

#MachineLearningForBiology #MachineLearning #ArtificialIntelligence #Bioinformatics #ComputationalBiology

---

## Number and Metric Traceability

*This post serves as the conceptual introduction to the series on forensic bioinformatics and label provenance.*
- **Historical Case Study:** Duke University pharmacogenomics clinical trial error analysis (Baggerly & Coombes, 2009 *Annals of Applied Statistics*).
- **Dataset / Environment Audit:** Verified in `envs/environment.yml` and `scripts/01_fetch_ncbi_metadata.sh` that data provenance pipeline is fully reproducible from source.
