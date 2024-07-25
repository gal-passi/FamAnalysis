# Family Analysis Pipeline

This repository contains code and tools developed for the paper Discovering predisposing genes for hereditary breast cancer using deep learning (https://academic.oup.com/bib/article/25/4/bbae346/7717952).

# Overview

The FamAnalysis pipeline processes large datasets of missense variants and ranks them based on their pathogenicity compared to other variants in the dataset.
We provide here a user-friendly command line interface of the DSRank methodology outlined in our [paper](https://academic.oup.com/bib/article/25/4/bbae346/7717952). We hope to extend the interface to include FamRank in the near future.
![pipeline_](https://github.com/user-attachments/assets/aff5dcc9-cfdd-4643-9fc3-e56a53375278)

**Updates to current version:** 
1. we have now replaced [FIRM](https://academic.oup.com/nar/article/47/13/6642/5523008) with [AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492) precomputed scores.
2. We have also replaced [ESM-1v](https://www.biorxiv.org/content/10.1101/2021.07.09.450648v2) inference with precomputed scores obtained from [ESM-1b](https://www.biorxiv.org/content/10.1101/2021.02.12.430858v2)
   as described in the study [Genome-wide prediction of disease variant effects with a deep protein language model](https://www.nature.com/articles/s41588-023-01465-0) - official GitHub repository (https://github.com/ntranoslab/esm-variants.git).
4. The current version can run on standard CPUs and does not require direct model inference.

**The pipeline consists of 3 stages:**
1. Create protein and mutations objects based on input data by pulling all relevant metadata from UniProtKB. 
2. Calculate pathogenicity scores using ESM, AlphaMissense and [EVEmodel](https://www.nature.com/articles/s41586-021-04043-8).
3. Normalize models scores and rank variants using an aggregated acore - DSrank (see paper).
