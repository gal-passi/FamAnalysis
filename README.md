   # Family Analysis Pipeline

This repository contains code and tools developed for the paper Discovering Predisposing Genes for Hereditary Breast Cancer Using Deep Learning (https://academic.oup.com/bib/article/25/4/bbae346/7717952).

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
3. Normalize models scores and rank variants using an aggregated score - DSrank (see paper).

# Setup

## Hardware requirements

The entire pipeline can run on standard CPU. To install the database and precomputed scores, make sure to have at least 23GB free on disk. 

## Software requirements

FamAnalysis requires python 3.10

## Installation 

After cloning the project, run:

`pip install pandas requests tqdm psutil click wget bio ` 

`python setup.py`

In `definition.py`, make sure to set `CONTACT` to your email address. This is required for UniProtKB requests.

**Note:** In case of an unstable internet connection, `setup.py` may raise a `ConnectionError` while downloading large files. 
Simply re-run the program from the method that crashed. The download will resume from the point where it stopped.

# Running FamAnalysis 
Run the following command to get help:

`python main.py --help`

Protein variants should be provided in the format shown in `example.csv`.

**Note:** DNA data (i.e., Chr, Start, End) is not required. If these fields are not available, please fill them with any arbitrary value.

## Initializing Database

To build protein and mutation objects for all variants run: 

`python main.py --data-path example.csv --action to-csv`

To limit CPUs used, run:

`python main.py --data-path example.csv --action to-csv --workers 16`

**Note:** When using > 16 CPUs or if internet connection is unstable it may be required to increase `TIMEOUT`, `WAIT_TIME` and `RETRIES` values in `definitions.py`

## Calculating Models Scores

To calculate model scores for all variants, run:

`python main.py --data-path example.csv --action score-AFM`

`python main.py --data-path example.csv --action score-ESM`

`python main.py --data-path example.csv --action score-EVE`

To re-calculate scores for variants that were previously scored, run:

`python main.py --data-path example.csv --action score-EVE --recalc 1`

To use only EVEmodel without using [CPT](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03024-6) imputation, run:

`python main.py --data-path example.csv --action score-EVE --use-cpt 0`

To set offset for mutation position, run:

`python main.py --data-path example.csv --action score-EVE --offset 1`

**Note** changing the offset from default is not recommended and should preferably be handled at the input file. 

## Ranking Scores

FamAnalysis uses DSRank to aggregate the model scores. The normalized scores of the three ML models’ predictions are averaged over all available scores per variant. 

To rank all variant, run:

`python main.py --data-path example.csv --action rank-DS `

use `--ds-thr` flag to set the threshold for the minimal number of model scores needed to calculate DSRank.

## Results

To create a CSV file summarizing all scores, run:

main.py --data-path example.csv --action to-csv --out-path results.csv 

use `--include-type 1` to include esm and eve score types in output file. 

## Reference

If you use this code, please cite our [paper](https://academic.oup.com/bib/article/25/4/bbae346/7717952):


