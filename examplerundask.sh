#!/bin/bash

#SBATCH --mem=32G
#SBATCH -c2
#SBATCH --time=24:0:0
#SBATCH --account="dina"

source ~/.bashrc
conda activate /cs/labs/dina/haimasree/condaenvs/famanalysis/

# pip install pandas requests tqdm psutil click wget bio torch fair-esm
# rm -rf DB
# python setup.py

time python runafmwithdask.py

conda deactivate
