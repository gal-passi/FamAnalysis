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

time python main.py --data-path input.csv --action score-EVE --use-cpt 0
time python main.py --data-path input.csv --action score-EVE --use-cpt 1
python main.py --action to-csv --out-path results_test_data_EVE.csv

time python main.py --data-path input.csv --action score-ESM
python main.py --action to-csv --out-path results_test_data_ESM1b.csv

time python main.py --data-path input.csv --action score-AFM
python main.py --action to-csv --out-path results_test_data_AFM.csv

time python main.py --action score-ESM3 --token <token>
python main.py --action to-csv --out-path results_test_data_ESM3.csv


conda deactivate
