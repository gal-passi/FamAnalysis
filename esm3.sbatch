#!/bin/bash
#SBATCH --job-name=esm3            
#SBATCH --output=job_%a.out     
#SBATCH --error=job_%a.err      
#SBATCH --time=01:00:00                 
#SBATCH --cpus-per-task=8      
#SBATCH --mem=16G                       
#SBATCH --array=1-$(python -c "import main; print(len(main.protein_names()))")%100

# Run the Python script with the current task ID
conda activate famanalysis
python main.py --data-path input.csv --action score-ESM3 --token <token> --recalc 1 --optimized 1 --task_id $SLURM_ARRAY_TASK_ID
