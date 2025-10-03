#!/bin/bash
#SBATCH --job-name=hcdist_filter
#SBATCH --cpus-per-task=24
#SBATCH --mem 240GB
#SBATCH --time=10:00:00
#SBATCH --output hcdist_filter.out
#SBATCH --error cdist_filter.err

module load R

echo -ne "Slurm submission: "; echo $(date)

bash bash_scripts/go_filter_data.sh hcdist_filter_param_list.txt

echo -ne "Slurm completed: "; echo $(date)
