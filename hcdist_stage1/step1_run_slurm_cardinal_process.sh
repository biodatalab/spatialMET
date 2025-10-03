#!/bin/bash
#SBATCH --job-name=cardinal_process
#SBATCH --cpus-per-task=24
#SBATCH --mem 240GB
#SBATCH --time=10:00:00
#SBATCH --output imzml_cardinal_process.out
#SBATCH --error imzml_cardinal_process.err

module load R

echo -ne "Slurm submission: "; echo $(date)

bash bash_scripts/go_cardinal_all.sh hcdist_imzml_param_list.txt

echo -ne "Slurm completed: "; echo $(date)
