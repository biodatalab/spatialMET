#!/bin/bash
#SBATCH --job-name=hcdist_pixel_tree
#SBATCH --cpus-per-task=32
#SBATCH --time=10:00:00
#SBATCH --mem=360GB
#SBATCH --output hcdist_pixel_tree_%A_%a.out
#SBATCH --output hcdist_pixel_tree_%A_%a.err
#SBATCH -a 1-2%2

args=$( head -n ${SLURM_ARRAY_TASK_ID} hcdist_pixel_tree_param_list.txt | tail -n1 )

arg1=$( echo "$args" | cut -f1 )
arg2=$( echo "$args" | cut -f2 )
arg3=$( echo "$args" | cut -f3 | tr -d '\n' | tr -d '\r' )

echo -ne "Slurm submission: "; echo $(date)

bash bash_scripts/go_detect_spatial_domains.sh "$arg1" "$arg2" "$arg3"

echo "End analysis: "
echo $(date)
