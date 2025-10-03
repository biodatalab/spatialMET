#!/bin/bash
#SBATCH --job-name=hcdist_result_files
#SBATCH --cpus-per-task=22
#SBATCH --time=10:00:00
#SBATCH --mem=120GB
#SBATCH --output hcdist_result_files_%A_%a.out
#SBATCH --error hcdist_result_files_%A_%a.err
#SBATCH -a 1-6%6

args=$( head -n ${SLURM_ARRAY_TASK_ID} hcdist_result_param_list.txt | tail -n1 )

arg1=$( echo "$args" | cut -f1 )
arg2=$( echo "$args" | cut -f2 )
arg3=$( echo "$args" | cut -f3 )
arg4=$( echo "$args" | cut -f4 )
arg5=$( echo "$args" | cut -f5 )
arg6=$( echo "$args" | cut -f6 | tr -d '\n' | tr -d '\r' )

echo -ne "Slurm submission: "; echo $(date)

bash bash_scripts/go_generate_hcdist_result_files.sh "$arg1" "$arg2" "$arg3" "$arg4" $arg5 $arg6

echo -ne "Slurm completed: "; echo $(date)
