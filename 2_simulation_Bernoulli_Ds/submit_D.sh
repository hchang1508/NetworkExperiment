#!/bin/bash
#SBATCH --job-name=sample
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haoge.chang@yale.edu
#SBATCH --partition=day
#SBATCH --time=24:00:00
#SBATCH --array=4
#SBATCH --requeue
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
conda deactivate
conda activate network_experiment

cd "/home/hc654/Unified/final_analysis"
#module load R/3.6.1-foss-2018b-X11-20180604

Rscript --vanilla extract_D.R 1 $SLURM_ARRAY_TASK_ID 
