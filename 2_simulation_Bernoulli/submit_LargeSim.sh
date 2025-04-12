#!/bin/bash
#SBATCH --job-name=sample
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haoge.chang@yale.edu
#SBATCH --partition=scavenge
#SBATCH --time=24:00:00
#SBATCH --array=1-3000
#SBATCH --requeue
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
conda deactivate
conda activate network_experiment

cd "/home/hc654/Unified/scripts_grace/final_analysis"
#module load R/3.6.1-foss-2018b-X11-20180604
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 4
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 5
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0

Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 1
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0.1
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0.2
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0.3
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0.4
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0.5
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0.6
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0.7

#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0.8
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 0.9
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 3 1 
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 7
#Rscript --vanilla LargeSim.R $SLURM_ARRAY_TASK_ID 1 2
