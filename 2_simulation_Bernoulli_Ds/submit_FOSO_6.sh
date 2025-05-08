#!/bin/bash
#SBATCH --job-name=sample
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haoge.chang@yale.edu
#SBATCH --partition=scavenge
#SBATCH --time=24:00:00
#SBATCH --array=1-10000
#SBATCH --requeue
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
conda deactivate
conda activate network_experiment

cd "home/hc654/NetworkExperiment/2_simulation_Bernoulli_Ds"
#module load R/3.6.1-foss-2018b-X11-20180604
Rscript --vanilla 0_compute_FOSO_all_clusters.R 3 4  $SLURM_ARRAY_TASK_ID 6
Rscript --vanilla 0_compute_FOSO_all_clusters.R 3 5  $SLURM_ARRAY_TASK_ID 6
Rscript --vanilla 0_compute_FOSO_all_clusters.R 3 6  $SLURM_ARRAY_TASK_ID 6
Rscript --vanilla 0_compute_FOSO_all_clusters.R 4 5  $SLURM_ARRAY_TASK_ID 6
Rscript --vanilla 0_compute_FOSO_all_clusters.R 4 6  $SLURM_ARRAY_TASK_ID 6
