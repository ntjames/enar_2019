#!/bin/bash
#SBATCH --job-name=cbmodsim
#SBATCH --mail-user=nathan.t.james@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:05:30
#SBATCH --mem=2G
#SBATCH --constraint=sandybridge|haswell
#SBATCH --array=1-5
#SBATCH --output=cbout/cbmod_sim_%a.out

module restore R_basic
R --version 

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

## run this first to compile STAN models
## Rscript --no-save cbmod_init.R

## run array job over sim ids
## there are 10800 total sim ids
Rscript --no-save cbmod.R $SLURM_ARRAY_TASK_ID
