#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=5GB
#SBATCH --job-name=producePlotsSlow
trap "echo recieved SIGUSR1;" SIGUSR1;
R CMD BATCH --no-save --no-restore -- producePlotsSlow.R producePlotsSlow.Rout.$SLURM_JOBID
