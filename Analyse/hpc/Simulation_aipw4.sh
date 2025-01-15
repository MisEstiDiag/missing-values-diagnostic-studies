#!/bin/bash
# Do not forget to select a proper partition if the default
# one is no fit for the job! You can do that either in the sbatch
# command line or here with the other settings.
#SBATCH --partition=std
#SBATCH --job-name=Sim_aipw4
#SBATCH --nodes=7
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
# Never forget that! Strange happenings ensue otherwise.
#SBATCH --export=NONE

set -e # Good Idea to stop operation on first error.

source /sw/batch/init.sh

# Load environment modules for your application here.
module switch env env/2021Q2-gcc-openmpi
module load R/4.1.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK          
export OMP_SCHEDULE=static                   
# export OMP_DISPLAY_ENV=verbose              
# export KMP_AFFINITY=verbose                  


# Actual work starting here. You might need to call
# srun or mpirun depending on your type of application
# for proper parallel work.	

Rscript /home/bbb2836/Simulation_aipw_idx4.R


exit