#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=20:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --partition=gpu    # gpu node(s)
#SBATCH --mail-user=emittman@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
export LD_LIBRARY_PATH="/usr:/usr/local:/usr/local/cuda-8.0:/usr/local/cuda-8.0/lib64"
module load R
cd /home/emittman/mm_simulations
R --vanilla CMD BATCH simulate_and_fit.R