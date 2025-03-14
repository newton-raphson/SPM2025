#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=10      # number of nodes
#SBATCH --ntasks-per-node=36  # 36 processor core(s) per node
###SBATCH --mem=100G       # maximum memory per node
#SBATCH --constraint=nova22
#SBATCH --job-name="plane"   # job name
##SBATCH --mail-user=samundra@iastate.edu # email address
##SBATCH --mail-type=BEGIN
##SBATCH --mail-type=END
##SBATCH --mail-type=FAIL
##SBATCH --reservation=mech-ai

#SBATCH -o plane_spm25.o%j
#SBATCH -e plane_spm25.e%j

# Load necessary modules
module purge
module load intel
module load intel-mkl

EXECUTABLE="/work/mech-ai-scratch/samundra/projects/SPM2025/build_3D/projects/Elasticity-SBM/sle-sbm"
comd="mpirun -np $SLURM_NTASKS ${EXECUTABLE}"
${comd} 