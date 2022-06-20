#!/bin/bash
#SBATCH --time=04:00:00                  # Job run time (hh:mm:ss)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=20             # Number of task (cores/ppn) per node
#SBATCH --job-name=overlaps          # Name of batch job
#SBATCH --partition=secondary-fdr           # Partition (queue)
#SBATCH --output=overlap-calculation.o%j            # Name of batch job output file
##SBATCH --error=overlap-calculation.e%j            # Name of batch job error file
##SBATCH --mail-user=wentzel4@illinois.edu  # Send email notifications
##SBATCH --mail-type=BEGIN,END           # Type of email notifications to send
#
###############################################################################
# Change to the directory from which the batch job was submitted
# Note: SLURM defaults to running jobs in the directory where
# they are submitted, no need for cd'ing to $SLURM_SUBMIT_DIR

#cd ${SLURM_SUBMIT_DIR}


# Load Mathematica module (Enable Mathematica in user environment)
module load mathematica/12

# Run Mathematica script
echo "Running job number " ${SLURM_ARRAY_TASK_ID}

math -noprompt -run < mech-em-overlap-solver.wls 1.0 0.01 "/home/wentzel4/GWs/GWs-Mechanical-Coupling/CoupledSpheres/radius1-thickness1percent/" "NbCavity-r=1.0-1.01.txt" ${SLURM_ARRAY_TASK_ID}


