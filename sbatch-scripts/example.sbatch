#!/bin/bash
#SBATCH --time=04:00:00                  # Job run time (hh:mm:ss)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=40             # Number of task (cores/ppn) per node
#SBATCH --job-name=GWsCalc          # Name of batch job
#SBATCH --partition=secondary           # Partition (queue)
#SBATCH --output=GWsCalc.o%j            # Name of batch job output file
##SBATCH --error=GWsCalc.e%j            # Name of batch job error file
##SBATCH --mail-user=wentzel4@illinois.edu  # Send email notifications
##SBATCH --mail-type=BEGIN,END           # Type of email notifications to send
#
###############################################################################
# Change to the directory from which the batch job was submitted
# Note: SLURM defaults to running jobs in the directory where
# they are submitted, no need for cd'ing to $SLURM_SUBMIT_DIR

#cd ${SLURM_SUBMIT_DIR}


# Load Mathematica module (Enable Mathematica in user environment)
module load mathematica

# Run Mathematica script
/usr/bin/time -v math -noprompt -run < ppf-units.wls
