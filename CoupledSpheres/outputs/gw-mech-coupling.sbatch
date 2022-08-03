#!/bin/bash
#SBATCH --time=04:00:00                  # Job run time (hh:mm:ss)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=20             # Number of task (cores/ppn) per node
#SBATCH --job-name=gwmech-coupling          # Name of batch job
#SBATCH --partition=secondary           # Partition (queue)
#SBATCH --output=gw-mech-coupling.o%j            # Name of batch job output file
#SBATCH --error=gw-mech-coupling.e%j            # Name of batch job error file
#SBATCH --mail-user=wentzel4@illinois.edu  # Send email notifications
#SBATCH --mail-type=BEGIN,END           # Type of email notifications to send
#
###############################################################################
# Change to the directory from which the batch job was submitted
# Note: SLURM defaults to running jobs in the directory where
# they are submitted, no need for cd'ing to $SLURM_SUBMIT_DIR

#cd ${SLURM_SUBMIT_DIR}


# Load Mathematica module (Enable Mathematica in user environment)
module load mathematica/12

# Run Mathematica script
echo "Running job number "

math -noprompt -run < gw-mech-overlap-solver.wls 0.2 0.005 "/home/wentzel4/GWs/GWs-Mechanical-Coupling/CoupledSpheres/radius0.2-thickness0.005/" "NbCavity-r=0.2-0.005.txt"


