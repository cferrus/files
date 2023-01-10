#!/bin/bash -
#SBATCH -J LensingQ3                   # Job Name
#SBATCH -o SpEC.stdout                # Output file name
#SBATCH -e SpEC.stderr                # Error file name
#SBATCH -n 128                  # Number of cores
#SBATCH --ntasks-per-node 128        # number of MPI ranks per node
#SBATCH -p cca                  # Queue name
#SBATCH -t 24:0:00   # Run time
#SBATCH -A cca                # Account name
#SBATCH --no-requeue
#SBATCH --constraint=ib

# DO NOT MANUALLY EDIT THIS FILE! See `MakeSubmit.py update -h`.
# This is for submitting a batch job on 'rusty'.
umask 0022
. bin/this_machine.env || echo 'Not using env file.'
set -x
export PATH=$(pwd -P)/bin:$PATH


EvolveGeodesicsWrapper -a="celestemelizef@gmail.com" -f=""
