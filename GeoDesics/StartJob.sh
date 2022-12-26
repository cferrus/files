#!/bin/bash
. bin/this_machine.env || echo 'Not using env file.'
# Only applicable for bbh2 runs:
# For PBandJ (Perform Branching after Junk) only EccRedLev should be started
# with StartJob.sh. ProhibitStartJobReruns.txt gets generated at Eccentricity
# Reduction restart, and safeguards against starting Levs != EccRedLev from t=0
# for PBandJ.
if test -f ProhibitStartJobReruns.txt; then
    echo
    echo "ERROR: ProhibitStartJobReruns.txt exists. See below:"
    cat ProhibitStartJobReruns.txt
    echo
    exit 1
fi

bin/DoMultipleRuns -L -c 'sbatch ./Submit.sh'
