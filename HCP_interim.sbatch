#!/bin/bash

#SBATCH --job-name=HCP_interim
#SBATCH --output=out/out.HCP
#SBATCH --error=error/err.HCP
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

TMPDIR=$LOCAL_SCRATCH
RESDIRAD="/scratch/users/jdurnez/interim/hcp_adaptive/u2/"
RESDIRNAD="/scratch/users/jdurnez/interim/hcp_nonadaptive/u2/"
HOMEDIR="/home/jdurnez/power_peak/"
DATADIR="/scratch/PI/russpold/data/HCP/"

SIM=$i
python $HOMEDIR/HCP.py "2.3" "15" "80" $SIM $HOMEDIR $DATADIR $TMPDIR $RESDIRNAD
python $HOMEDIR/HCP_interim.py "2.3" "15" "80" $SIM $HOMEDIR $DATADIR $TMPDIR $RESDIRAD
