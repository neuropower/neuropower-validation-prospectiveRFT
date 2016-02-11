#!/bin/bash

#SBATCH --job-name=HCP_power
#SBATCH --output=out/out.HCP
#SBATCH --error=error/err.HCP
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold

TMPDIR=$LOCAL_SCRATCH
#TMPDIR="/local-scratch/jdurnez/7316018/"
RESDIR10="/scratch/users/jdurnez/power_peak_HCP_10/"
RESDIR15="/scratch/users/jdurnez/power_peak_HCP_15/"
HOMEDIR="/home/jdurnez/power_peak/"
DATADIR="/scratch/PI/russpold/data/HCP/"

SIM=$i
#for j in $(seq 1 3);
#do
#  SIM=$(((i-1)*3+j))
python $HOMEDIR/HCP.py "15" "65" $SIM $HOMEDIR $DATADIR $TMPDIR $RESDIR15
python $HOMEDIR/HCP.py "10" "70" $SIM $HOMEDIR $DATADIR $TMPDIR $RESDIR10
#done
