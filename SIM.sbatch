#!/bin/bash

#SBATCH --job-name=SIM_power
#SBATCH --output=out/out.SIM
#SBATCH --error=error/err.SIM
#SBATCH --time=20:00:00
#SBATCH --nodes=1
########SBATCH --partition=russpold

DATADIR10="/scratch/users/jdurnez/power_peak_SIM_10/"
DATADIR15="/scratch/users/jdurnez/power_peak_SIM_15/"
TMPDIR=$LOCAL_SCRATCH
HOMEDIR="/home/jdurnez/power_peak/"

for j in $(seq 1 10);
do
  SIM=$(((i-1)*10+j))
  python $HOMEDIR/SIM.py 10 30 $SIM $HOMEDIR $DATADIR10 $TMPDIR
  python $HOMEDIR/SIM.py 15 35 $SIM $HOMEDIR $DATADIR15 $TMPDIR
done
