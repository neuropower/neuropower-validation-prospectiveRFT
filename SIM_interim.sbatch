#!/bin/bash

#SBATCH --job-name=FPR
#SBATCH --output=out/out.SIM
#SBATCH --error=error/err.SIM
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

#TMPDIR="/local-scratch/jdurnez/7606565/"
NADIR="/scratch/users/jdurnez/interim/nonadaptive/"
ADIR="/scratch/users/jdurnez/interim/adaptive/"
TMPDIR=$LOCAL_SCRATCH
HOMEDIR="/home/jdurnez/power_peak/"

for j in $(seq 1 3);
do
  SIM=$(((i-1)*3+j))
  #SIM=$i
  python $HOMEDIR/SIM_interim.py 15 61 $SIM $HOMEDIR $ADIR $TMPDIR
  python $HOMEDIR/SIM.py 15 61 $SIM $HOMEDIR $NADIR $TMPDIR
done
