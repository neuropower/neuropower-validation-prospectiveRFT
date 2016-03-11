#!/bin/bash

#SBATCH --job-name=SIM_interim
#SBATCH --output=out/out.SIM
#SBATCH --error=error/err.SIM
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

#TMPDIR="/local-scratch/jdurnez/7606565/"
DATADIR="/scratch/users/jdurnez/interim_SIM/"
TMPDIR=$LOCAL_SCRATCH
HOMEDIR="/home/jdurnez/power_peak/"

for j in $(seq 1 10);
do
  SIM=$(((i-1)*10+j))
  #SIM=$i
  python $HOMEDIR/SIM_interim.py 15 60 $SIM $HOMEDIR $DATADIR $TMPDIR
done
