#!/bin/bash

#SBATCH --job-name=biascorrection
#SBATCH --output=out/out.SIM
#SBATCH --error=error/err.SIM
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

#TMPDIR="/local-scratch/jdurnez/7786729/"
DATADIR="/scratch/users/jdurnez/interim/bias/"
HOMEDIR="/home/jdurnez/power_peak/"
TMPDIR=$LOCAL_SCRATCH

for j in $(seq 1 30);
do
  SIM=$(((i-1)*30+j))
  python $HOMEDIR/SIM_interim_biascorr.py 15 100 $SIM "0.8" 0 $HOMEDIR $DATADIR $TMPDIR
  python $HOMEDIR/SIM_interim_biascorr.py 15 100 $SIM "0.8" 1 $HOMEDIR $DATADIR $TMPDIR
done
