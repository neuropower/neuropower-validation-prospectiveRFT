#!/bin/bash

#SBATCH --job-name=SIM_interim
#SBATCH --output=out/out.SIM
#SBATCH --error=error/err.SIM
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

#TMPDIR="/local-scratch/jdurnez/7733199/"
DATADIR10="/scratch/users/jdurnez/interim_SIM_biascorr_10/"
DATADIR15="/scratch/users/jdurnez/interim_SIM_biascorr_15/"
HOMEDIR="/home/jdurnez/power_peak/"
TMPDIR=$LOCAL_SCRATCH

for j in $(seq 1 50);
do
  SIM=$(((i-1)*10+j))
  python $HOMEDIR/SIM_interim_biascorr.py 10 100 $SIM "0.8" $HOMEDIR $DATADIR10 $TMPDIR
  python $HOMEDIR/SIM_interim_biascorr.py 15 100 $SIM "0.8" $HOMEDIR $DATADIR15 $TMPDIR
done
