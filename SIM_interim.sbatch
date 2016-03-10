#!/bin/bash

#SBATCH --job-name=SIM_interim
#SBATCH --output=out/out.SIM
#SBATCH --error=error/err.SIM
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

#TMPDIR="/local-scratch/jdurnez/7606565/"
DATADIR10="/scratch/users/jdurnez/interim_SIM_10/"
DATADIR15="/scratch/users/jdurnez/interim_SIM_15/"
TMPDIR=$LOCAL_SCRATCH
HOMEDIR="/home/jdurnez/power_peak/"

#for j in $(seq 1 10);
#do
#  SIM=$(((i-1)*10+j))
SIM=$i
python $HOMEDIR/SIM_interim.py 10 30 $SIM $HOMEDIR $DATADIR10 $TMPDIR
python $HOMEDIR/SIM_interim.py 15 35 $SIM $HOMEDIR $DATADIR15 $TMPDIR
#done
