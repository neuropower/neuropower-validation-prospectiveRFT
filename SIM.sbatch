#!/bin/bash

#SBATCH --job-name=SIM_power
#SBATCH --output=out/out.SIM
#SBATCH --error=error/err.SIM
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

DATADIR10="/scratch/users/jdurnez/power/SIM_10/"
DATADIR15="/scratch/users/jdurnez/power/SIM_15/"
DATADIR20="/scratch/users/jdurnez/power/SIM_20/"
DATADIR25="/scratch/users/jdurnez/power/SIM_25/"
TMPDIR=$LOCAL_SCRATCH
HOMEDIR="/home/jdurnez/power_peak/"

SIM=$i
for j in $(seq 1 2);
do
  SIM=$(((i-1)*2+j))
  python $HOMEDIR/SIM.py "2.3" 10 60 $SIM $HOMEDIR $DATADIR10 $TMPDIR
  python $HOMEDIR/SIM.py "2.3" 15 60 $SIM $HOMEDIR $DATADIR15 $TMPDIR
  python $HOMEDIR/SIM.py "2.3" 20 60 $SIM $HOMEDIR $DATADIR20 $TMPDIR
  python $HOMEDIR/SIM.py "2.3" 25 60 $SIM $HOMEDIR $DATADIR25 $TMPDIR
done
