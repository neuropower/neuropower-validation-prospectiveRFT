#!/bin/bash

#SBATCH --job-name=SIM.PRE.CS
#SBATCH --output=error/out.SIM.PRE.CS
#SBATCH --error=error/err.SIM.PRE.CS
#SBATCH --time=40:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

. ./config.sh

PILOT=15
FINAL=40
EXC="2.3"
ADAPTIVE="predictive"
MODEL="CS"

OUTDIR=$(echo $RESDIR\SIM_$ADAPTIVE\_$PILOT\_$EXC\_$MODEL)

if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi

#for j in $(seq 1 3);
#do
#SIM=$(((i-1)*3+j))
SEED=$i
python $SCRIPTDIR/SIM.py $EXC $PILOT $FINAL $SEED $SIMFILEDIR $OUTDIR $TMPDIR $ADAPTIVE $MODEL
#done
