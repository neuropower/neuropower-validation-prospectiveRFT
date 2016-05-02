#!/bin/bash

#SBATCH --job-name=SIM.predictive
#SBATCH --output=error/out.SIM.predictive
#SBATCH --error=error/err.SIM.predictive
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

. ./config.sh

PILOT=10
FINAL=65
U="2.3"

SCRIPTDIR=$(echo $HOMEDIR\PythonScripts/)
FILEDIR=$(echo $HOMEDIR\SimulationFiles/)
OUTDIR=$(echo $RESDIR\SIM_predictive_$PILOT\_$U/)

if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi

#for j in $(seq 1 3);
#do
#SIM=$(((i-1)*3+j))
SIM=$i
python $SCRIPTDIR/SIM.py $U $PILOT $FINAL $SIM $FILEDIR $RESDIR $TMPDIR "predictive"
#done
