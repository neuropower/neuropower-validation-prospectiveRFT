#!/bin/bash

#SBATCH --job-name=SIM.interim
#SBATCH --output=error/out.SIM.interim
#SBATCH --error=error/err.SIM.interim
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

. ./config.sh

PILOT=15
FINAL=65
U="2.3"

SCRIPTDIR=$(echo $HOMEDIR\PythonScripts/)
FILEDIR=$(echo $HOMEDIR\SimulationFiles/)
OUTDIR_A=$(echo $RESDIR\SIM_adaptive_$PILOT\_$U/)
OUTDIR_P=$(echo $RESDIR\SIM_predictive_$PILOT\_$U/)

if [ ! -d "$OUTDIR_A" ]; then mkdir $OUTDIR_A; fi
if [ ! -d "$OUTDIR_P" ]; then mkdir $OUTDIR_P; fi

for j in $(seq 1 5);
do
  SIM=$(((i-1)*5+j))
  SIM=$i
  python $SCRIPTDIR/SIM.py $U $PILOT $FINAL $SIM $FILEDIR $OUTDIR_P $TMPDIR "predictive"
  python $SCRIPTDIR/SIM.py $U $PILOT $FINAL $SIM $FILEDIR $OUTDIR_A $TMPDIR "adaptive"
done
