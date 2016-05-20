#!/bin/bash

#SBATCH --job-name=HCP.interim
#SBATCH --output=error/out.HCP.interim
#SBATCH --error=error/err.HCP.interim
#SBATCH --time=25:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

. ./config.sh

PILOT=15
FINAL=65
U="2.3"

OUTDIR_A=$(echo $RESDIR\hcp_adaptive_$PILOT\_$U/)
OUTDIR_P=$(echo $RESDIR\hcp_predictive_$PILOT\_$U/)

if [ ! -d "$OUTDIR_A" ]; then mkdir $OUTDIR_A; fi
if [ ! -d "$OUTDIR_P" ]; then mkdir $OUTDIR_P; fi

SIM=$i
python $SCRIPTDIR/HCP.py $U $PILOT $FINAL $SIM $FILEDIR $HCPDIR $TMPDIR $OUTDIR_A "adaptive"
python $SCRIPTDIR/HCP.py $U $PILOT $FINAL $SIM $FILEDIR $HCPDIR $TMPDIR $OUTDIR_P "predictive"
