#!/bin/bash

#SBATCH --job-name=HCP.prospective
#SBATCH --output=error/out.HCP.prospective
#SBATCH --error=error/err.HCP.prospective
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

. ./config.sh

PILOT=10
FINAL=65
U="2.3"

SCRIPTDIR=$(echo $HOMEDIR\PythonScripts/)
FILEDIR=$(echo $HOMEDIR\HcpFiles/)

OUTDIR=$(echo $RESDIR\HCP_predictive_$PILOT\_$U/)

if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi

SIM=$i
python $SCRIPTDIR/HCP.py $U $PILOT $FINAL $SIM $FILEDIR $HCPDIR $TMPDIR $OUTDIR "predictive"
