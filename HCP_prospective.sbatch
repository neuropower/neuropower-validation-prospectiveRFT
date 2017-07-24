#!/bin/bash

#SBATCH --output=error/POWER.%a.txt
#SBATCH --error=error/POWER.%a.txt
##SBATCH --qos=russpold
##SBATCH -p russpold
#SBATCH --time=35:00:00
#SBATCH --mail-user=joke.durnez@gmail.com
#SBATCH --mail-type=ALL   # email me when the job starts

. ./config_tacc.sh

# Load modules
module use /scratch/PI/russpold/modules
source /share/PI/russpold/software/setup_all.sh

# Parameters for software
export PATH=$PATH:$FSLDIR
source $FSLDIR/etc/fslconf/fsl.sh
export FSLOUTPUTTYPE='NIFTI_GZ'

# Parameters for analysis type
PILOT=15
FINAL=61
EXC="2.5"
MODEL="RFT"
ADAPTIVE="predictive"

export SCRIPTDIR=$(echo $HOMEDIR\PythonScripts/)
export FILEDIR=$(echo $HOMEDIR\HcpFiles/)
export OUTDIR=$(echo $RESDIR\HCP_$ADAPTIVE\_$PILOT\_$EXC\_$MODEL/)

if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi

export SEED=$SLURM_ARRAY_TASK_ID

# Loop over jobs, loop over contrasts
starts=(0)
ends=(47)

start=0
TEMPORARY=$(echo $TMPDIR\tempdir\_$SEED\_$start)
mkdir $TEMPORARY
python $SCRIPTDIR/HCP.py $EXC $PILOT $FINAL $SEED $ADAPTIVE $MODEL ${starts[start]} ${ends[start]} $TEMPORARY

echo "Analysis finished"
