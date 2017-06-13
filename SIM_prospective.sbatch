#!/bin/bash

#SBATCH --job-name=SIM.NAD.RFT
#SBATCH --output=error/out.SIM.NAD.RFT
#SBATCH --error=error/err.SIM.NAD.RFT
#SBATCH --qos=russpold
#SBATCH -p russpold
#SBATCH --time=24:00:00

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

export PILOT=15
export FINAL=61
export MODEL="RFT"
export SEED=$SLURM_ARRAY_TASK_ID

# exc 2.0

export EXC="2.0"
export ADAPTIVE="predictive"
export OUTDIR=$(echo $RESDIR\SIM_$ADAPTIVE\_$PILOT\_$EXC\_$MODEL)
if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi
python $SCRIPTDIR/SIM.py

# exc 2.5

export EXC="2.5"
export ADAPTIVE="predictive"
export OUTDIR=$(echo $RESDIR\SIM_$ADAPTIVE\_$PILOT\_$EXC\_$MODEL)
if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi
python $SCRIPTDIR/SIM.py

# exc 3.0

export EXC="3.0"
export ADAPTIVE="predictive"
export OUTDIR=$(echo $RESDIR\SIM_$ADAPTIVE\_$PILOT\_$EXC\_$MODEL)
if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi
python $SCRIPTDIR/SIM.py


echo "Analysis finished"
