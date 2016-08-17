#!/bin/bash

#SBATCH --job-name=SIM.NAD.RFT
#SBATCH --output=error/out.SIM.NAD.RFT
#SBATCH --error=error/err.SIM.NAD.RFT
#SBATCH --time=3:00:00
#SBATCH -n 25
#SBATCH -N 1
#SBATCH -p normal

#SBATCH --mail-user=joke.durnez@gmail.com
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes


. ./config_tacc.sh

# Load modules
module load fsl/5.0.9
module load python/2.7.11

# Parameters for software
export FSLDIR=/work/01329/poldrack/lonestar/software_lonestar5/fsl-5.0.8
export PATH=$PATH:$FSLDIR
source $FSLDIR/etc/fslconf/fsl.sh
export FSLOUTPUTTYPE='NIFTI_GZ'

PILOT=10
FINAL=40
EXC="2.3"
MODEL="RFT"
BATCH=$i

simsperbatch=25
#indices=($(seq 2 226))
#indices=("${indices[@]}" $(seq 251 275))
#echo ${#indices[@]}

for k in $(seq 1 $simsperbatch);
do
  id=$(((i-1)*$simsperbatch+k))
  SEED=$id
  #SEED=${indices[$id]}

ADAPTIVE="predictive"
export OUTDIR=$(echo $RESDIR\SIM_$ADAPTIVE\_$PILOT\_$EXC\_$MODEL)
if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi
python $SCRIPTDIR/SIM.py $EXC $PILOT $FINAL $SEED $ADAPTIVE $MODEL &

done

wait
