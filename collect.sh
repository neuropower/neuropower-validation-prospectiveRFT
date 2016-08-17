#!/bin/bash

#SBATCH --job-name=HCP.coll
#SBATCH --output=error/out.HCP.coll
#SBATCH --error=error/err.HCP.coll
#SBATCH --time=2:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p normal

#SBATCH --mail-user=joke.durnez@gmail.com
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes

module load fsl/5.0.9
module load python/2.7.11
module load Rstats/3.2.1

export FSLDIR=/work/01329/poldrack/lonestar/software_lonestar5/fsl-5.0.8
export PATH=$PATH:$FSLDIR/bin:$AFNIDIR
source $FSLDIR/etc/fslconf/fsl.sh
export FSLOUTPUTTYPE='NIFTI_GZ'

PILOT=15
FINAL=40 #SIM
#FINAL=65 #HCP
EXC="2.3"
SIMS=500
MODALITY='SIM'
ADAPTIVE='predictive'
MODEL='RFT'

. ./config_tacc.sh

export OUTDIR=$(echo $RESDIR$MODALITY\_$ADAPTIVE\_$PILOT\_$EXC\_$MODEL)

# for i in $(seq 1 $SIMS)
# do
#   rm $OUTDIR/estimation_$MODALITY\_$i.csv
# done
#
# for i in $(seq 1 $SIMS)
# do
#   echo $i
#   prt=($(ls $OUTDIR/estimation_$MODALITY\_$i\_*.csv | sort -n -t _ -k 3))
#   for k in "${prt[@]}"
#   do
#     more $k >> $OUTDIR/estimation_$MODALITY\_$i.csv
#   done
# done

python $SCRIPTDIR/aggregate_estimation.py $PILOT $FINAL $SIMS $MODALITY $ADAPTIVE $EXC $MODEL

Rscript $HOMEDIR\Figures/SIM_figures_NIMG.R $TABDIR $HOMEDIR $FIGDIR
#Rscript $HOMEDIR\Figures/HCP_figures_NIMG.R $TABDIR $HOMEDIR $FIGDIR
