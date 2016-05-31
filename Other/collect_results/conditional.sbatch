#!/bin/bash

#SBATCH --job-name=HCP_power
#SBATCH --output=out/out.cond.HCP
#SBATCH --error=error/err.cond.HCP
#SBATCH --time=8:00:00
#SBATCH --partition=russpold
#SBATCH --qos=russpold

BASEFOLDER="/scratch/users/jdurnez/"
OUTDIR="/scratch/users/jdurnez/interim_results/cond_summ/"

l=$i
steps=2
simstart=$((($l-1)*$steps+1))
simstop=$((simstart+$steps-1))
python estimation_conditional_hcp.py 15 61 $simstart $simstop 16 'hcp' 'adaptive' 'u2' $BASEFOLDER $OUTDIR
#python estimation_conditional_hcp.py 15 61 $simstart $simstop 16 'hcp' 'adaptive' 'u3' $BASEFOLDER $OUTDIR
#python estimation_conditional_hcp.py 15 61 $simstart $simstop 16 'hcp' 'nonadaptive' 'u2' $BASEFOLDER $OUTDIR
#python estimation_conditional_hcp.py 15 61 $simstart $simstop 16 'hcp' 'nonadaptive' 'u3' $BASEFOLDER $OUTDIR
