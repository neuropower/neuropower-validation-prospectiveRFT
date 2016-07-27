#!/bin/bash

#SBATCH -J HCP.NAD.RFT
#SBATCH -o error/out.HCP.NAD.RFT
#SBATCH -e error/err.HCP.NAD.RFT
#SBATCH -t 30:00:00
#SBATCH -N 1
#SBATCH -n 25
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

# Parameters for analysis type
PILOT=15
FINAL=65
EXC="2.3"
MODEL="RFT"
ADAPTIVE="predictive"

export SCRIPTDIR=$(echo $HOMEDIR\PythonScripts/)
export FILEDIR=$(echo $HOMEDIR\HcpFiles/)
export OUTDIR=$(echo $RESDIR\HCP_$ADAPTIVE\_$PILOT\_$EXC\_$MODEL/)

if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi

# How many jobs are there for each script
JOBSPERSCRIPT=5

# Loop over jobs, loop over contrasts
for job in $(seq 1 $JOBSPERSCRIPT)
do
  BASE=$(($i-1))
  SEED=$(($BASE*$JOBSPERSCRIPT+$job))

  starts=(0 10 20 30 40)
  ends=(10 20 30 40 47)

    for start in {0..5};
    do
      echo ${starts[start]}
      echo ${ends[start]}
    echo /usr/bin/time -v python $SCRIPTDIR/HCP.py $EXC $PILOT $FINAL $SEED $ADAPTIVE $MODEL ${starts[start]} ${ends[start]} &
    sleep 2m

    done

done

wait
