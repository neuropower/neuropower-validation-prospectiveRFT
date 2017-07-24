#!/bin/bash

# This file contains the directory structure used for the
# validations. To re-run all analyses, please adjust this file.

# This is the github folder with all analysis files
export HOMEDIR="/home/jdurnez/power_revision/"

# Subdirectories of home
export SCRIPTDIR=$(echo $HOMEDIR\PythonScripts/)
export HCPFILEDIR=$(echo $HOMEDIR\HcpFiles/)
export SIMFILEDIR=$(echo $HOMEDIR\SimulationFiles/)

# SCRATCHDIR --> for large i/o and large files, ideally HCP-environment
export SCRATCHDIR="/scratch/users/jdurnez/power_revision/"

# Subdirectories of scratch
export RESDIR=$(echo $SCRATCHDIR\Results/)
export PEAKDIR=$(echo $SCRATCHDIR\Peaks/)
export TMPDIR=$(echo $SCRATCHDIR\tmp/)
export TABDIR=$(echo $SCRATCHDIR\tables/)
export FIGDIR=$(echo $SCRATCHDIR\figures/)

# Create subdirectories if they don't exist already
if [ ! -d "$RESDIR" ]; then mkdir $RESDIR; fi
if [ ! -d "$TMPDIR" ]; then mkdir $TMPDIR; fi
if [ ! -d "$TABDIR" ]; then mkdir $TABDIR; fi
if [ ! -d "$FIGDIR" ]; then mkdir $FIGDIR; fi
if [ ! -d "$PEAKDIR" ]; then mkdir $PEAKDIR; fi

# Connectome-in-a-box folder
export HCPDIR="/scratch/PI/russpold/data/HCP/"
