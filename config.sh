#!/bin/bash

# This file contains the directory structure used for the
# validations. To re-run all analyses, please adjust this file.

# This is the github folder with all analysis files
HOMEDIR="/home/jdurnez/power/"

# Subdirectories of home
SCRIPTDIR=$(echo $HOMEDIR\PythonScripts/)
HCPFILEDIR=$(echo $HOMEDIR\HcpFiles/)
SIMFILEDIR=$(echo $HOMEDIR\SimulationFiles/)

# SCRATCHDIR --> for large i/o and large files, ideally HCP-environment
SCRATCHDIR="/scratch/users/jdurnez/power/"

# Subdirectories of scratch
RESDIR=$(echo $SCRATCHDIR\Results/)
TMPDIR=$(echo $SCRATCHDIR\tmp/)
TABDIR=$(echo $SCRATCHDIR\tables/)
FIGDIR=$(echo $SCRATCHDIR\figures/)

# Connectome-in-a-box folder
HCPDIR="/scratch/PI/russpold/data/HCP/"
