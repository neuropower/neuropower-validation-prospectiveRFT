# neuropower-validation

This repository contains all code to produce the methods and results in our paper _'Power and sample size calculations for fMRI studies based on the prevalence of active peaks'_.  Everything was run on a SLURM system.  The structure of this repository is as follows:

Folders:
- **Example/** contains all code used for the example in the paper
- **PythonScripts/**: contains a few pythonscripts for the HCP or simulated validations
- **SimulationFiles/**: contains files used by the simulations (activation maps,...)
- **HcpFiles/**: contains files used by the HCP validations (userID's,...)
- **Functions/**: python functions called by the python scripts
- **Figures/**: contains scripts producing the figures in the paper

Scripts:
- **SIM_tacc.sbatch**: the batch script that executes all simulations
- **HCP-prospective.sbatch**: the batch script that executes all HCP validations
- **collect.sh**: a batch script (meant for interactive use) that combines all data
- **config_tacc.sh**: some parameters about local machine structure, see below


### Workflow

##### Reproduce
To reproduce the results on a local machine, please customize **config_tacc.sh** with a personal folder structure.
To reproduce the results in [this paper](http://biorxiv.org/content/early/2016/04/20/049429), follow these pipelines:

##### Example
All files are inside folder **Example/**
Running **Example.sh** calls **Example.py**, while using **Mask.nii.gz** and **Zstat1.nii.gz**, which produces the figures.

##### Validations
1. **SIM_prospective.sbatch** performs the simulations and exports for each simulation (1) a file with the model estimation, (2) a file with the predicted power and (3) a file with the true power.  Calls **PythonScripts/SIM.py** while using files in the folder **SimulationFiles/**.  Should be run 20 times with environment variable `i` ranging from 0 to 20.
2. **HCP_prospective.sbatch** performs the resampling scheme using the [Human Connectome Project](http://www.humanconnectome.org/) data.  Calls **PythonScripts/HCP.py** while using files in the folder **HcpFiles/**.  Should be run 25 times with environment variable `i` ranging from 0 to 20.
3. **collect.sh**: depending on whether the parameter `MODALITY` is set to `HCP` or `SIM`, this script combines and reshapes all the data in 4 files:
  - true_$MODALITY_predictive_2.3_RFT.csv: true values for power
  - prediction_$MODALITY_predictive_2.3_RFT.csv: predicted values for power
  - estimation_$MODALITY_predictive_2.3_RFT.csv: model estimation parameters
  - conditional_$MODALITY_predictive_2.3_RFT.csv: predicted sample size versus true sample size for 0.8 power
  This script calls **PythonScripts/aggregate_estimation.py**

### Computing time
These analyses should be performed on a high-performance computing environment.  We performed the analysis on [Texas Advanced Computing Center](https://www.tacc.utexas.edu/): 120 general compute nodes with dual socket Intel(R) Xeon(R) CPU E5-2650 v2 @ 2.60GHz (8 core/socket); 64 GB 1866 MHz DDR3.  
- HCP analyses: computing information for 1 resample, 1 contrast:
```
User time (seconds): 557.53
System time (seconds): 52.56
Percent of CPU this job got: 84%
Elapsed (wall clock) time (h:mm:ss or m:ss): 12:04.25
Maximum resident set size (kbytes): 1821596
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 150
Minor (reclaiming a frame) page faults: 21774042
Voluntary context switches: 2223540
Involuntary context switches: 1269
File system inputs: 8738958
File system outputs: 7781542
Page size (bytes): 4096
Exit status: 0
```
- SIM analyses: computing information for 1 simulation, 1 condition:
```
User time (seconds): 199.64
System time (seconds): 7.72
Percent of CPU this job got: 84%
Elapsed (wall clock) time (h:mm:ss or m:ss): 4:04.51
Maximum resident set size (kbytes): 408764
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 43
Minor (reclaiming a frame) page faults: 1583111
Voluntary context switches: 75497
Involuntary context switches: 543
File system inputs: 3132563
File system outputs: 4279842
Page size (bytes): 4096
Exit status: 0
```
