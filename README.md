# neuropower-validation

This repository contains all code to validate the methodology behind [Neuropower](www.neuropowertools.org).  

## Reproduce
To reproduce the results on a local machine, please customize config.sh with a personal folder structure.

## Workflow
To reproduce the results in [This paper](http://biorxiv.org/content/early/2016/04/20/049429), follow these pipelines:
### Example
Running **Example.sh**, produces the figures.

### Simulations
1. **SIM_prospective.sbatch** performs the simulations and exports for each simulation (1) a file with the model estimation, (2) a file with the predicted power and (3) a file with the true power.

## Computing time
These analyses should be performed on a high-performance computing environment.  The subsample procedure with the data from the [Human Connectome Project](http://www.humanconnectome.org/) take about 12 minutes per contrast per iteration (+/- 10 hours for all contrasts per iteration).  The simulations average at 9 minutes per condition (+/- 1.5 hours for all conditions per simulation).
