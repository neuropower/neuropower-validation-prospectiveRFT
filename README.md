# neuropower-validation

This repository contains all code to validate the methodology behind [Neuropower](www.neuropowertools.org).  

### Workflow

##### Reproduce
To reproduce the results on a local machine, please customize config.sh with a personal folder structure.
To reproduce the results in [this paper](http://biorxiv.org/content/early/2016/04/20/049429), follow these pipelines:

##### Example
Running **Example.sh**, produces the figures.

##### Simulations
1. **SIM_prospective.sbatch** performs the simulations and exports for each simulation (1) a file with the model estimation, (2) a file with the predicted power and (3) a file with the true power.

### Computing time
These analyses should be performed on a high-performance computing environment.  We performed the analysis on [Stanford's Sherlock](http://sherlock.stanford.edu/mediawiki/index.php/Main_Page): 120 general compute nodes with dual socket Intel(R) Xeon(R) CPU E5-2650 v2 @ 2.60GHz (8 core/socket); 64 GB 1866 MHz DDR3.  The subsample procedure with the data from the [Human Connectome Project](http://www.humanconnectome.org/) take about 12 minutes per contrast per iteration (+/- 10 hours for all contrasts per iteration).  The simulations average at 9 minutes per condition (+/- 1.5 hours for all conditions per simulation).
