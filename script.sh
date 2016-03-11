# from home to sherlock

scp ~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Validation/* jdurnez@sherlock.stanford.edu:/home/jdurnez/power_peak/
scp ~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Functions/* jdurnez@sherlock.stanford.edu:/home/jdurnez/power_peak/

# from sherlock to home

scp jdurnez@sherlock.stanford.edu:/scratch/users/jdurnez/interim_biascorr.tar.gz ~/Downloads/
scp jdurnez@sherlock.stanford.edu:/scratch/users/jdurnez/power_peak_HCP.tar.gz ~/Downloads/

srun -p russpold --qos=russpold --time=10:00:00 --x11 -n1 --pty bash
srun --time=4:00:00 --x11 --pty bash

i=2
export i
sbatch HCP_interim.sbatch
sbatch SIM_interim.sbatch

for i in {1..100}
do
  export i
  sbatch SIM_interim.sbatch
done

for i in {1..40}
do
  export i
  sbatch SIM_interim_biascorr.sbatch
done

# check number of entries in result files
for j in $(seq 1 500);
do
  more estimation_hcp_$j.csv | wc
done
