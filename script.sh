# from home to sherlock

scp ~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Validation/* jdurnez@sherlock.stanford.edu:/home/jdurnez/power_peak/
scp ~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Functions/* jdurnez@sherlock.stanford.edu:/home/jdurnez/power_peak/
scp ~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Validation/SIM_interim_biascorr.py jdurnez@sherlock.stanford.edu:/home/jdurnez/power_peak/

# from sherlock to home

scp jdurnez@sherlock.stanford.edu:/scratch/users/jdurnez/power_SIM.tar.gz ~/Downloads/
scp jdurnez@sherlock.stanford.edu:/scratch/users/jdurnez/interim.tar.gz ~/Downloads/

srun -p russpold --qos=russpold --time=10:00:00 --x11 -n1 --pty bash
srun --time=4:00:00 --x11 --pty bash

i=1
export i
sbatch SIM_interim.sbatch
sbatch HCP_interim.sbatch

for i in {251..500}
do
  export i
  sbatch HCP_interim.sbatch
done

for i in {1..20}
do
  export i
  sbatch SIM_interim_biascorr.sbatch
done

for i in {260..500}
do
  export i
  sbatch HCP_interim.sbatch
done


# check number of entries in result files
for j in $(seq 1 500);
do
  more estimation_hcp_$j.csv | wc
done



sim_adaptive / u2 : 750
sim_adaptive / u3 : 1196
sim_nonadaptive / u2 : 750
sim_nonadaptive / u3 : 1186
