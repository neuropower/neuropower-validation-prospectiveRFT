# from home to sherlock

scp ~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Validation/* jdurnez@sherlock.stanford.edu:/home/jdurnez/power_peak/
scp ~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Functions/* jdurnez@sherlock.stanford.edu:/home/jdurnez/power_peak/
scp ~/Documents/Onderzoek/Studie_4_propow/ProspectivePower-Validation/SIM_interim_biascorr.py jdurnez@sherlock.stanford.edu:/home/jdurnez/power_peak/

scp ~/Documents/Onderzoek/Studie_meta_analysis/EffectSizeEstimation/* jdurnez@sherlock.stanford.edu:/home/jdurnez/
scp jdurnez@sherlock.stanford.edu:/home/jdurnez/clusterstability* ~/Downloads/
scp jdurnez@sherlock.stanford.edu:/scratch/users/jdurnez/clusterstab/minss.nii.gz ~/Downloads/
# from sherlock to home

scp jdurnez@sherlock.stanford.edu:/scratch/users/jdurnez/power.tar.gz ~/Downloads/
scp jdurnez@sherlock.stanford.edu:/scratch/users/jdurnez/interim.tar.gz ~/Downloads/

srun -p russpold --qos=russpold --time=10:00:00 --x11 -n1 --pty bash
srun --time=4:00:00 --x11 --pty bash

i=1
export i
sbatch SIM_interim.sbatch
sbatch HCP_interim.sbatch

for i in {1..250}
do
  export i
  sbatch SIM.sbatch
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
hcp_adaptive / u2 : ?!
hcp_adaptive / u3 : 500
hcp_nonadaptive / u2: 250
hcp_nonadaptive / u3: 500


8.41 --> 8.51
10 min * 16 = 2u30
2u30 * 4 = 10u
