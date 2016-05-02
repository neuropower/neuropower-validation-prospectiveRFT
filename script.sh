srun -p russpold --qos=russpold --time=10:00:00 --x11 -n1 --pty bash

i=1
export i
sbatch HCP_interim.sbatch
sbatch HCP_prospective.sbatch
sbatch SIM_interim.sbatch
sbatch SIM.sbatch

for i in {1..250}
do
  export i
  sbatch SIM.sbatch
done

for i in {1..125}
do
  export i
  sbatch conditional.sbatch
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
