srun -p russpold --qos=russpold --time=10:00:00 --x11 -n1 --pty bash

for i in {1..200}
do
  export i
  sbatch SIM_prospective.sbatch
done

for i in {1..100}
do
  export i
  sbatch HCP_interim.sbatch
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
