srun -p russpold --qos=russpold --time=10:00:00 --x11 -n1 --pty bash

for i in {1..200}
do
  export i
  sbatch SIM_prospective.sbatch
done


i=1
export i
sbatch SIM_prospective.sbatch


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

rsync -azP jdurnez@sherlock.stanford.edu:/scratch/users/jdurnez/power/tables/ /Users/Joke/Documents/Onderzoek/ProjectsOngoing/power/ValidationResults/
