module load julia/1.2.0
sbatch slurm-pgcm-ridge-0.0_and_0.2 equil 10 100
sbatch slurm-pgcm-ridge-0.4_and_0.8 equil 10 100
sbatch slurm-pgcm-ridge-0.6_and_a-sens equil 10 100
sbatch slurm-pgcm-ridge-0.6_r-sens equil 10 100

