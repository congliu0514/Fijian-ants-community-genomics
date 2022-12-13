#!/bin/bash
#SBATCH --mail-user="liu.cong@oist.jp"
#SBATCH --job-name=fiji_phylogeny
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mem=500G
#SBATCH -p largemem
#SBATCH -c 10
#SBATCH -n 4
#SBATCH -t 60-20
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

source activate python3
ipyrad -p params-phylogeny_outgroup.txt -s1234567 -f -b -d -c 40
