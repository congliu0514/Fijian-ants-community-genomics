#!/bin/bash
#SBATCH --mail-user="liu.cong@oist.jp"
#SBATCH --job-name=Beast_strumigenys_fix
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mem=480G
#SBATCH -p largemem
#SBATCH -c 24
#SBATCH -n 1
#SBATCH -t 60-20     
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module use /apps/unit/EconomoU/.modulefiles
module load beast2/v2.4.4
module load Beagle/1260

beast -threads 24 beast3.xml
