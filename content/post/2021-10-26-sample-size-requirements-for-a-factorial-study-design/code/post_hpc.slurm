#!/bin/bash
#SBATCH --job-name=ppst_hpc
#SBATCH --partition=cpu_short
#SBATCH --time=12:00:00                           # Time limit hrs:min:sec
#SBATCH --mem-per-cpu=8G
#SBATCH --output=post_hpc.out                     # Standard output and error log
#SBATCH --mail-type=END,FAIL                      # send email when job ends
#SBATCH --mail-user=keith.goldfeld@nyulangone.org

module load r/4.0.3
module load cmdstan/2.25.0
cd /gpfs/data/troxellab/ksg/r

Rscript --vanilla post_hpc.R