#!/bin/bash
#SBATCH --job-name=sim
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=keith.goldfeld@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=12:00:00                              # Time limit hrs:min:sec
#SBATCH --output=simulation.out                     # Standard output and error log

module load r/4.1.1

cd /gpfs/data/troxellab/ksg/r
Rscript --vanilla simulation.R
