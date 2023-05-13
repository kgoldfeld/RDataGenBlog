#!/bin/bash
#SBATCH --job-name=ss
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=keith.goldfeld@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=12:00:00                              # Time limit hrs:min:sec
#SBATCH --output=ss.out                              # Standard output and error log

module load r/4.2.2

cd /gpfs/data/troxellab/ksg/r
Rscript --vanilla sample_size.R

