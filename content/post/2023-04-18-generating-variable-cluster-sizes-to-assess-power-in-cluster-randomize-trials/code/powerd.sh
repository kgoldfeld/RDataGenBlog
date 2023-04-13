#!/bin/bash
#SBATCH --job-name=powerd
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=keith.goldfeld@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --mem=32G
#SBATCH --time=12:00:00                              # Time limit hrs:min:sec
#SBATCH --output=powerd.out                            # Standard output and error log

module load r/4.2.2

cd /gpfs/data/troxellab/ksg/r
Rscript --vanilla powerd.R
