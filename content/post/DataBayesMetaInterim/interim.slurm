#!/bin/bash
#SBATCH --job-name=interim_main
#SBATCH --mail-type=END,FAIL                      # send email if the job end or fail
#SBATCH --mail-user=keith.goldfeld@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=4:00:00                            # Time limit hrs:min:sec
#SBATCH --output=interim.out                      # Standard output and error log

module load r/3.6.3
cd /gpfs/data/troxellab/ksg/r

Rscript --vanilla interim.R