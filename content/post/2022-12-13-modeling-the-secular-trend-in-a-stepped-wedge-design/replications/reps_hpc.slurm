#!/bin/bash
#SBATCH --job-name=sw
#SBATCH --partition=cpu_short
#SBATCH --time=3:00:00                              # Time limit hrs:min:sec
#SBATCH --output=sw_hpc.out                  # Standard output and error log
#SBATCH --mail-type=END,FAIL                        # send email when job ends
#SBATCH --mail-user=keith.goldfeld@nyulangone.org

module load r/4.1.1
cd /gpfs/data/troxellab/ksg/r

Rscript --vanilla reps_hpc.R
