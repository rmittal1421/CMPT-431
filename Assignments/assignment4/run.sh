#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --mem=10G

srun ./triangle_counting_parallel --strategy 1
