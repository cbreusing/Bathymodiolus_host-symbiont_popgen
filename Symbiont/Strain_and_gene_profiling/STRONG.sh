#!/bin/bash
#SBATCH -J STRONG
#SBATCH -t 96:00:00
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -p bigmem
#SBATCH --mem=1t
#SBATCH -o STRONG.out
#SBATCH -e STRONG.err

source activate STRONG

STRONG --config config.yaml --threads 32 --verbose STRONG
