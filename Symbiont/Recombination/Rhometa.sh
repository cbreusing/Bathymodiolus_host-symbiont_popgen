#!/bin/bash
#SBATCH -J rhometa
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o rhometa.out
#SBATCH -e rhometa.err

# Example is shown for Rhometa analyses across populations
source activate rhometa

nextflow run /gpfs/data/rbeinart/Software/rhometa/theta_est.nf --bam ../Thiodubiliella.merged.norm.bam --fa ../Thiodubiliella_pangenome.fasta
nextflow run /gpfs/data/rbeinart/Software/rhometa/lookup_table_gen.nf --theta 0.00021508486573214717 
nextflow run /gpfs/data/rbeinart/Software/rhometa/rho_est.nf --bam ../Thiodubiliella.merged.norm.bam --fa ../Thiodubiliella_pangenome.fasta
