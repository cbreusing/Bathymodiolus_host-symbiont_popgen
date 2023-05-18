#!/bin/bash
#SBATCH -J MITObim
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -o MITObim_J2432_B%a.out
#SBATCH -e MITObim_J2432_B%a.err
#SBATCH --array=238-257

module load perl
module load mira
module load bbmap

# Examples are shown based on Tow Cam population J2432
reformat.sh in1=J2432_B${SLURM_ARRAY_TASK_ID}_R1_clean.fastq in2=J2432_B${SLURM_ARRAY_TASK_ID}_R2_clean.fastq out=J2432_B${SLURM_ARRAY_TASK_ID}_interleaved.fastq
mkdir J2432_B${SLURM_ARRAY_TASK_ID}_mtDNA
cd J2432_B${SLURM_ARRAY_TASK_ID}_mtDNA
/gpfs/data/rbeinart/Software/MITObim/MITObim.pl -start 1 -end 30 -sample J2432_B${SLURM_ARRAY_TASK_ID} -ref B_septemdierum_mitogenome -readpool ../J2432_B${SLURM_ARRAY_TASK_ID}_interleaved.fastq --quick ../B_septemdierum_mitogenome.fasta --clean






