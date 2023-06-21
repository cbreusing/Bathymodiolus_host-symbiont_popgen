#!/bin/bash
#SBATCH -J MITObim
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -o MITObim_%a.out
#SBATCH -e MITObim_%a.err
#SBATCH --array=238-257

module load perl
module load mira
module load bbmap

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filelist.txt)

reformat.sh in1=${FILE}_R1_clean.fastq in2=${FILE}_R2_clean.fastq out=${FILE}_interleaved.fastq
mkdir ${FILE}_mtDNA
cd ${FILE}_mtDNA
/gpfs/data/rbeinart/Software/MITObim/MITObim.pl -start 1 -end 30 -sample ${FILE} -ref B_septemdierum_mitogenome -readpool ../${FILE}_interleaved.fastq --quick ../B_septemdierum_mitogenome.fasta --clean






