#!/bin/bash
#SBATCH -J ReadFilter
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o ReadFilter_%a.out
#SBATCH -e ReadFilter_%a.err
#SBATCH --array=1-105

module load samtools
module load seqtk
module load bowtie2/2.3.5.1
module load bbmap

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filelist.txt)

java -jar /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/trinity-plugins/Trimmomatic/trimmomatic-0.36.jar PE -threads 24 -phred33 ${FILE}_1.fq.gz ${FILE}_2.fq.gz ${FILE}_R1_paired.fq ${FILE}_R1_unpaired.fq ${FILE}_R2_paired.fq ${FILE}_R2_unpaired.fq ILLUMINACLIP:/gpfs/data/rbeinart/cbreusing/Adapters/AllAdaptors.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:75
bowtie2 -p 24 -x /gpfs/data/rbeinart/Databases/contaminants -1 ${FILE}_R1_paired.fq -2 ${FILE}_R2_paired.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${FILE}.bowtie2.cont.sorted.bam
samtools view -@ 24 -f12 ${FILE}.bowtie2.cont.sorted.bam > ${FILE}.cont.unmapped.sam
cut -f1 ${FILE}.cont.unmapped.sam | sort | uniq > ${FILE}.cont.unmapped_ids.lst
seqtk subseq ${FILE}_R1_paired.fq ${FILE}.cont.unmapped_ids.lst > ${FILE}_R1_clean.fastq
seqtk subseq ${FILE}_R2_paired.fq ${FILE}.cont.unmapped_ids.lst > ${FILE}_R2_clean.fastq



