#!/bin/bash
#SBATCH -J ReadFilter
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o ReadFilter_J2432_B%a.out
#SBATCH -e ReadFilter_J2432_B%a.err
#SBATCH --array=238-257

module load samtools
module load seqtk
module load bowtie2/2.3.5.1
module load bbmap

# Examples are shown based on Tow Cam population J2432
java -jar /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/trinity-plugins/Trimmomatic/trimmomatic-0.36.jar PE -threads 24 -phred33 J2432_B${SLURM_ARRAY_TASK_ID}_1.fq.gz J2432_B${SLURM_ARRAY_TASK_ID}_2.fq.gz J2432_B${SLURM_ARRAY_TASK_ID}_R1_paired.fq J2432_B${SLURM_ARRAY_TASK_ID}_R1_unpaired.fq J2432_B${SLURM_ARRAY_TASK_ID}_R2_paired.fq J2432_B${SLURM_ARRAY_TASK_ID}_R2_unpaired.fq ILLUMINACLIP:/gpfs/data/rbeinart/cbreusing/Adapters/AllAdaptors.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:75
bowtie2 -p 24 -x /gpfs/data/rbeinart/Databases/contaminants -1 J2432_B${SLURM_ARRAY_TASK_ID}_R1_paired.fq -2 J2432_B${SLURM_ARRAY_TASK_ID}_R2_paired.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > J2432_B${SLURM_ARRAY_TASK_ID}.bowtie2.cont.sorted.bam
samtools view -@ 24 -f12 J2432_B${SLURM_ARRAY_TASK_ID}.bowtie2.cont.sorted.bam > J2432_B${SLURM_ARRAY_TASK_ID}.cont.unmapped.sam
cut -f1 J2432_B${SLURM_ARRAY_TASK_ID}.cont.unmapped.sam | sort | uniq > J2432_B${SLURM_ARRAY_TASK_ID}.cont.unmapped_ids.lst
seqtk subseq J2432_B${SLURM_ARRAY_TASK_ID}_R1_paired.fq J2432_B${SLURM_ARRAY_TASK_ID}.cont.unmapped_ids.lst > J2432_B${SLURM_ARRAY_TASK_ID}_R1_clean.fastq
seqtk subseq J2432_B${SLURM_ARRAY_TASK_ID}_R2_paired.fq J2432_B${SLURM_ARRAY_TASK_ID}.cont.unmapped_ids.lst > J2432_B${SLURM_ARRAY_TASK_ID}_R2_clean.fastq



