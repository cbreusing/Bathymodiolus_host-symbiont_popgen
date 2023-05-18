#!/bin/bash
#SBATCH -J mtVariantCalling
#SBATCH -t 24:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o mtVariantCalling_J2432_B%a.out
#SBATCH -e mtVariantCalling_J2432_B%a.err
#SBATCH --array=238-257

source activate lofreq
module load bowtie2/2.3.5.1
module load samtools
module load gatk/4.1.6.0
module load picard-tools/2.17.11
module load R
module load vcftools
module load bcftools

# Examples are shown based on Tow Cam population J2432
bowtie2-build J2432_B${SLURM_ARRAY_TASK_ID}_mitogenome.fasta J2432_B${SLURM_ARRAY_TASK_ID}_mitogenome.fasta
samtools faidx J2432_B${SLURM_ARRAY_TASK_ID}_mitogenome.fasta

bowtie2 --very-sensitive -p 24 -x J2432_B${SLURM_ARRAY_TASK_ID}_mitogenome.fasta -1 ../J2432_B${SLURM_ARRAY_TASK_ID}_R1_clean.fastq -2 ../J2432_B${SLURM_ARRAY_TASK_ID}_R2_clean.fastq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > J2432_B${SLURM_ARRAY_TASK_ID}.mt.sorted.bam
java -Xmx2g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar MarkDuplicates I=J2432_B${SLURM_ARRAY_TASK_ID}.mt.sorted.bam O=J2432_B${SLURM_ARRAY_TASK_ID}.mt.dedup.bam M=J2432_B${SLURM_ARRAY_TASK_ID}.mt.metrics.txt REMOVE_DUPLICATES=true
samtools index J2432_B${SLURM_ARRAY_TASK_ID}.mt.dedup.bam
lofreq viterbi -f J2432_B${SLURM_ARRAY_TASK_ID}_mitogenome.fasta -k J2432_B${SLURM_ARRAY_TASK_ID}.mt.dedup.bam | samtools sort -@ 24 - > J2432_B${SLURM_ARRAY_TASK_ID}.mt.realigned.bam
lofreq indelqual -f J2432_B${SLURM_ARRAY_TASK_ID}_mitogenome.fasta --dindel -o J2432_B${SLURM_ARRAY_TASK_ID}.mt.indelqual.bam J2432_B${SLURM_ARRAY_TASK_ID}.mt.realigned.bam
samtools index J2432_B${SLURM_ARRAY_TASK_ID}.mt.indelqual.bam
samtools view -bS -h -F4 J2432_B${SLURM_ARRAY_TASK_ID}.mt.indelqual.bam > J2432_B${SLURM_ARRAY_TASK_ID}.filt.mt.dedup.bam
java -Xmx2g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar AddOrReplaceReadGroups I=J2432_B${SLURM_ARRAY_TASK_ID}.filt.mt.dedup.bam O=J2432_B${SLURM_ARRAY_TASK_ID}.filt.mt.dedup.RG.bam RGID=J2432_B${SLURM_ARRAY_TASK_ID} RGLB=LIB_J2432_B${SLURM_ARRAY_TASK_ID} RGPL=ILLUMINA RGPU=FLOWCELL1 RGSM=J2432_B${SLURM_ARRAY_TASK_ID} VALIDATION_STRINGENCY=SILENT

lofreq call --call-indels -q 20 -Q 20 -m 30 -f J2432_B${SLURM_ARRAY_TASK_ID}_mitogenome.fasta -o J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.vcf J2432_B${SLURM_ARRAY_TASK_ID}.filt.mt.dedup.RG.bam

cat J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.vcf | vcf-sort -c > J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.sorted.vcf
bgzip -c J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.sorted.vcf > J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.vcf.gz
tabix -p vcf J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.vcf.gz
bcftools filter -g 5 -o J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.filtered.vcf -O v J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.vcf.gz
vcftools --vcf J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.filtered.vcf --max-meanDP 1000 --recode --recode-INFO-all --out J2432_B${SLURM_ARRAY_TASK_ID}.mt.LoFreq.FINAL

