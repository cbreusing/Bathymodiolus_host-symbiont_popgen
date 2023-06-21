#!/bin/bash
#SBATCH -J mtVariantCalling
#SBATCH -t 24:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o mtVariantCalling_%a.out
#SBATCH -e mtVariantCalling_%a.err
#SBATCH --array=1-105

source activate lofreq
module load bowtie2/2.3.5.1
module load samtools
module load gatk/4.1.6.0
module load picard-tools/2.17.11
module load R
module load vcftools
module load bcftools

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filelist.txt)

bowtie2-build ${FILE}_mitogenome.fasta ${FILE}_mitogenome.fasta
samtools faidx ${FILE}_mitogenome.fasta

bowtie2 --very-sensitive -p 24 -x ${FILE}_mitogenome.fasta -1 ../${FILE}_R1_clean.fastq -2 ../${FILE}_R2_clean.fastq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${FILE}.mt.sorted.bam
java -Xmx2g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar MarkDuplicates I=${FILE}.mt.sorted.bam O=${FILE}.mt.dedup.bam M=${FILE}.mt.metrics.txt REMOVE_DUPLICATES=true
samtools index ${FILE}.mt.dedup.bam
lofreq viterbi -f ${FILE}_mitogenome.fasta -k ${FILE}.mt.dedup.bam | samtools sort -@ 24 - > ${FILE}.mt.realigned.bam
lofreq indelqual -f ${FILE}_mitogenome.fasta --dindel -o ${FILE}.mt.indelqual.bam ${FILE}.mt.realigned.bam
samtools index ${FILE}.mt.indelqual.bam
samtools view -bS -h -F4 ${FILE}.mt.indelqual.bam > ${FILE}.filt.mt.dedup.bam
java -Xmx2g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar AddOrReplaceReadGroups I=${FILE}.filt.mt.dedup.bam O=${FILE}.filt.mt.dedup.RG.bam RGID=${FILE} RGLB=LIB_${FILE} RGPL=ILLUMINA RGPU=FLOWCELL1 RGSM=${FILE} VALIDATION_STRINGENCY=SILENT

lofreq call --call-indels -q 20 -Q 20 -m 30 -f ${FILE}_mitogenome.fasta -o ${FILE}.mt.LoFreq.vcf ${FILE}.filt.mt.dedup.RG.bam

cat ${FILE}.mt.LoFreq.vcf | vcf-sort -c > ${FILE}.mt.LoFreq.sorted.vcf
bgzip -c ${FILE}.mt.LoFreq.sorted.vcf > ${FILE}.mt.LoFreq.vcf.gz
tabix -p vcf ${FILE}.mt.LoFreq.vcf.gz
bcftools filter -g 5 -o ${FILE}.mt.LoFreq.filtered.vcf -O v ${FILE}.mt.LoFreq.vcf.gz
vcftools --vcf ${FILE}.mt.LoFreq.filtered.vcf --max-meanDP 1000 --recode --recode-INFO-all --out ${FILE}.mt.LoFreq.FINAL

