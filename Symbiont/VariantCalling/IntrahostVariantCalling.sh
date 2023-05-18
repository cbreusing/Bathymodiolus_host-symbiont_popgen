#!/bin/bash
#SBATCH -J Freebayes
#SBATCH -t 100:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o Freebayes.out
#SBATCH -e Freebayes.err

module load bowtie2/2.3.5.1
module load samtools
module load gatk/4.1.6.0
module load picard-tools/2.17.11
module load R
module load vcftools
module load bcftools

for file in `cat subset.list`
do

freebayes -f Thiodubiliella_pangenome.fasta -b ${file}.Thiodubiliella.subsampled.bam -v ${file}.normalized.vcf -F 0.00 -C 1 -p 1 --pooled-continuous -g 1000 -m 30 -q 20 --min-coverage 10 --haplotype-length 0 --report-monomorphic

cat ${file}.normalized.vcf | vcf-sort -c > ${file}.normalized.sorted.vcf
bgzip -c ${file}.normalized.sorted.vcf > ${file}.normalized.vcf.gz
tabix -p vcf ${file}.normalized.vcf.gz
bcftools filter -g 5 -o ${file}.normalized.filtered.vcf -O v ${file}.normalized.vcf.gz
gatk VariantsToTable -V ${file}.normalized.filtered.vcf -O ${file}.normalized.notFiltered.AF.txt -ASGF AD

done
