#!/bin/bash
#SBATCH -J Freebayes
#SBATCH -t 100:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o Freebayes.out
#SBATCH -e Freebayes.err

source activate lofreq
module load bowtie2/2.3.5.1
module load samtools
module load gatk/4.1.6.0
module load picard-tools/2.17.11
module load R
module load vcftools
module load bcftools

bowtie2-build Thiodubiliella_pangenome.fasta Thiodubiliella_pangenome.fasta
samtools faidx Thiodubiliella_pangenome.fasta

for file in `cat filelist.txt`
do 
# Map reads to reference
bowtie2 --very-sensitive-local -p 32 -x Thiodubiliella_pangenome.fasta -1 ../${file}_R1_clean.fastq -2 ../${file}_R2_clean.fastq | samtools view -bS -h -@ 32 - | samtools sort -@ 32 - > ${file}.Thiodubiliella.sorted.bam
# Mark and remove duplicates
java -Xmx8g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar MarkDuplicates I=${file}.Thiodubiliella.sorted.bam O=${file}.Thiodubiliella.dedup.bam M=${file}.Thiodubiliella.metrics.txt REMOVE_DUPLICATES=true
samtools index ${file}.Thiodubiliella.dedup.bam
# Indel realignment
lofreq viterbi -f Thiodubiliella_pangenome.fasta -k ${file}.Thiodubiliella.dedup.bam | samtools sort -@ 32 - > ${file}.Thiodubiliella.realigned.bam
# Base recalibration
lofreq indelqual -f Thiodubiliella_pangenome.fasta --dindel -o ${file}.Thiodubiliella.indelqual.bam ${file}.Thiodubiliella.realigned.bam
samtools index ${file}.Thiodubiliella.indelqual.bam
samtools view -bS -h -F4 ${file}.Thiodubiliella.indelqual.bam > ${file}.Thiodubiliella.filt.dedup.bam
# Add read groups to BAM files
samtools index ${file}.Thiodubiliella.filt.dedup.bam
java -Xmx8g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar AddOrReplaceReadGroups I=${file}.Thiodubiliella.filt.dedup.bam O=${file}.Thiodubiliella.filt.dedup.RG.bam RGID=${file} RGLB=LIB_${file} RGPL=ILLUMINA RGPU=FLOWCELL1 RGSM=${file} VALIDATION_STRINGENCY=SILENT
samtools index ${file}.Thiodubiliella.filt.dedup.RG.bam
# Softclip reads - this is necessary to overcome some bugs in Freebayes that occur when reads overlap the alignment ends
bam trimBam ${file}.Thiodubiliella.filt.dedup.RG.bam - -c | samtools sort -@ 32 -n - > ${file}.Thiodubiliella.notNormalized.softclipped.bam
samtools fixmate ${file}.Thiodubiliella.notNormalized.softclipped.bam - | samtools sort -@ 32 - > ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam
# Check number of alignments for each BAM file, select the lowest number of alignments for downsampling
samtools view -c ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam >> alignment_counts.txt
done

num=869296

for file in `cat filelist.txt`
do
count=`samtools view -c ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam`
if [ $num -le $count ]
    then
    frac=`bc -l <<< $num/$count`
    samtools view -@ 32 -h -bs $frac ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam > ${file}.Thiodubiliella.subsampled.bam
fi
samtools index ${file}.Thiodubiliella.subsampled.bam
done

ls *.subsampled.bam > bam.list

freebayes -f Thiodubiliella_pangenome.fasta -L bam.list -v Thiodubiliella.Freebayes.vcf -F 0.01 -C 1 -p 1 --pooled-continuous -g 1000 -m 30 -q 20 --min-coverage 10 --haplotype-length 0 --report-monomorphic

cat Thiodubiliella.Freebayes.vcf | vcf-sort -c > Thiodubiliella.Freebayes.sorted.vcf
bgzip -c Thiodubiliella.Freebayes.sorted.vcf > Thiodubiliella.Freebayes.vcf.gz
tabix -p vcf Thiodubiliella.Freebayes.vcf.gz
bcftools filter -g 5 -i 'REF!="N" && SRP > 5 && SAP > 5 && EPP > 5 && QUAL > 20 && FORMAT/DP > 5' -o Thiodubiliella.Freebayes.filtered.vcf -O v Thiodubiliella.Freebayes.vcf.gz
vcftools --vcf Thiodubiliella.Freebayes.filtered.vcf --max-missing 0.75 --remove-indv J2432_B240 --remove-indv J2432_B243 --remove-indv J2432_B246 --remove-indv J2432_B253 --recode --recode-INFO-all --out Thiodubiliella.Freebayes
vcftools --vcf Thiodubiliella.Freebayes.recode.vcf --missing-indv --out Thiodubiliella.Freebayes.recode
vcftools --vcf Thiodubiliella.Freebayes.recode.vcf --site-mean-depth --out Thiodubiliella.Freebayes.recode
vcftools --vcf Thiodubiliella.Freebayes.recode.vcf --max-meanDP 10 --recode --recode-INFO-all --out Thiodubiliella.Freebayes.FINAL
vcftools --vcf Thiodubiliella.Freebayes.FINAL.recode.vcf --extract-FORMAT-info GT --out Thiodubiliella.Freebayes.FINAL
gatk VariantsToTable -V Thiodubiliella.Freebayes.FINAL.recode.vcf -O Thiodubiliella.Freebayes.FINAL.RD.FORMAT -F CHROM -F POS -F REF -ASGF RO
gatk VariantsToTable -V Thiodubiliella.Freebayes.FINAL.recode.vcf -O Thiodubiliella.Freebayes.FINAL.AD.FORMAT -F CHROM -F POS -F ALT -ASGF AO --split-multi-allelic

