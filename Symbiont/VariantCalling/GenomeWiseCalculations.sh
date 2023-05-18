#!/bin/bash
#SBATCH -J GenomeWiseCalculations
#SBATCH -t 100:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o GenomeWiseCalculations.out
#SBATCH -e GenomeWiseCalculations.err

module load samtools
module load R
module load vcftools
module load bcftools
module load perl

# Fst and Pst calculations between individual samples
cat Thiodubiliella.Freebayes.FINAL.recode.vcf | vcf-sort -c > Thiodubiliella.Freebayes.FINAL.recode.sorted.vcf
bgzip -c Thiodubiliella.Freebayes.FINAL.recode.sorted.vcf > Thiodubiliella.Freebayes.FINAL.recode.sorted.vcf.gz
tabix -p vcf Thiodubiliella.Freebayes.FINAL.recode.sorted.vcf.gz

# Create sample-specific VCF files and coverage information
for file in `cat subset.list`
do
bcftools view -e 'REF="N"' -m2 -M2 -v snps -U -s ${file} -Ov -o ${file}.FINAL.reformat.vcf Thiodubiliella.Freebayes.FINAL.recode.sorted.vcf.gz
grep -wvf regions.txt ${file}.FINAL.reformat.vcf > tmp
mv tmp ${file}.FINAL.reformat.vcf
samtools coverage -o ${file}.coverage.txt ${file}.Thiodubiliella.subsampled.bam
sed -i "s/#rname.*//g" ${file}.coverage.txt
sed -i '1d' ${file}.coverage.txt
perl get_core_coverage.pl core_genes.txt ${file}.coverage.txt ${file}.core_coverage.txt
awk '{ total += $7; count++ } END { print total/count }' ${file}.core_coverage.txt > ${file}.expected_coverage.txt
sort -n ${file}.core_coverage.txt | awk '{ a[i++]=$7; } END { print a[int(i/2)]; }' > ${file}.expected_median_coverage.txt 
perl make_coverage_list.pl ${file}.coverage.txt ${file}.expected_median_coverage.txt ${file}.gene_coverage.list
done

# Scripts used from Picazo et al. (2019, 2022); vcf_to_mergedcounts.py needs to be slightly modified to use Freebayes VCF input:
# dp4 = info[9].split(':')
# alt_values = int(dp4[5]) 
# ref_values = int(dp4[3])

python /gpfs/data/rbeinart/cbreusing/Scripts/vcf_to_mergedcounts.py .
Rscript /gpfs/data/rbeinart/cbreusing/Scripts/structure.r Thiodubiliella.txt
python /gpfs/data/rbeinart/cbreusing/Scripts/genome_wise_calculations.py Thiodubiliella.txtFst_pos.txt Thiodubiliella_pangenome.fasta Thiodubiliella

Rscript /gpfs/data/rbeinart/cbreusing/Scripts/pangenome_structure.r Thiodubiliella.pan.txt
python /gpfs/data/rbeinart/cbreusing/Scripts/pangenome_wise_calculations.py Thiodubiliella.pan.txtFst_pos.txt Thiodubiliella_pangenome.fasta Thiodubiliella.pan

# Pst calculations between vent locations
awk '{print $1"\t"$2"\t"$3}' J2432_B238.gene_coverage.list > TC.gene_coverage.list
awk '{print $1"\t"$2"\t"$3}' BS_691.gene_coverage.list > THM.gene_coverage.list
awk '{print $1"\t"$2"\t"$3}' BS_443.gene_coverage.list > ABE.gene_coverage.list
awk '{print $1"\t"$2"\t"$3}' BS_241.gene_coverage.list > TM.gene_coverage.list

for file in `cat TC.list`
do awk '{print $4"\t"$5}' ${file}.gene_coverage.list | paste TC.gene_coverage.list - > tmp
mv tmp TC.gene_coverage.list
done

for file in `cat THM.list`
do awk '{print $4"\t"$5}' ${file}.gene_coverage.list | paste THM.gene_coverage.list - > tmp
mv tmp THM.gene_coverage.list
done

for file in `cat ABE.list`
do awk '{print $4"\t"$5}' ${file}.gene_coverage.list | paste ABE.gene_coverage.list - > tmp
mv tmp ABE.gene_coverage.list
done

for file in `cat TM.list`
do awk '{print $4"\t"$5}' ${file}.gene_coverage.list | paste TM.gene_coverage.list - > tmp
mv tmp TM.gene_coverage.list
done

perl ../../Scripts/count_cov_entries.pl TC.gene_coverage.list > TC.gene_coverage.merged.list
perl ../../Scripts/count_cov_entries.pl THM.gene_coverage.list > THM.gene_coverage.merged.list
perl ../../Scripts/count_cov_entries.pl ABE.gene_coverage.list > ABE.gene_coverage.merged.list
perl ../../Scripts/count_cov_entries.pl TM.gene_coverage.list > TM.gene_coverage.merged.list

Rscript /gpfs/data/rbeinart/cbreusing/Scripts/pangenome_structure.r Thiodubiliella.pan.merged.txt
python /gpfs/data/rbeinart/cbreusing/Scripts/pangenome_wise_calculations.py Thiodubiliella.pan.merged.txtFst_pos.txt Thiodubiliella_pangenome.fasta Thiodubiliella.pan.merged

