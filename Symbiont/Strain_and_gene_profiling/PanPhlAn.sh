#!/bin/bash
#SBATCH -J PanPhlAn
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o PanPhlAn.out
#SBATCH -e PanPhlAn.err

module load samtools
module load bowtie2/2.3.5.1
module load python/3.7.4
module load bbmap
module load seqtk

for file in `cat subset.list`
do
samtools view -h -@ 24 -o ${file}.Thiodubiliella.subsampled.sam ${file}.Thiodubiliella.subsampled.bam
cut -f1 ${file}.Thiodubiliella.subsampled.sam | sort | uniq > ${file}_ids.lst
seqtk subseq ../${file}_R1_clean.fastq ${file}_ids.lst > ${file}_R1.fastq
seqtk subseq ../${file}_R2_clean.fastq ${file}_ids.lst > ${file}_R2.fastq
reformat.sh in=${file}_R1.fastq in2=${file}_R2.fastq out=${file}.fastq overwrite=t
panphlan_map.py -i ${file}.fastq -o map_results/${file}_mapping.csv --nproc 24 --indexes panphlan/Thiodubiliella -p panphlan/Thiodubiliella_pangenome.tsv
done

panphlan_profiling.py -p panphlan/Thiodubiliella_pangenome.tsv -v -i map_results --o_matrix gene_presence_absence.csv --min_coverage 1 --left_max 2.05 --right_min 0.30


