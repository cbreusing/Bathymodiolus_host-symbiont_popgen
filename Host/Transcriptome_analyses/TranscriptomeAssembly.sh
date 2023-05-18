#!/bin/bash
#SBATCH -J TranscriptomeAssembly
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=150g
#SBATCH -o TranscriptomeAssembly.out
#SBATCH -e TranscriptomeAssembly.err

module load fastqc
module load perl
module load samtools
module load bowtie2/2.3.5.1
module load bbmap
module load jellyfish/2.2.10
module load python/3.7.4
module load cdhit
module load transdecoder/5.4.0
module load blast/2.6.0+ hmmer/3.1b2
module load seqtk

# This is an initial quality check
fastqc -t 24 /gpfs/data/rbeinart/Shared/Raw_reads/Bathymodiolus_septemdierum_RNAseq_reads/*.fastq.gz
 
for file in R1928-5_241 R1931-6_461 R1935-3_710
do

# Quality and adaptor trimming
java -jar /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/trinity-plugins/Trimmomatic/trimmomatic-0.36.jar PE -threads 24 -phred33 /gpfs/data/rbeinart/Shared/Raw_reads/Bathymodiolus_septemdierum_RNAseq_reads/${file}_1.fastq.gz /gpfs/data/rbeinart/Shared/Raw_reads/Bathymodiolus_septemdierum_RNAseq_reads/${file}_2.fastq.gz ${file}_R1_paired.fastq ${file}_R1_unpaired.fastq ${file}_R2_paired.fastq ${file}_R2_unpaired.fastq ILLUMINACLIP:/gpfs/data/rbeinart/cbreusing/Adapters/Illumina.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:75

# These commands remove contaminants from the data, such as human DNA and the PhiX standard
bowtie2 -p 24 -x /gpfs/data/rbeinart/Databases/contaminants -1 ${file}_R1_paired.fastq -2 ${file}_R2_paired.fastq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.bowtie2.sorted.bam
samtools view -@ 24 -f12 ${file}.bowtie2.sorted.bam > ${file}.unmapped.sam
cut -f1 ${file}.unmapped.sam | sort | uniq > ${file}.unmapped_ids.lst
seqtk subseq ${file}_R1_paired.fastq ${file}.unmapped_ids.lst > ${file}_R1_clean.fastq
seqtk subseq ${file}_R2_paired.fastq ${file}.unmapped_ids.lst > ${file}_R2_clean.fastq

# The contaminant free samples are then merged and used in sortmerna to remove rRNA
/gpfs/data/rbeinart/Software/sortmerna-2.1b/scripts/merge-paired-reads.sh ${file}_R1_clean.fastq ${file}_R2_clean.fastq ${file}_merged.fastq

sortmerna --ref /gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-bac-16s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-bac-23s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-arc-16s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-arc-23s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-euk-18s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-euk-28s:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/rfam-5s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/rfam-5.8s-db --reads ${file}_merged.fastq --sam --best 1 --min_lis 2 --fastx --paired_in --aligned ${file}_rRNA --other ${file}_non_rRNA --log -a 24 -v 

# Separate host and symbiont reads with BBmap
bbmap.sh -Xmx50g ref=/gpfs/data/rbeinart/cbreusing/Bathymodiolus_septemdierum_metagenomics/GCA_013416635.1_ASM1341663v1_genomic.fna nodisk in=${file}_non_rRNA.fastq ambiguous=best outm=${file}_symbiont_#.fq outu=${file}_host_#.fq

done

# Provide a sample list of separated host reads
for file in R1928-5_241_host R1931-6_461_host R1935-3_710_host
do

# These commands correct sequencing errors and remove flagged overrepresented sequences
perl /gpfs/data/rbeinart/Software/rcorrector/run_rcorrector.pl -1 ${file}_1.fq -2 ${file}_2.fq -t 24
python /gpfs/data/rbeinart/Software/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 ${file}_1.cor.fq -2 ${file}_2.cor.fq -s ${file}
fastqc -t 24 unfixrm_${file}_1.cor.fq unfixrm_${file}_2.cor.fq
unzip unfixrm_${file}_1.cor_fastqc.zip
unzip unfixrm_${file}_2.cor_fastqc.zip
python /gpfs/data/rbeinart/Software/TranscriptomeAssemblyTools/RemoveFastqcOverrepSequenceReads.py -1 unfixrm_${file}_1.cor.fq -2 unfixrm_${file}_2.cor.fq -fql unfixrm_${file}_1.cor_fastqc/fastqc_data.txt -fqr unfixrm_${file}_2.cor_fastqc/fastqc_data.txt

done

# Assembly
/gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/Trinity --seqType fq --max_memory 150G --left rmoverrep_unfixrm_R1928-5_241_host_1.cor.fq,rmoverrep_unfixrm_R1931-6_461_host_1.cor.fq,rmoverrep_unfixrm_710_host_1.cor.fq --right rmoverrep_unfixrm_R1928-5_241_host_2.cor.fq,rmoverrep_unfixrm_R1931-6_461_host_2.cor.fq,rmoverrep_unfixrm_R1935-3_710_host_2.cor.fq --SS_lib_type RF --output B_septemdierum.PasaFly.Trinity --bfly_algorithm PASAFLY --CPU 24 --full_cleanup

# Transcript clustering to remove redundancies
cd-hit-est -i B_septemdierum.PasaFly.Trinity.fasta -o B_septemdierum.Trinity.merged95.fasta -c 0.95 -n 10 -r 1 -T 24 -d 0 -M 50000 -g 1 -aS 0.95

# This command is just for quality control to see how well the assembly reflects the original input reads
bowtie2 -p 24 -q --no-unal -x B_septemdierum.Trinity.merged95.fasta -1 rmoverrep_unfixrm_R1928-5_241_host_1.cor.fq,rmoverrep_unfixrm_R1931-6_461_host_1.cor.fq,rmoverrep_unfixrm_710_host_1.cor.fq -2 rmoverrep_unfixrm_R1928-5_241_host_2.cor.fq,rmoverrep_unfixrm_R1931-6_461_host_2.cor.fq,rmoverrep_unfixrm_R1935-3_710_host_2.cor.fq 2>align_stats.txt

# ORF prediction
TransDecoder.LongOrfs -t B_septemdierum.Trinity.merged95.fasta

blastp -query B_septemdierum.Trinity.merged95.fasta.transdecoder_dir/longest_orfs.pep -db /gpfs/data/rbeinart/Databases/uniref90.fasta -max_target_seqs 1 -outfmt "6 std stitle" -evalue 1e-5 -num_threads 24 > blastp.outfmt6
/gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/util/misc/blast_outfmt6_group_segments.pl blastp.outfmt6 B_septemdierum.Trinity.merged95.fasta.transdecoder_dir/longest_orfs.pep /gpfs/data/rbeinart/Databases/uniref90.fasta > blastp.outfmt6.grouped
hmmscan --cpu 24 --domtblout pfam.domtblout /gpfs/data/rbeinart/Databases/Pfam-A.hmm B_septemdierum.Trinity.merged95.fasta.transdecoder_dir/longest_orfs.pep

TransDecoder.Predict -t B_septemdierum.Trinity.merged95.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

# Removal of remaining bacterial contaminants with BlobTools
cat rmoverrep_unfixrm*_host.1.cor.fq > B_septemdierum_1.fq
cat rmoverrep_unfixrm*_host.2.cor.fq > B_septemdierum_2.fq

bowtie2-build B_septemdierum.Trinity.merged95.fasta.transdecoder.cds B_septemdierum.Trinity.merged95.fasta.transdecoder.cds
bowtie2 --very-sensitive -p 24 -x B_septemdierum.Trinity.merged95.fasta.transdecoder.cds -1 B_septemdierum_1.fq -2 B_septemdierum_2.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > B_septemdierum.sorted.bam
samtools view -h -@ 24 -o B_septemdierum.sam B_septemdierum.sorted.bam 

blobtools create -i B_septemdierum.Trinity.merged95.fasta.transdecoder.cds -o B_septemdierum -s B_septemdierum.sam -t B_septemdierum_hits.txt
blobtools view -i B_septemdierum.blobDB.json -o B_septemdierum -r all
blobtools plot -i B_septemdierum.blobDB.json --notitle -r superkingdom --format pdf --colours colors.txt
grep "Eukaryota" B_septemdierum.blobDB.table.txt > seqs_to_keep.txt
grep "no-hit" B_septemdierum.blobDB.table.txt >> seqs_to_keep.txt
perl -anle 'print $F[0]' seqs_to_keep.txt > euk_seqs.txt
seqtk subseq B_septemdierum.Trinity.merged95.fasta.transdecoder.cds euk_seqs.txt > B_septemdierum.Trinity.merged95.filtered.cds

# Assembly statistics
perl /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/util/TrinityStats.pl B_septemdierum.Trinity.merged95.filtered.cds > TrinityStats.txt

# Assembly completeness
source activate busco

busco -i B_septemdierum.Trinity.merged95.filtered.cds -l mollusca_odb10 -o busco_mollusca_merged -m tran -c 24

