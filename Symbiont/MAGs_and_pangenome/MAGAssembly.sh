#!/bin/bash
#SBATCH -J MAGReconstruction
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=250g
#SBATCH -o MAGReconstruction_J2432_B%a.out
#SBATCH -e MAGReconstruction_J2432_B%a.err
#SBATCH --array=238-257

module load samtools/1.12
module load bowtie2/2.3.5.1
module load hmmer
module load idba
module load prodigal
module load R
module load ruby

# Examples are shown based on Tow Cam population J2432
metaspades.py -1 J2432_B${SLURM_ARRAY_TASK_ID}_R1_clean.fastq -2 J2432_B${SLURM_ARRAY_TASK_ID}_R2_clean.fastq -k 21,31,41,51,61,71,81,91,101,111,121 -t 24 -o metaSpades_J2432_B${SLURM_ARRAY_TASK_ID}
cd metaSpades_J2432_B${SLURM_ARRAY_TASK_ID}
bowtie2-build scaffolds.fasta scaffolds.fasta
bowtie2 -p 24 -x scaffolds.fasta -1 ../J2432_B${SLURM_ARRAY_TASK_ID}_R1_clean.fastq -2 ../J2432_B${SLURM_ARRAY_TASK_ID}_R2_clean.fastq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > J2432_B${SLURM_ARRAY_TASK_ID}.bowtie2.sorted.bam
source activate metabat2
jgi_summarize_bam_contig_depths --outputDepth depth.txt J2432_B${SLURM_ARRAY_TASK_ID}.bowtie2.sorted.bam
metabat2 -i scaffolds.fasta -a depth.txt -m 1500 -o metabat/J2432_B${SLURM_ARRAY_TASK_ID} --unbinned -t 24
mkdir maxbin
run_MaxBin.pl -contig scaffolds.fasta -reads ../J2432_B${SLURM_ARRAY_TASK_ID}_R1_clean.fastq -reads2 ../J2432_B${SLURM_ARRAY_TASK_ID}_R2_clean.fastq -out maxbin/J2432_B${SLURM_ARRAY_TASK_ID} -thread 24 -min_contig_length 500
mkdir maxbin/bins
mv maxbin/*fasta maxbin/bins/.
mv maxbin/J2432_B${SLURM_ARRAY_TASK_ID}.noclass maxbin/bins/J2432_B${SLURM_ARRAY_TASK_ID}.noclass.fasta
mv maxbin/J2432_B${SLURM_ARRAY_TASK_ID}.tooshort maxbin/bins/J2432_B${SLURM_ARRAY_TASK_ID}.tooshort.fasta
source deactivate
source activate metawrap-env
metawrap bin_refinement -o metawrap -A metabat -B maxbin/bins -t 24 -m 250 -c 70 -x 10
source deactivate
/gpfs/data/rbeinart/Software/DAS_Tool/src/Fasta_to_Scaffolds2Bin.sh -i maxbin/bins -e fasta > maxbin.scaffolds2bin.tsv
/gpfs/data/rbeinart/Software/DAS_Tool/src/Fasta_to_Scaffolds2Bin.sh -i metabat -e fa > metabat.scaffolds2bin.tsv
ln -s metawrap/metawrap_70_10_bins.contigs metawrap.scaffolds2bin.tsv
DAS_Tool -i maxbin.scaffolds2bin.tsv,metabat.scaffolds2bin.tsv,metawrap.scaffolds2bin.tsv -l MaxBin2,metaBAT2,metaWRAP --score_threshold 0 -c scaffolds.fasta -o DAS_Tool/J2432_B${SLURM_ARRAY_TASK_ID} --write_bins 1 -t 24

