#!/bin/bash
#SBATCH -J ClonalFrameML
#SBATCH -t 200:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o ClonalFrameML.out
#SBATCH -e ClonalFrameML.err

module unload java
module load raxml
module load mafft
module load trimal
module load perl
module load R

progressiveMauve --output=Thiodubiliella.xmfa --output-guide-tree=Thiodubiliella.tree --backbone-output=Thiodubiliella.backbone J2426_B105.fasta R1931-6_463.fasta J2432_B243.fasta J2432_B256.fasta R1935-3_712.fasta R1935-3_691.fasta R1928-5_239.fasta J2428_B189.fasta
stripSubsetLCBs Thiodubiliella.xmfa Thiodubiliella.xmfa.bbcols Thiodubiliella_core_genes.xmfa 500 8
perl /gpfs/data/rbeinart/cbreusing/Scripts/xmfa2fasta.pl --align --file Thiodubiliella_core_genes.xmfa > Thiodubiliella_concatenated_core_genes.fa
perl /gpfs/data/rbeinart/cbreusing/Scripts/RenameContigsv2.pl Thiodubiliella_concatenated_core_genes.fa Thiodubiliella.txt Thiodubiliella_concatenated_core_genes-renamed.fa
mafft --auto --thread 24 Thiodubiliella_concatenated_core_genes-renamed.fa > Thiodubiliella_realigend_core_genes.fa
trimal -in Thiodubiliella_realigend_core_genes.fa -out Thiodubiliella_cleaned_core_genes.fa -htmlout Thiodubiliella_trimal.html -resoverlap 0.75 -seqoverlap 80
python /gpfs/data/rbeinart/cbreusing/Scripts/fasta2phy.py Thiodubiliella_cleaned_core_genes.fa Thiodubiliella_cleaned_core_genes.phy
raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 - 100 -s Thiodubiliella_cleaned_core_genes.phy -n Thiodubiliella_core_genes.tre
ClonalFrameML RAxML_bestTree.Thiodubiliella_core_genes.tre Thiodubiliella_cleaned_core_genes.fa Thiodubiliella_clonalframe_standard.out -kappa 4.496380932 -emsim 100
ClonalFrameML RAxML_bestTree.Thiodubiliella_core_genes.tre Thiodubiliella_cleaned_core_genes.fa Thiodubiliella_clonalframe_per_branch.out -kappa 4.496380932 -embranch true -embranch_dispersion 0.1 -initial_values "0.436884 0.0157872 0.0859513"
run_gubbins.py Thiodubiliella_cleaned_core_genes.fa --bootstrap 100 --starting-tree RAxML_bestTree.Thiodubiliella_core_genes.tre --prefix Thiodubiliella --verbose --threads 24
Rscript /gpfs/data/rbeinart/cbreusing/Scripts/cfml_results.R Thiodubiliella_clonalframe_standard.out
Rscript /gpfs/data/rbeinart/cbreusing/Scripts/cfml_results.R Thiodubiliella_clonalframe_per_branch.out

