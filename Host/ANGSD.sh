#!/bin/bash
#SBATCH -J SNPCalling
#SBATCH -t 200:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o SNPCalling.out
#SBATCH -e SNPCalling.err

source activate lofreq
module load bowtie2/2.3.5.1
module load samtools
module load gatk/4.1.6.0
module load picard-tools/2.17.11
module load R
module load fastme

bowtie2-build B_septemdierum.Trinity.merged95.filtered.fasta B_septemdierum.Trinity.merged95.filtered.fasta
samtools faidx B_septemdierum.Trinity.merged95.filtered.fasta

for file in `cat filelist.txt`
do
bowtie2 --very-sensitive-local -p 24 -x B_septemdierum.Trinity.merged95.filtered.fasta -1 ../${file}_R1_clean.fastq -2 ../${file}_R2_clean.fastq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.sorted.bam
java -Xmx10g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar MarkDuplicates I=${file}.sorted.bam O=${file}.dedup.bam M=${file}.metrics.txt
samtools index ${file}.dedup.bam
lofreq viterbi -f B_septemdierum.Trinity.merged95.filtered.fasta -k ${file}.dedup.bam | samtools sort -@ 24 - > ${file}.realigned.bam
lofreq indelqual -f B_septemdierum.Trinity.merged95.filtered.fasta --dindel -o ${file}.indelqual.bam ${file}.realigned.bam
samtools index ${file}.indelqual.bam
done

angsd -P 24 -bam bam.list -ref B_septemdierum.Trinity.merged95.filtered.fasta -out ANGSD/B_septemdierum.qc -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 30 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 1000
Rscript /gpfs/data/rbeinart/Software/ngsTools/Scripts/plotQC.R ANGSD/B_septemdierum.qc

for POP in TC THM ABE TM;
do
ALL=`echo $(cat ${POP}.list | wc -l)`
FRAC=`echo $(cat ${POP}.list | wc -l)*0.75 | bc -l`
IND=`printf "%.*f\n" 0 ${FRAC}`
angsd -P 24 -bam ${POP}.list -ref B_septemdierum.Trinity.merged95.filtered.fasta -gl 1 -baq 1 -C 50 -minInd ${IND} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 260 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-6 -doGlf 3 -minMaf 0.01 -skipTriallelic 1 -out ANGSD/${POP}
NSITES=`zcat ANGSD/${POP}.mafs.gz | tail -n+2 | wc -l`
zcat ANGSD/${POP}.glf.gz > ANGSD/${POP}.glf
/gpfs/data/rbeinart/Software/ngsTools/ngsF/ngsF.sh --n_ind ${ALL} --n_sites ${NSITES} --glf ANGSD/${POP}.glf --out ANGSD/${POP}.indF
angsd -P 24 -bam ${POP}.list -ref B_septemdierum.Trinity.merged95.filtered.fasta -anc B_septemdierum.Trinity.merged95.filtered.fasta -gl 1 -baq 1 -C 50 -minInd ${IND} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 260 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 4 -doSaf 2 -indF ANGSD/${POP}.indF -out ANGSD/${POP}
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/${POP}.saf.idx -tole 1e-6 -maxIter 5000 -P 24 -fold 1 > ANGSD/$POP.sfs
angsd -P 24 -bam ${POP}.list -ref B_septemdierum.Trinity.merged95.filtered.fasta -anc B_septemdierum.Trinity.merged95.filtered.fasta -gl 1 -baq 1 -C 50 -minInd ${IND} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 260 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 4 -doSaf 2 -indF ANGSD/${POP}.indF -doThetas 1 -pest ANGSD/${POP}.sfs -out ANGSD/${POP}
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS saf2theta ANGSD/${POP}.saf.idx -outname ANGSD/${POP} -sfs ANGSD/${POP}.sfs -fold 1
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/thetaStat do_stat ANGSD/${POP}.thetas.idx
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/thetaStat do_stat ANGSD/${POP}.thetas.idx -win 500 -step 100 -outnames ANGSD/${POP}.thetas
done

cat ANGSD/TC.indF ANGSD/THM.indF ANGSD/ABE.indF ANGSD/TM.indF > ANGSD/B_septemdierum.indF 

# Analysis was performed with a minimum individual proportion of 66% and 75%
angsd -P 24 -bam bam.list -ref B_septemdierum.Trinity.merged95.filtered.fasta -gl 1 -baq 1 -C 50 -minInd 70 -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 260 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-6 -dosnpstat 1 -doHWE 1 -sb_pval 0.05 -hetbias_pval 0.05 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGeno 2 -minMaf 0.01 -indF ANGSD/B_septemdierum.indF -skipTriallelic 1 -out ANGSD/B_septemdierum_66

/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/TC.saf.idx ANGSD/THM.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/TC-THM.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/TC.saf.idx ANGSD/ABE.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/TC-ABE.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/TC.saf.idx ANGSD/TM.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/TC-TM.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/THM.saf.idx ANGSD/ABE.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/THM-ABE.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/THM.saf.idx ANGSD/TM.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/THM-TM.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/ABE.saf.idx ANGSD/TM.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/ABE-TM.folded.sfs

/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/TC.saf.idx ANGSD/THM.saf.idx -sfs ANGSD/TC-THM.folded.sfs -fold 1 -fstout ANGSD/TC-THM -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/TC-THM.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/TC-THM.fst.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats2 ANGSD/TC-THM.fst.idx -win 500 -step 100 -whichFST 1 -maxIter 5000 > ANGSD/TC-THM.fst.window.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/TC.saf.idx ANGSD/ABE.saf.idx -sfs ANGSD/TC-ABE.folded.sfs -fold 1 -fstout ANGSD/TC-ABE -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/TC-ABE.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/TC-ABE.fst.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats2 ANGSD/TC-ABE.fst.idx -win 500 -step 100 -whichFST 1 -maxIter 5000 > ANGSD/TC-ABE.fst.window.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/TC.saf.idx ANGSD/TM.saf.idx -sfs ANGSD/TC-TM.folded.sfs -fold 1 -fstout ANGSD/TC-TM -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/TC-TM.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/TC-TM.fst.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats2 ANGSD/TC-TM.fst.idx -win 500 -step 100 -whichFST 1 -maxIter 5000 > ANGSD/TC-TM.fst.window.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/THM.saf.idx ANGSD/ABE.saf.idx -sfs ANGSD/THM-ABE.folded.sfs -fold 1 -fstout ANGSD/THM-ABE -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/THM-ABE.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/THM-ABE.fst.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats2 ANGSD/THM-ABE.fst.idx -win 500 -step 100 -whichFST 1 -maxIter 5000 > ANGSD/THM-ABE.fst.window.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/THM.saf.idx ANGSD/TM.saf.idx -sfs ANGSD/THM-TM.folded.sfs -fold 1 -fstout ANGSD/THM-TM -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/THM-TM.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/THM-TM.fst.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats2 ANGSD/THM-TM.fst.idx -win 500 -step 100 -whichFST 1 -maxIter 5000 > ANGSD/THM-TM.fst.window.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/ABE.saf.idx ANGSD/TM.saf.idx -sfs ANGSD/ABE-TM.folded.sfs -fold 1 -fstout ANGSD/ABE-TM -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/ABE-TM.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/ABE-TM.fst.txt 
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats2 ANGSD/ABE-TM.fst.idx -win 500 -step 100 -whichFST 1 -maxIter 5000 > ANGSD/ABE-TM.fst.window.txt



