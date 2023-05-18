#!/bin/bash
#SBATCH -J mergebam
#SBATCH -t 100:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o mergebam.out
#SBATCH -e mergebam.err

module load bamtools/2.4.1 
module load samtools
module load gatk/4.1.6.0
module load picard-tools/2.17.11

num=798471

for file in `cat ABE_norm.list`
do
count=`samtools view -c ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam`
if [ $num -le $count ]
    then
    frac=`bc -l <<< $num/$count`
    samtools view -h -bs $frac ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam > ${file}.Thiodubiliella.recomb.ABE.norm.bam
fi
ls ${file}.Thiodubiliella.recomb.ABE.norm.bam >> bam.ABE.norm.list
done

bamtools merge -list bam.ABE.norm.list -out Thiodubiliella.merged.ABE.norm.bam

for file in `cat THM_norm.list`
do
count=`samtools view -c ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam`
if [ $num -le $count ]
    then
    frac=`bc -l <<< $num/$count`
    samtools view -h -bs $frac ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam > ${file}.Thiodubiliella.recomb.THM.norm.bam
fi
ls ${file}.Thiodubiliella.recomb.THM.norm.bam >> bam.THM.norm.list
done

bamtools merge -list bam.THM.norm.list -out Thiodubiliella.merged.THM.norm.bam

for file in `cat TM_norm.list`
do
count=`samtools view -c ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam`
if [ $num -le $count ]
    then
    frac=`bc -l <<< $num/$count`
    samtools view -h -bs $frac ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam > ${file}.Thiodubiliella.recomb.TM.norm.bam
fi
ls ${file}.Thiodubiliella.recomb.TM.norm.bam >> bam.TM.norm.list
done

bamtools merge -list bam.TM.norm.list -out Thiodubiliella.merged.TM.norm.bam

for file in `cat TC_norm.list`
do
count=`samtools view -c ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam`
if [ $num -le $count ]
    then
    frac=`bc -l <<< $num/$count`
    samtools view -h -bs $frac ${file}.Thiodubiliella.notNormalized.softclipped.sorted.bam > ${file}.Thiodubiliella.recomb.TC.norm.bam
fi
ls ${file}.Thiodubiliella.recomb.TC.norm.bam >> bam.TC.norm.list
done

bamtools merge -list bam.TC.norm.list -out Thiodubiliella.merged.TC.norm.bam

ls J2423_B41.Thiodubiliella.recomb.ABE.norm.bam J2426_B104.Thiodubiliella.recomb.ABE.norm.bam BS_443.Thiodubiliella.recomb.ABE.norm.bam J2432_B238.Thiodubiliella.recomb.TC.norm.bam J2432_B251.Thiodubiliella.recomb.TC.norm.bam J2432_B244.Thiodubiliella.recomb.TC.norm.bam BS_714.Thiodubiliella.recomb.THM.norm.bam BS_699.Thiodubiliella.recomb.THM.norm.bam BS_716.Thiodubiliella.recomb.THM.norm.bam J2428_B195.Thiodubiliella.recomb.TM.norm.bam J2428_B200.Thiodubiliella.recomb.TM.norm.bam BS_249.Thiodubiliella.recomb.TM.norm.bam > bam.recomb.norm.list
bamtools merge -list bam.recomb.norm.list -out Thiodubiliella.merged.norm.bam


