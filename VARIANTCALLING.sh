#!/usr/bin/env bash

#SBATCH --ntasks=1 --cpus-per-task=1 --mem=30g
#SBATCH --time=24:00:00
#SBATCH -o %x.o
#SBATCH -e %x.e

### RE-Setting the environment:
module --force purge
module load calcua/2020a

## MODULES:
module load GATK
module load BCFtools
module load Python
module load fastp
module load BWA
module load Java
module load BioTools

REF=TriTrypDB-54_LaethiopicaL147_Genome.fasta
#inputfile= list of isolate names (in loop when submitting script)

## Indexing the reference genome
bwa index ${REF}

## Mapping paired end reads against reference genome
if find ${inputfile}_R2.fastq.gz;
then 
	bwa mem ${REF} ${inputfile}_R1.fastq.gz ${inputfile}_R2.fastq.gz > ${inputfile}.sam
	gzip -cd ${inputfile}_R1.fastq.gz | wc -l > ${inputfile}_R1.fastq.lines.txt
	gzip -cd ${inputfile}_R2.fastq.gz | wc -l > ${inputfile}_R2.fastq.lines.txt;
else bwa mem ${REF} ${inputfile}_R1.fastq.gz ${inputfile}_R3.fastq.gz > ${inputfile}.sam
	gzip -cd ${inputfile}_R1.fastq.gz | wc -l > ${inputfile}_R1.fastq.lines.txt
	gzip -cd ${inputfile}_R3.fastq.gz | wc -l > ${inputfile}_R3.fastq.lines.txt;
fi

## Conversion to binary file .bam
samtools view -b -o ${inputfile}.bam ${inputfile}.sam #sam is seq allignment file; bam file binary allignment; first output then input
samtools sort -o ${inputfile}.sort.bam ${inputfile}.bam 		
samtools index ${inputfile}.sort.bam 

rm ${inputfile}.bam | rm ${inputfile}.sam

gatk AddOrReplaceReadGroups -I ${inputfile}.sort.bam -O ${inputfile}.sort.RG.bam -SORT_ORDER coordinate -RGLB bar -RGPL illumina -RGSM ${inputfile} -CREATE_INDEX True -RGPU unit1 -RGID foo

## Marking duplicate reads
gatk MarkDuplicates -I ${inputfile}.sort.RG.bam -O ${inputfile}.sort.markdups.RG.bam -M ${inputfile}.sort.markdups.RG_metrics.txt
samtools index ${inputfile}.sort.markdups.RG.bam
samtools flagstat ${inputfile}.sort.markdups.RG.bam > ${inputfile}.sort.markdups.flagstat

rm ${inputfile}.sort.bam |rm ${inputfile}.sort.RG.bam

#need for fai reference file to work with gatk
samtools faidx ${REF} #unzip file if needed

# Variant calling with local re-de-novo assembly
gatk HaplotypeCaller -R ${REF} -I ${inputfile}.sort.markdups.RG.bam -O ${inputfile}.vcf.gz -ERC GVCF

#Depth files
samtools depth -a ${inputfile}.\sort\.markdups.RG.bam > ${inputfile}.depth
./depth2window_unzipped.py ${inputfile}.2000 ${inputfile}.depth 

#Read quality files 
samtools depth -Q 25 -q 25 -d 0 ${inputfile}.sort.markdups.RG.bam |awk '$3 > 30 {print $1"\t"$2"\t"$3}' | grep 'LaeL147_'| wc -l > ${inputfile}_25_25_30.txt
samtools depth -d 0 ${inputfile}.sort.markdups.RG.bam |awk '$3 > 30 {print $1"\t"$2"\t"$3}' | grep 'LaeL147_'| wc -l > ${inputfile}_30.txt
samtools depth -d 0 ${inputfile}.sort.markdups.RG.bam |awk '$3 > 20 {print $1"\t"$2"\t"$3}' | grep 'LaeL147_'| wc -l > ${inputfile}_20.txt
