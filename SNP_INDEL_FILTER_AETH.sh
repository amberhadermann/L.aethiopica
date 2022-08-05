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
module load SAMtools
module load fastp
module load BWA
module load Java

REF= TriTrypDB-54_LaethiopicaL147_Genome.fasta

# Combining gVCF files --> all samples into one vcf file
gatk CombineGVCFs -R ${REF} -V 1123-81.vcf.gz -V 117-82.vcf.gz -V 130-83.vcf.gz -V 1464-85.vcf.gz -V 1561-80.vcf.gz -V 169-83.vcf.gz -V 32-83.vcf.gz -V 678-82.vcf.gz -V 68-83.vcf.gz -V 85-83.vcf.gz -V GEREcl7.vcf.gz -V L100.vcf.gz -V L100cl1.vcf.gz -V L127.vcf.gz -V L86.vcf.gz -V LEM2357.vcf.gz -V LEM2358cl3.vcf.gz -V LEM3464.vcf.gz -V LEM3469.vcf.gz -V LEM3469cl1.vcf.gz -V LEM3469cl5.vcf.gz -V LEM3469cl7.vcf.gz -V LEM3469cl8.vcf.gz -V LEM3469cl9.vcf.gz -V LEM3497.vcf.gz -V LEM3498.vcf.gz -V WANDERA.vcf.gz -O esembled_Aeth.vcf.gz

# Genotyping vcf file:
gatk GenotypeGVCFs -R ${REF} -V esembled_Aeth.vcf.gz -O esembled_Aeth.GENO.vcf.gz

# Separating SNPs and INDELs from each other:
gatk SelectVariants -R ${REF} -V esembled_Aeth.GENO.vcf.gz -select-type SNP -O esembled_Aeth.GENO.SNP.vcf.gz
gatk SelectVariants -R ${REF} -V esembled_Aeth.GENO.vcf.gz -select-type INDEL -O esembled_Aeth.GENO.INDEL.vcf.gz


# Hard filtering
##--> GATK recommendations for SNPs and INDELs
##--> filter for SNPs in SNP clusters: (e.g. clustersize=3 clusterwindow=10)
##--> Filter additionally on QUAL --> 100, 200, ... (or to extra high quality 1500 )
##--> filter additionally on format DP (individual variant read depth) (<10)
##--> filter additionally on format GQ (ind. genotype quality) (,10, <25, <50)


# filter step wise and check different quality statistics
#SNP
gatk VariantFiltration  -R ${REF} -V esembled_Aeth.GENO.SNP.vcf.gz --cluster-window-size 10 --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'GATKrecommended' -O esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.vcf.gz

gzip -cd esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.vcf.gz | grep -v 'GATKrecommended' | grep -v 'SnpCluster' | grep -v '\./\.' | bgzip > esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.vcf.gz

bcftools filter -e 'FMT/DP<5' esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.vcf.gz | bcftools filter -e 'FMT/GQ<40' - | bgzip > esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP5.GQ40.vcf.gz

#INDEL
gatk VariantFiltration -R ${REF} -V esembled_Aeth.GENO.INDEL.vcf.gz --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filter-name 'GATKrecommended'  -O esembled_Aeth.GENO.INDEL.GATKrecom.vcf.gz

#Filter check
bcftools query -f '%MQ\n' esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.vcf.gz > esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ.txt
bcftools query -f '%QUAL\n' esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.vcf.gz > esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.QUAL.txt
bcftools query -f '[%GQ\t]\n' esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.vcf.gz> esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ.txt
bcftools query -f '[%DP\t]\n' esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.vcf.gz > esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP.txt
bcftools query -f '%DP\n' esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.vcf.gz> esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP.txt

