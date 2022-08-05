#!/usr/bin/env bash

#SBATCH --ntasks=1 --cpus-per-task=1 --mem=30g
#SBATCH --time=48:00:00
#SBATCH -o %x.o
#SBATCH -e %x.e
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=amber.hadermann@hotmail.be

### RE-Setting the environment:
module --force purge
module load calcua/2020a

## MODULES:
#module load GATK
#module load BCFtools
module load Python
#module load SAMtools
module load fastp
#module load BWA
module load Java
module load BioTools

#REF=TriTrypDB-54_LaethiopicaL147_Genome.fasta

# Combining gVCF files --> all samples into one vcf file
gatk CombineGVCFs -R ${REF} -V 103-83.vcf.gz -V 1123-81.vcf.gz -V 117-82.vcf.gz -V 130-83.vcf.gz -V 1464-85.vcf.gz -V 1561-87.vcf.gz -V 169-83.vcf.gz -V 32-83.vcf.gz -V 678-82.vcf.gz -V 68-83.vcf.gz -V 85-83.vcf.gz -V ERR1913337.vcf.gz -V ERR205814.vcf.gz -V ERR205815.vcf.gz -V ERR205816.vcf.gz -V ERR205817.vcf.gz -V ERR207776.vcf.gz -V ERR304754.vcf.gz -V ERR304759.vcf.gz -V ERR304762.vcf.gz -V ERR304763.vcf.gz -V ERR304766.vcf.gz -V GEREcl7.vcf.gz -V HUSSEN.vcf.gz -V L100.vcf.gz -V L100cl1.vcf.gz -V L127.vcf.gz -V L86.vcf.gz -V LEM2358cl3.vcf.gz -V LEM3464.vcf.gz -V LEM3469.vcf.gz -V LEM3469cl1.vcf.gz -V LEM3469cl5.vcf.gz -V LEM3469cl7.vcf.gz -V LEM3469cl8.vcf.gz -V LEM3469cl9.vcf.gz -V LEM3497.vcf.gz -V LEM3498.vcf.gz -V SRR1028158.vcf.gz -V ERR311388.vcf.gz -V WANDERA.vcf.gz -O esembled_ALL.vcf.gz

# Genotyping vcf file:
gatk GenotypeGVCFs -R ${REF} -V esembled_ALL.vcf.gz -O esembled_ALL.GENO.vcf.gz

# Separating SNPs and INDELs from each other:
gatk SelectVariants -R ${REF} -V esembled_ALL.GENO.vcf.gz -select-type SNP -O esembled_ALL.GENO.SNP.vcf.gz
gatk SelectVariants -R ${REF} -V esembled_ALL.GENO.vcf.gz -select-type INDEL -O esembled_ALL.GENO.INDEL.vcf.gz


# Hard filtering
##--> GATK recommendations for SNPs and INDELs
##--> filter for SNPs in SNP clusters: (e.g. clustersize=3 clusterwindow=10)
##--> Filter additionally on QUAL --> 100, 200, ... (or to extra high quality 1500 )
##--> filter additionally on format DP (individual variant read depth) (<10)
##--> filter additionally on format GQ (ind. genotype quality) (,10, <25, <50)

# filter step wise and check different quality statistics
### filtering based on gatk recommendations:
#SNP
gatk VariantFiltration  -R ${REF} -V esembled_ALL.GENO.SNP.vcf.gz --cluster-window-size 10 --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'GATKrecommended' -O esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.vcf.gz

gzip -cd esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.vcf.gz | grep -v ‘GATKrecommended’ | grep -v ‘SnpCluster’ | grep -v ‘\./\.’ | bgzip > esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.vcf.gz

bcftools filter -e ‘FMT/DP<5’ esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.vcf.gz | bcftools filter -e ‘FMT/GQ<20’ - | bgzip > esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP5.GQ20.vcf.gz

#INDEL
gatk VariantFiltration -R ${REF} -V esembled_ALL.GENO.INDEL.vcf.gz --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filter-name 'GATKrecommended'  -O esembled_ALL.GENO.INDEL.GATKrecom.vcf.gz

#filter check
bcftools query -f '%MQ\n' esembled_ALL.GENO.SNP.GATKrecom.own.PASS.vcf.gz > esembled_ALL.GENO.SNP.GATKrecom.own.PASS.MQ.txt
bcftools query -f '%QUAL\n' esembled_ALL.GENO.SNP.GATKrecom.own.PASS.vcf.gz > esembled_ALL.GENO.SNP.GATKrecom.own.PASS.QUAL.txt
bcftools query -f '[%GQ\t]\n' esembled_ALL.GENO.SNP.GATKrecom.own.PASS.vcf.gz> esembled_ALL.GENO.SNP.GATKrecom.own.PASS.GQ.txt
bcftools query -f '[%DP\t]\n' esembled_ALL.GENO.SNP.GATKrecom.own.PASS.vcf.gz > esembled_ALL.GENO.SNP.GATKrecom.own.PASS.FMTDP.txt
bcftools query -f '%DP\n' esembled_ALL.GENO.SNP.GATKrecom.own.PASS.vcf.gz> esembled_ALL.GENO.SNP.GATKrecom.own.PASS.DP.txt

#indexing VCF files
tabix -p vcf Esembled/esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP5.GQ20.vcf.gz

#No hybrids
bcftools view -a -s ^L86,LEM3469,LEM3469cl1,LEM3469cl5,LEM3469cl7,LEM3469cl8,LEM3469cl9 esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP5.GQ20.vcf.gz | bcftools view -e ‘ALT==“.”’ - | bgzip > esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP5.GQ20.NOHYBRIDS.vcf.gz

##Filtering per chromosome
for I in 06\:200000\-600000 09 10 11 15 18\:520000\-750000 20 22\:0\-300000 24;
do
bcftools view -r LaeL147_${I} Esembled/esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP5.GQ20.vcf.gz | bgzip > esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP5.GQ20.CHR${I}.vcf.gz
done
