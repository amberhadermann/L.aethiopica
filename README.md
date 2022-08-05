# Genome diversity of _Leishmania aethiopica_ from Ethiopia

This github repository contains all scripts written in the light of this study.

## Contents
1. [VARIANTCALLING.sh][VC]: A script designed to run on the Vlaams Supercomputer Centrum (VSC(r); sbatch format) to start the analysis by preparing the read files and mapping them against the reference genome. 
2. [SNP_INDEL_FILTER_ALL.sh][ALL]: A script designed to run on the Vlaams Supercomputer Centrum (VSC(r); sbatch format) to combine VCF files of all isolates (ALL41) and filter Single Nucleotide Polymorphisms (SNPs) and INsertions and DELetions (INDELs). In addition to creating all intermediate VCF files needed for further analyses on the ALL41 dataset.
3. [SNP_INDEL_FILTER_AETH.sh][AETH]: A script designed to run on the Vlaams Supercomputer Centrum (VSC(r); sbatch format) to combine VCF files of the _L. aethiopica_ isolates (AETH20) and filter Single Nucleotide Polymorphisms (SNPs) and INsertions and DELetions (INDELs). 
4. [POPULATION_STRUCTURE.sh](https://github.com/amberhadermann/L.aethiopica/blob/main/POPULATION_STRUCTURE.sh):  A script designed to run locally as to perform LDdecay and ADMIXTURE analysis. 
5. [Thesis_Rscript.R](https://github.com/amberhadermann/L.aethiopica/blob/main/Thesis_Rscript.R): An R script designed to be run through Rstudio to do all statistical analyses and vizualisation used during this project.

## Run instructions per script
##### [VARIANTCALLING.sh][VC]
##### Input needed
Make sure to include the following items in your workingdirectory:

- Reference genome: Here, we used [TriTrypDB-54_LaethiopicaL147_Genome.fasta][TriTrypDB-54_LaethiopicaL147_Genome.fasta]. 
- Sequence read files of all isolates/samples in the format: *_R1.fastq.gz, *_R2.fastq.gz or *_R1.fastq.gz.
- [depth2window.py][dw] script, created by [Frederick Van Den Broeck](https://github.com/FreBio). 
##### Running on the commandline
Run this script in loop on the command line as for example: ```for i in sample1 sample2; do sbatch -J ${i}_variantcalling --export=${inputfile}=$i VARIANTCALLING.sh; done```
The output of this example will include a sorted and marked bam-file, a VCF-file, depth files and read-quality files of sample1 and sample2. Additionally, you will find error and input files, marked by _variantcalling.e and _variantcalling.o, per isolate listed.
 => **NOTE**: replace sample1 and sample2 with a list of all isolates you want processed.

##### [SNP_INDEL_FILTER_ALL.sh][ALL]
##### Input needed
Make sure to include the following items in your workingdirectory:

- Reference genome: Here, we used [TriTrypDB-54_LaethiopicaL147_Genome.fasta][TriTrypDB-54_LaethiopicaL147_Genome.fasta]. 
- All isolates/samples zipped VCF files of the ALL41 dataset.
 
##### Running on the commandline
Run this script on the command line as for example: ```sbatch -J SNP_INDEL_FILTER_ALL SNP_INDEL_FILTER_ALL.sh```
The output of this example will include VCF files of the entire ALL41 dataset, the dataset without hybrids and one per interesting chromosome found during the alternate allele frequencies analysis. Furthermore, producing files per statistic used while determining filter conditions. Additionally, you will find error and input files, marked by SNP_INDEL_FILTER_ALL.e and SNP_INDEL_FILTER_ALL.o.
 
##### [SNP_INDEL_FILTER_AETH.sh][AETH]
##### Input needed
Make sure to include the following items in your workingdirectory:

- Reference genome: Here, we used [TriTrypDB-54_LaethiopicaL147_Genome.fasta][TriTrypDB-54_LaethiopicaL147_Genome.fasta]. 
- All isolates/samples zipped VCF files of the AETH20 dataset.
 
##### Running on the commandline
Run this script on the command line as for example: ```sbatch -J SNP_INDEL_FILTER_AETH SNP_INDEL_FILTER_AETH.sh```
The output of this example will include VCF files of the entire AETH20 dataset and files per statistic used while determining filter conditions. Additionally, you will find error and input files, marked by SNP_INDEL_FILTER_AETH.e and SNP_INDEL_FILTER_AETH.o.

##### [POPULATION_STRUCTURE.sh](https://github.com/amberhadermann/L.aethiopica/blob/main/POPULATION_STRUCTURE.sh)
##### Input needed
Make sure to include the following items in your workingdirectory:

- Final AETH20 VCF file (esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP5.GQ40.vcf.gz): output of [SNP_INDEL_FILTER_AETH.sh][AETH].
- Locally installed [ADMIXTURE](https://dalexander.github.io/admixture/download.html) and [PopLDdecay](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjz8IrWsq_5AhWFIMUKHVQ1AY4QFnoECBEQAQ&url=https%3A%2F%2Fgithub.com%2FBGI-shenzhen%2FPopLDdecay&usg=AOvVaw2kKoF_9WuHoGlW930inm50).

##### Running on the commandline
Run this script on the command line by running the file: ```./POPULATION_STRUCTURE.sh```.

##### [Thesis_Rscript.R](https://github.com/amberhadermann/L.aethiopica/blob/main/Thesis_Rscript.R)
##### Input needed
Make sure to set your working directory containing the output from all previously mentioned scripts in Rstudio.
However, for some analyses you will need 012 files, that are created on the commandline using ```vcftools --gzvcf ${input}.vcf.gz --012 --out ${input}```. The ${input} must be either looped or replaced with the names given in the R code. As to run vcftools, you need a locally installed [BCFtools](http://www.htslib.org/download/).


[TriTrypDB-54_LaethiopicaL147_Genome.fasta]: <https://tritrypdb.org/common/downloads/release-54/LaethiopicaL147>
[VC]: <https://github.com/amberhadermann/L.aethiopica/blob/main/VARIANTCALLING.sh>
[dw]: <https://github.com/FreBio/mytools/blob/master/depth2window.py>
[ALL]: <https://github.com/amberhadermann/L.aethiopica/blob/main/SNP_INDEL_FILTER_ALL.sh>
[AETH]: <https://github.com/amberhadermann/L.aethiopica/blob/main/SNP_INDEL_FILTER_AETH.sh>

