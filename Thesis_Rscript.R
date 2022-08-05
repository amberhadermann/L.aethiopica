#########################################
#####        QUALITY PLOTTING       #####
#####     Leishmania AETHIOPICA     #####
#########################################

##Glossery
#Format Depth (FMTDP)
#Genotype quality (GQ)
#Mapping quality (MQ)
#Total depth across all genotypes (DP)
#Overall SNP quality (QUAL)

#######################
### 1. ALL isolates ###
#######################

#  import packages:
library(ggplot2)

#Data import/prep:
esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.QUAL <- read.table("esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.QUAL.txt", quote="\"", comment.char="")
esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ <- read.table("esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ.txt", quote="\"", comment.char="")
  esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ$V1 <- as.numeric(esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ$V1) #creating a nummeric dataframe
esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ <- read.table("esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ.txt", quote="\"", comment.char="")
  for (i in c(1:41)){esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ[,i] <- as.numeric(esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ[,i])} #creating a nummeric dataframe
esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP <- read.table("esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP.txt", quote="\"", comment.char="")
  for (i in c(1:41)){esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP[,i] <- as.numeric(esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP[,i])} #creating a nummeric dataframe
esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP <- read.table("esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP.txt", quote="\"", comment.char="")

#Plotting:
pdf("QUAL_ALL.pdf")
print(hist(esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.QUAL$V1, main = "Overall SNP quality ALL41",xlab= "QAUL", breaks = 500000, xlim = c(0,50000)))
dev.off()

pdf("MQ_ALL.pdf")
print(hist(esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ$V1, main = "Mapping quality ALL41",xlab= "MQ", breaks = 500000, xlim = c(30,70), ylim = c(0,1000)))
dev.off()

pdf("GQ_ALL.pdf")
print(for (i in c(1:41)) {hist(esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ[,i],main = paste("Genotype quality ALL41", i),xlab= "GQ",breaks = 50, xlim = c(0,100), ylim= c(0,10000))})
dev.off()

pdf("FMTDP_ALL.pdf")
print(for (i in c(1:41)) {hist(esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP[,i],main = paste("Format Depth ALL41", i),xlab= "FMTDP",breaks = 50000, xlim = c(0,200))})
dev.off()

pdf("DP_ALL.pdf")
print(hist(esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP$V1, breaks = 50000, xlim = c(0,5000), main = "Total depth across all genotypes ALL41",xlab= "DP"))
dev.off()

########################
### 2. AETH isolates ###
########################

#Data import/prep:
esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.QUAL <- read.table("esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.QUAL.txt", quote="\"", comment.char="")
esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ <- read.table("esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ.txt", quote="\"", comment.char="")
  esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ$V1 <- as.numeric(esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ$V1) #creating a nummeric dataframe
esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ <- read.table("esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ.txt", quote="\"", comment.char="")
  for (i in c(1:20)){esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ[,i] <- as.numeric(esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ[,i])} #creating a nummeric dataframe
esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP <- read.table("esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP.txt", quote="\"", comment.char="")
  for (i in c(1:20)){esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP[,i] <- as.numeric(esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP[,i])}
esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP <- read.table("esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP.txt", quote="\"", comment.char="")

#Plotting:
QUAL_Aeth <- hist(esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.QUAL$V1,main = paste("QUAL AETH ", i), breaks = 50000, xlim = c(0,50000))
pdf("QUAL_Aeth.pdf")
print(hist(esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.QUAL$V1,main ="Overall SNP quality AETH20",xlab= "QAUL", breaks = 500000, xlim = c(0,50000)))
dev.off()

pdf("MQ_Aeth.pdf")
print(hist(esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.MQ$V1,main ="Mapping quality AETH20",xlab= "MQ", breaks = 500000, xlim = c(30,70), ylim = c(0,1000)))
dev.off()

pdf("GQ_Aeth.pdf")
print(for (i in c(1:20)) {hist(esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.GQ[,i],main = paste("Genotype quality AETH20 ", i),xlab= "GQ",breaks = 50, xlim = c(0,100), ylim= c(0,10000))})
dev.off()

pdf("FMTDP_Aeth.pdf")
print(for (i in c(1:20)) {hist(esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.FMTDP[,i],main = paste("Format Depth AETH20 ", i),xlab= "FMTDP",breaks = 50000, xlim = c(0,200))})
dev.off()

pdf("DP_Aeth.pdf")
print(hist(esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP$V1, breaks = 50000, xlim = c(0,5000), main = "Total depth across all genotypes AETH20",xlab= "DP"))
dev.off()

#########################################
#####        HET/HOM PLOTTING       #####
#####     Leishmania AETHIOPICA     #####
#########################################

#Data packages needed:
library(data.table)
library(ggplot2)
library(scales)
library(ggrepel)

###########################
### 1. data import/prep ###
###########################

#function for reading in 012 files into R:
read.geno <- function(file) {
  require(data.table)
  geno <- fread(file, data.table = F, header = F)[,-1]
  rownames(geno) <- as.character(read.table(paste(file, 'indv', sep='.'))[,1])
  genopos <- read.table(paste(file, 'pos', sep='.'))
  colnames(geno) <- as.character(paste(genopos[,1], genopos[,2], sep=';'))
  return(as.data.frame(geno))
}

data_all <- read.geno('esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP5.GQ20.012')
data_aeth <- read.geno('esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP5.GQ40.012')

###########################
### 2. Counting HOM/HET ###
###########################
hetsites_aeth <- apply(data_aeth, 1, function(x) sum(x==1)) #counting heterozygous sites
  hetsites_TABLE_aeth <- as.data.frame(hetsites_aeth)
homsites_aeth <- apply(data_aeth, 1, function(x) sum(x==2)) #counting homozygous sites
  homsites_TABLE_aeth <- as.data.frame(homsites_aeth)

hetsites_all <- apply(data_all, 1, function(x) sum(x==1)) #counting heterozygous sites
  hetsites_TABLE_all <- as.data.frame(hetsites_all)
homsites_all <- apply(data_all, 1, function(x) sum(x==2)) #counting homozygous sites
  homsites_TABLE_all <- as.data.frame(homsites_all)

Freq<-cbind.data.frame(homsites,hetsites) #creating 1 dataframe including both heterozygous and homozygous site counts
setDT(Freq, keep.rownames = "SampleName")
#categorizing all includes isolates according to species
Freq$Type <- with(Freq, ifelse(startsWith(Freq$SampleName,"L100"),"L. aethiopica",
                               ifelse(startsWith(Freq$SampleName,"ERR2"), "L. donovani complex",
                                      ifelse(startsWith(Freq$SampleName,"ERR30"),"L. donovani complex",
                                             ifelse(startsWith(Freq$SampleName,"ERR30"),"L. donovani complex",
                                                    ifelse(startsWith(Freq$SampleName,"HUS"),"L. donovani complex",
                                                           ifelse(startsWith(Freq$SampleName,"ERR1"),"L. donovani complex",
                                                                  ifelse(startsWith(Freq$SampleName,"ERR31"),"L. tropica",
                                                                         ifelse(startsWith(Freq$SampleName,"LEM3469"),"L. aethiopica/L. donovani hybrid",
                                                                                ifelse(startsWith(Freq$SampleName,"L86"),"L. aethiopica/L. tropica hybrid",
                                                                                       ifelse(startsWith(Freq$SampleName,"S"),"L. major","L. aethiopica")))))))))))

Freq$Type <- as.factor(Freq$Type)

####################
### 3. Plotting  ###
####################

#total plot with all the different strains
TotalPlot<-ggplot(data = Freq, aes(x = Freq$homsites, y = Freq$hetsites, group=Freq$Type,label=Freq$SampleName))+
  geom_point(aes(shape=Freq$Type),size=4)+
  geom_text_repel(size=3,force = 10,max.overlaps = 17)+
  scale_shape_manual(values=c(5, 1, 17, 6,19,15, 0))+
  labs(x="Number of homozygous SNPs",y="Number of heterozygous SNPs", 
       title="Frequency plot of L. aethiopica, presumed hybrids and other species")+
  scale_x_continuous(breaks=seq(0,1000000,100000),limits=c(0,500000))+
  scale_y_continuous(breaks=seq(0,1000000,10000),limits=c(0,430000),
                     labels = scientific)+
  theme_bw()+
  theme(legend.title=element_blank(),legend.position = c(0.79, 0.85),
        legend.text = element_text(size=10))


pdf("TotalPlot.pdf")
print(TotalPlot)
dev.off()

#########################################
#####        FIXED SNP COUNT        #####
#####     Leishmania AETHIOPICA     #####
#########################################

###########################
### 1. data import/prep ###
###########################

#catergorizing all samples according to species
Laeth <- c("103-83","1123-81","117-82","130-83","1464-85","1561-87","169-83","32-83","678-82","68-83","85-83","GEREcl7","L100","L100cl1","L127","LEM2358cl3","LEM3464","LEM3497","LEM3498","WANDERA")
Lhyb1 <- c("LEM3469")#c("LEM3469","LEM3469cl1","LEM3469cl5","LEM3469cl7","LEM3469cl8","LEM3469cl9")  #L. aethiopica/L. donovani hybrid (and clones)
Lhyb2 <- c("L86") #L. aethiopica/L. tropica hybrid
Lmajo <- "SRR1028158"
Ltrop <- "ERR311388"
Ldono <- c("ERR1913337","ERR205814","ERR205815","ERR205816","ERR205817","ERR207776","ERR304754","ERR304759","ERR304762","ERR304763","ERR304766","HUSSEN")

data <- read.geno('esembled_ALL.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP5.GQ20.012')
freqs<- matrix(ncol = ncol(data), nrow=6)
freqs[1,] <- colSums(data[Laeth,])/(2*length(Laeth))
freqs[2,] <- colSums(data[Lhyb1,])/(2*length(Lhyb1))
freqs[3,] <- colSums(data[Lhyb2,])/(2*length(Lhyb2))
freqs[4,] <- colSums(data[Lmajo,])/(2*length(Lmajo))
freqs[5,] <- colSums(data[Ltrop,])/(2*length(Ltrop))
freqs[6,] <- colSums(data[Ldono,])/(2*length(Ldono))
rownames(freqs) <- c('Laeth','Lhyb1','Lhyb2','Lmajo','Ltrop','Ldono')

#################################
### 2. Fixed SNP counting WGS ###
#################################
## Fixed SNPs per species
Ltrop.fixed <- which(freqs['Ltrop',] == 1 & freqs['Ldono',] == 0 & freqs['Lmajo',] == 0 & freqs['Laeth',] == 0)
#52,958
Ldono.fixed <- which(freqs['Ltrop',] == 0 & freqs['Ldono',] == 1 & freqs['Lmajo',] == 0 & freqs['Laeth',] == 0)
#294,890
Lmajo.fixed <- which(freqs['Ltrop',] == 0 & freqs['Ldono',] == 0 & freqs['Lmajo',] == 1 & freqs['Laeth',] == 0)
#253,458

## Heterozygous SNPs in L86
L86het <- which(freqs['Lhyb2',] == 0.5) ## 164,566 heterozygous SNPs
    ##SNP fraction only present in fixed species
length(which(Ltrop.fixed %in% L86het))/length(Ltrop.fixed) ## 0.8503531
length(which(Ldono.fixed %in% L86het))/length(Ldono.fixed) ## 0.002516192
length(which(Lmajo.fixed %in% L86het))/length(Lmajo.fixed) ## 0.002651327

## Heterozygous and homozygous SNPs in LEM3469
LEMhet <- which(freqs['Lhyb1',] == 0.5) ##  416,672 heterozygous SNPs
LEMhom <- which(freqs['Lhyb1',] == 1) ##  10,229 homozygous SNPs
    ##SNP fraction only present in fixed species
length(which(Ltrop.fixed %in% LEMhet))/length(Ltrop.fixed) #5.664866e-05
length(which(Ldono.fixed %in% LEMhet))/length(Ldono.fixed) #0.9096612
length(which(Ldono.fixed %in% LEMhom))/length(Ldono.fixed) #6.78219e-05
length(which(Lmajo.fixed %in% LEMhet))/length(Lmajo.fixed) #0.0001262537


#####################################
### 2. Fixed SNP counting per chr ###
#####################################

## SNPs per chromosome
library(stringr)
listchromosomes <- str_split(colnames(data), ';', simplify = T)[,1]
chromosomes <- unique(listchromosomes)
chromosomes <- chromosomes[1:36]

countresultsHUSSEN <- countresultsERR205817 <- countresultsL127 <- countresultsLEM3469 <- countresultsL86 <- matrix(ncol=2,nrow=length(chromosomes))

for (i in 1:length(chromosomes)) {
  data.sub <- data['LEM3469',which(listchromosomes == chromosomes[i])]
  countresultsLEM3469[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['L127',which(listchromosomes == chromosomes[i])]
  countresultsL127[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['L86',which(listchromosomes == chromosomes[i])]
  countresultsL86[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['ERR205817',which(listchromosomes == chromosomes[i])]
  countresultsERR205817[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['HUSSEN',which(listchromosomes == chromosomes[i])]
  countresultsHUSSEN[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
}


countresultsLEM3469cl9 <-countresultsLEM3469cl8 <-countresultsLEM3469cl7 <-countresultsLEM3469cl5 <- countresultsLEM3469cl1 <- countresultsLEM3469 <- countresultsL86 <- matrix(ncol=2,nrow=length(chromosomes))


#Hybrids:
for (i in 1:length(chromosomes)) {
  data.sub <- data['LEM3469',which(listchromosomes == chromosomes[i])]
  countresultsLEM3469[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['LEM3469cl1',which(listchromosomes == chromosomes[i])]
  countresultsLEM3469cl1[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['L86',which(listchromosomes == chromosomes[i])]
  countresultsL86[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['LEM3469cl5',which(listchromosomes == chromosomes[i])]
  countresultsLEM3469cl5[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['LEM3469cl7',which(listchromosomes == chromosomes[i])]
  countresultsLEM3469cl7[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['LEM3469cl8',which(listchromosomes == chromosomes[i])]
  countresultsLEM3469cl8[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
  data.sub <- data['LEM3469cl9',which(listchromosomes == chromosomes[i])]
  countresultsLEM3469cl9[i,] <- table(as.character(data.sub))[2:3]#/sum(table(as.character(data.sub))[2:3])
}

#Plotting:
pdf("HOMHET_CHR.pdf")
{
  par(mar=c(6,6,2,0.5))
  par(mfrow=c(1,4))
  barplot(t(countresultsLEM3469), names.arg = chromosomes, las = 2, legend.text = c('HET', 'HOM'), horiz = T, args.legend = list(x = "bottomright", inset = c(4, 0.10)))
  title(main = "L.aethiopica/L. donovani hybrid", font.main = 4)
  title(main = "B", adj  = 0 )
  title(sub = "LEM3469")
  barplot(t(countresultsERR205817), las = 2, legend.text = c('HET', 'HOM'), horiz = T, args.legend = list(x = "bottomright", inset = c(4, 0.10)), xlim=c(0,40000))
  title(main = "L. donovani", font.main = 4)
  title(sub = "ERR205817")
  barplot(t(countresultsL86), las = 2, legend.text = c('HET', 'HOM'), horiz = T, args.legend = list(x = "bottomright", inset = c(4, 0.10)))
  title(main = "L.aethiopica", font.main = 2)
  title(sub = "L127")
  barplot(t(countresultsL86), las = 2, legend.text = c('HET', 'HOM'), horiz = T, args.legend = list(x = "bottomright", inset = c(-0.5, 0.10), cex=1, bty = "n"), xlim=c(0,40000))
  title(main = "L.aethiopica/L. tropica hybrid", font.main = 4)
  title(sub = "L86")
}
dev.off() 



#########################################
#####        SNP DIFFERENCE         #####
#####     Leishmania AETHIOPICA     #####
#########################################

#Counting the amount of shared SNPs between two isolates
done <- vector()
results <- matrix(ncol = nrow(geno), nrow=nrow(geno))
for (i in 1:nrow(geno)) {
  done <- append(x = done, values = i)
  for (j in 1:nrow(geno)) {
    if (j %in% done) {
      results[i,j] <- NA
    } else {
      results[i,j] <- length(c(which(geno[i,]==2 & geno[j,]==0),which(geno[i,]==0 & geno[j,]==2)))
    }
  }
}

geno_2 <- read.geno('esembled_Aeth.18.NOL100cl1.NO678_82.012')
results_2 <- matrix(ncol = nrow(geno_2), nrow=nrow(geno_2))
for (i in 1:nrow(geno_2)) {
  done <- append(x = done, values = i)
  for (j in 1:nrow(geno_2)) {
    if (j %in% done) {
      results_2[i,j] <- NA
    } else {
      results_2[i,j] <- length(c(which(geno_2[i,]==2 & geno_2[j,]==0),which(geno_2[i,]==0 & geno_2[j,]==2)))
    }
  }
}

geno_3 <- read.geno('esembled_Aeth.pop3.012')
results_3 <- matrix(ncol = nrow(geno_3), nrow=nrow(geno_3))
for (i in 1:nrow(geno_3)) {
  done <- append(x = done, values = i)
  for (j in 1:nrow(geno_3)) {
    if (j %in% done) {
      results_3[i,j] <- NA
    } else {
      results_3[i,j] <- length(c(which(geno_3[i,]==2 & geno_3[j,]==0),which(geno_3[i,]==0 & geno_3[j,]==2)))
    }
  }
}

#install.packages("writexl")
library(writexl)
results<- as.data.frame(result)
col.names(results_3) <- pop_3
row.names(results) <- info.samples$`WORK NAME`
write_xlsx(results,"DIFF_SNPS_AETH.xlsx")

#########################################
#####           ADMIXTURE           #####
#####     Leishmania AETHIOPICA     #####
#########################################

###  import packages:
library(ggplot2); library(readxl); library(forcats); 
library(tidyr); library(ggpubr); library(dplyr); library(tidyverse)


######################
### 1. CV - ERROR ####
######################

#Data import/prep:
CV_NODUP <- read.table("CV_NODUP.txt", quote="\"", comment.char="")
CV_NODUP$K <- c(1:10)
CV_NODUP

#GGplot 
CV_NODUP_plot <- {ggplot(CV_NODUP, aes(K, V4))+ 
    geom_line(linetype = "dashed")+
    geom_point()+
    labs(x = "K", y = "CV") +
    theme_minimal()+
    scale_x_continuous(breaks=seq(1, 10, 1))}

#######################
### 2. POP VISUAL  ####
#######################
#Defining population used:
pop_2 <- c("103-83","1123-81","117-82","130-83","1464-85","1561-87","169-83","32-83","68-83","85-83","GEREcl7","L100","L127","LEM2358cl3","LEM3464","LEM3497","LEM3498","WANDERA")

#K=2: Data import/prep + GGplot
{K2 <- read.table('esembled_Aeth.18.NOL100cl1.NO678_82.2.Q', header = F, sep = ' ')
  rownames(K2) <- pop_2
  K2.ord <- K2[order(K2$V1, K2$V2),]
  K2$Sample <- rownames(K2)
  K2 <- as_tibble(K2)
  plot_data <- K2 %>% 
    mutate(id = Sample) %>% 
    #arrange(Sample,K5$Sample) %>%
    #arrange(match(id, K5$Sample)) %%
    gather('pop', 'prob', V1:V2) %>% 
    group_by(id) %>% 
    mutate(likely_assignment = pop[which.max(prob)],
           assingment_prob = max(prob)) %>% 
    arrange(match(Sample, plot_data_K5$Sample)) %>%#,likely_assignment, desc(assingment_prob)) %>% 
    ungroup() %>% 
    mutate(id = forcats::fct_inorder(factor(id)))
  K2_plot <- ggplot(plot_data, aes(id, prob, fill = pop))+
    geom_col(color = "gray", size = 0.1) +scale_fill_manual(values=my_cols)+
    theme_minimal() + labs(x = "Individuals", y = "Ancestry", title = "K = 2") +
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
      panel.grid = element_blank()
    ) 
}

#K=5: : Data import/prep + GGplot
{K5 <- read.table('esembled_Aeth.18.NOL100cl1.NO678_82.5.Q', header = F, sep = ' ')
  rownames(K5) <-  pop_2
  K5$Sample <- pop_2
  K5 <- K5[order(K5$V1, K5$V2,K5$V3,K5$V4,K5$V5 ),]
  K5 <- as_tibble(K5)
  plot_data_K5 <- K5 %>% 
    mutate(id = Sample) %>% 
    gather('pop', 'prob', V1:V5) %>% 
    group_by(id) %>% 
    mutate(likely_assignment = pop[which.max(prob)],
           assingment_prob = max(prob)) %>% 
    arrange(likely_assignment, desc(assingment_prob)) %>% 
    ungroup() %>% 
    mutate(id = forcats::fct_inorder(factor(id)))
  K5_plot <- ggplot(plot_data_K5, aes(id, prob, fill = pop))+
    geom_col(color = "gray", size = 0.1) +scale_fill_manual(values=my_cols)+
    theme_minimal() + labs(x = "Individuals", y = "Ancestry", title = "K = 5") +
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
      panel.grid = element_blank()
    ) 
}

#Final plot
pdf('GGarrage_K2&5.pdf')
ggarrange(CV_NODUP_plot, ggarrange(K2_plot, K5_plot, labels = c("B", "C"),ncol = 1, nrow = 2), labels = c("A"), ncol= 2, nrow=1)
dev.off()


#########################################
##### Principal Component Analysis  #####
#####     Leishmania AETHIOPICA     #####
#########################################

##  import packages:
library(data.table)
library(stringr)
library(readxl)
library(adegenet)
library(ggplot2)

##  Functions required to load our data:
#Function to read in genotype files:
read.geno <- function(file) {
  geno <- fread(file, data.table = F, header = F)[,-1]
  rownames(geno) <- as.character(read.table(paste(file, 'indv', sep='.'))[,1]) #indv names as rownames
  genopos <- read.table(paste(file, 'pos', sep='.'))
  colnames(geno) <- as.character(paste(genopos[,1], genopos[,2], sep=';')) #Chr.pos as colnames
  return(as.data.frame(geno))
}
#Function to convert genotype file to genlight format:
geno2gl <- function(geno) {
  list <- as.list(as.data.frame(t(geno)))
  loci <- colnames(geno)
  positions <- as.character(lapply(str_split(loci,';'), function(x) x[2]))
  chromosomes <- as.character(lapply(str_split(loci,';'), function(x) x[1]))
  gl <- new('genlight', as.list(as.data.frame(t(geno))))
  gl@chromosome <- as.factor(chromosomes)
  gl@position <- as.factor(positions)
  gl@loc.names <- loci
  return(gl)
}
#Function to get Explained variance of PC axes:
get_explained_var_glPCA <- function(glPCA_obj){
  eig <- glPCA_obj$eig
  rbind(
    SD = sqrt(eig),
    Prop_of_Variance = eig/sum(eig),
    Cumulativa_Prop_of_Variance = cumsum(eig)/sum(eig))
}

######################
### 1. data import ###
######################
#Importing .012 file
geno <- read.geno(file = 'esembled_Aeth.GENO.SNP.GATKrecom.SNPCLUSTER.PASS.NOMISSING.DP5.GQ40.012')

#Importing META-DATA:
info.samples <- read_excel('SAMPELS_DATASET_2005_01.xlsx', sheet = 2)

################################
### 2. data conversion/ prep ###
################################
#Conversion to Genlight objects:
geno.gl <- geno2gl(geno)
#geno.gl_2 <- geno2gl(geno_2)

#META-DATA prep:
#info.samples$YEAR <- as.factor(info.samples$YEAR)
#info.samples$CITY <- as.factor(info.samples$CITY )
info.samples$POP <- as.factor(info.samples$POP_K5)
info.samples$POP2 <- as.factor(info.samples$POP_K2)

#Setting a grouping variable in the genlight format:
#geno.gl.YEAR <- geno.gl; geno.gl.YEAR$pop <- info.samples$YEAR
#geno.gl.CITY <- geno.gl; geno.gl.CITY$pop <- info.samples$CITY
geno.gl.POP <- geno.gl; geno.gl.POP$pop <- info.samples$POP
geno.gl.POP2 <- geno.gl; geno.gl.POP2$pop <- info.samples$POP2

#######################################
### 3. Principal Component Analysis ###
#######################################
pca.geno.gl <- glPca(geno.gl); pca.geno.gl
10
#pca.geno.gl.YEAR <- glPca(geno.gl.YEAR); pca.geno.gl.YEAR
#10
#pca.geno.gl.CITY <- glPca(geno.gl.CITY); pca.geno.gl.CITY
#10
pca.geno.gl.POP <- glPca(geno.gl.POP); pca.geno.gl.POP
10
pca.geno.gl.POP2 <- glPca(geno.gl.POP2); pca.geno.gl.POP2
10

#Get explained variance of PC axes:
get_explained_var_glPCA(pca.geno.gl)

#######################
### 4. Plotting PCA ###
#######################

#Quick n Dirty scatter plots:
par(mfrow= c(1,2))
scatter(pca.geno.gl)
scatter(pca.geno.gl,xax = 1, yax = 3)
scatter(pca.geno.gl,xax = 2, yax = 3)

#Some scatter plots with more colour in life:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D","#666666",'cornflowerblue', 'black')
par(mfrow=c(3,1)) #3rowsx1column scatterplots
par(mfrow=c(1,1))

#Defining color scheme
unique(geno.gl.POP$pop)
my_cols_2 <- c( "Black","#1B9E77","#D95F02", "#7570B3","#E7298A","chocolate4","#66A61E")
pal <- function(col=c("goldenrod","mediumorchid","coral"), border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
pal()
pal(my_cols_3)

my_cols_3 <- c( "Black","#1B9E77", "LIGHTBLUE","#D95F02")

plot(pca.geno.gl$scores[,1], pca.geno.gl$scores[,2],
     cex=2, pch=20, col=alpha(my_cols_3[geno.gl.POP2$pop],0.6), main='K=2 PC1/PC2', xlab= "", ylab="", xlim= c(-250,150))


#Final plots: different layouts
pdf('POP.PCA.SCATTER.pdf')
{layout(matrix(c(1,2,3,4,4,4), ncol=3, byrow=TRUE), heights=c(5, 1))
  par(mai=rep(0.5, 4))
  ##  PC1-2 
  plot(pca.geno.gl$scores[,1], pca.geno.gl$scores[,2],
       cex=3, pch=20, col=alpha(my_cols_2[geno.gl.POP$pop],0.6), main='K=5 PC1/PC2',
       xlab= "", ylab="", xlim= c(-250,100))
  text(pca.geno.gl$scores[,1], pca.geno.gl$scores[,2] , labels = info.samples$`WORK NAME`, cex = 1, pos = 2, col = "black")
  title(main = "A", adj  = 0)
  title(xlab= "PC 1 (31%)", ylab="PC 2 (18%)", line=2, cex.lab=1.2)
  ## PC1-3
  plot(pca.geno.gl$scores[,1], pca.geno.gl$scores[,3],
       cex=3, pch=20, col=alpha(my_cols_2[geno.gl.POP$pop],0.6), main='K=5 PC1/PC3', xlab= "", ylab="", xlim= c(-250,100))
  text(pca.geno.gl$scores[,1], pca.geno.gl$scores[,3] , labels = info.samples$`WORK NAME`, cex = 1, pos = 2, col = "black")
  title(main = "B", adj  = 0)
  title(xlab= "PC 1 (31%)", ylab="PC 3 (12%)", line=2, cex.lab=1.2)
  ## PC 2-3
  plot(pca.geno.gl$scores[,1], pca.geno.gl$scores[,2],
       cex=3, pch=20, col=alpha(my_cols_3[geno.gl.POP2$pop],0.6), main='K=2 PC1/PC2', xlab= "", ylab="", xlim= c(-250,150))
  text(pca.geno.gl$scores[,1], pca.geno.gl$scores[,2] , labels = info.samples$`WORK NAME`, cex = 1, pos = 2, col = "black")
  title(main = "C", adj  = 0)
  title(xlab="PC 1 (31%)", ylab="PC 2 (18%)", line=2, cex.lab=1.2)
  par(mai=c(0,0,0,0))
  plot.new()
  legend(x=0.18,y=0.9,legend=unique(geno.gl.POP$pop),col= my_cols_2[unique(geno.gl.POP$pop)],pch=20,ncol=4,cex=1.3,pt.cex=4,xpd=TRUE)
  legend(x=0.8,y=0.9,legend=unique(geno.gl.POP2$pop),col= my_cols_3[unique(geno.gl.POP2$pop)],pch=20,ncol=2,cex=1.3,pt.cex=4,xpd=TRUE)
}
dev.off()



pdf('POP.PCA.SCATTER.layout2.pdf')
{layout(matrix(c(1,2,3,4,4,4), ncol=2, nrow=2, byrow=TRUE), heights=c(2, 2))
  par(mai=rep(0.5, 4))
  ##  PC1-2 
  plot(pca.geno.gl$scores[,1], pca.geno.gl$scores[,2],
       cex=2, pch=20, col=alpha(my_cols_2[geno.gl.POP$pop],0.6), main='POP PC1/PC2',
       xlab= "", ylab="", xlim= c(-250,100))
  text(pca.geno.gl$scores[,1], pca.geno.gl$scores[,2] , labels = info.samples$`WORK NAME`, cex = 0.4, pos = 2, col = "black")
  title(xlab= "PC 1 (31%)", ylab="PC 2 (18%)", line=1.9, cex.lab=1.2)
  ## PC1-3
  plot(pca.geno.gl$scores[,1], pca.geno.gl$scores[,3],
       cex=2, pch=20, col=alpha(my_cols_2[geno.gl.POP$pop],0.6), main='POP PC1/PC3', xlab= "", ylab="", xlim= c(-250,100))
  text(pca.geno.gl$scores[,1], pca.geno.gl$scores[,3] , labels = info.samples$`WORK NAME`, cex = 0.4, pos = 2, col = "black")
  title(xlab= "PC 1 (31%)", ylab="PC 3 (12%)", line=1.9, cex.lab=1.2)
  ## PC 2-3
  plot(pca.geno.gl$scores[,1], pca.geno.gl$scores[,2],
       cex=2, pch=20, col=alpha(my_cols_3[geno.gl.POP2$pop],0.6), main='K=2 PC1/PC2', xlab= "", ylab="", xlim= c(-250,150))
  text(pca.geno.gl$scores[,1], pca.geno.gl$scores[,2] , labels = info.samples$`WORK NAME`, cex = 0.8, pos = 2, col = "black")
  title(xlab="PC 1 (31%)", ylab="PC 2 (18%)", line=2, cex.lab=1.2)
  plot(pca.geno.gl$scores[,1], pca.geno.gl$scores[,3],
       cex=2, pch=20, col=alpha(my_cols_3[geno.gl.POP2$pop],0.6), main='K=2 PC1/PC', xlab= "", ylab="", xlim= c(-250,150))
  text(pca.geno.gl$scores[,1], pca.geno.gl$scores[,3] , labels = info.samples$`WORK NAME`, cex = 0.8, pos = 2, col = "black")
  title(xlab="PC 1 (31%)", ylab="PC 3 (12%)", line=2, cex.lab=1.2)
  #par(mai=c(0,0,0,0))
  # plot.new()
  #legend(x = "topleft", legend = unique(geno.gl.POP$pop), pch=20, pt.cex=2 ,cex=1,col = my_cols_2[unique(geno.gl.POP$pop)])
}
dev.off()









#########################################
#####        FIS: inbreeding        #####
#####     Leishmania AETHIOPICA     #####
#########################################

###########################
### 1. data import/prep ###
###########################

geno_3 <- read.geno('esembled_Aeth.pop3.012')
pop_3 <- plot_data_K5$id[plot_data_K5$likely_assignment=="V3"]

###################
### 2. FUNCTION ###
###################
#function to calculate the FIS based on basic population genetic formulas
fis <- function(geno,pop) {
  p <- apply(geno[pop,],2,sum)/(2*length(pop))
  p2 <- p[which(p != 0 & p != 1)]
  Ho <- apply(geno[pop,], 2, function(x) sum(x==1))/length(pop)
  Ho <- Ho[which(p != 0 & p != 1)]
  He <- (2*p2*(1-p2))
  return(1-(Ho/He))}

####################################################
### 3. Calculating FIS for subpopulation V3: K=5 ###
####################################################
fis_AETH_3 <-fis(geno_3,pop_3)
fis_AETH_dataframe_3 <- as.data.frame(fis_AETH_3)
fis_AETH_dataframe_3$fis_AETH_3 <- as.numeric(fis_AETH_dataframe_3$fis_AETH_3)
fis_AETH_dataframe_3$names <- row.names(fis_AETH_dataframe_3)
write_xlsx(fis_AETH_dataframe_3,"FIS_AETH_3.xlsx")

mean(fis_AETH_dataframe_3$fis_AETH_3) 

####################
### 4. Plotting  ###
####################
pdf("Fis_Aeth_3.pdf")
hist(fis_AETH_dataframe_3$fis_AETH_3, xlab = "FIS", main = "FIS pop 3")
dev.off()


#########################################
#####        LDdecay plotting       #####
#####     Leishmania AETHIOPICA     #####
#########################################
esembled_Aeth.pop3 <- read.delim("esembled_Aeth.pop3.stat", header=FALSE, comment.char="#")

#Plotting with different limits
pdf("LDdecay_10000_pop3.pdf")
plot(esembled_Aeth.pop3$V2~esembled_Aeth.pop3$V1, type = "l",xlim = c(0, 10000), main= "LDdecay pop 3", xlab = "Distance (b)", ylab = "r^2" )
dev.off()

pdf("LDdecay_2000_pop3.pdf")
plot(esembled_Aeth.pop3$V2~esembled_Aeth.pop3$V1, type = "l",xlim = c(0, 2000), main= "LDdecay pop 3", xlab = "Distance (b)", ylab = "r^2" )
dev.off()

pdf("LDdecay_1000_pop3.pdf")
plot(esembled_Aeth.pop3$V2~esembled_Aeth.pop3$V1, type = "l",xlim = c(0, 1000), main= "LDdecay pop 3", xlab = "Distance (b)", ylab = "r^2" )
dev.off()



#########################################
#####        Depths to window       #####
#####     Leishmania AETHIOPICA     #####
#########################################
library(R.utils)
library(stringr)
library(gplots)
library(plyr)

###########################
### 1. data import/prep ###
###########################

files_10 <- list.files(path = 'depths_aeth', pattern = '_window_depths.txt', full.names = T)
files_10
myfilelist_10 <- lapply(files_10, read.table, header = T)
names(myfilelist_10) <- list.files(path = 'depths_aeth', pattern = '_window_depths.txt', full.names=FALSE)
#names(myfilelist_10) <- sub('depths_aeth',' ',files_10)

windows <- list()
extracted_cols <- lapply(myfilelist_10, function(x) x[,4]) #note the added comma!
windows <- do.call(cbind, extracted_cols) #no quotes needed around cbind
windows <- as.data.frame(windows)
colnames(windows) <- colnames(files_10)

#defining the matrix dims
depths <- matrix(ncol = 2+2*length(myfilelist_10), nrow=nrow(myfilelist_10$"103-83.2000_window_depths.txt"))
for (i in 1:length(myfilelist_10)) {
  depths[,2+i] <- myfilelist_10[[i]]$N_ZERO_DEPTH
  depths[,2+length(myfilelist_10)+i] <- myfilelist_10[[i]]$MEDIANHAPDEPTH
}

depths[,1] <- myfilelist_10$'103-83.2000_window_depths.txt'$CHR
depths[,2] <- myfilelist_10$'103-83.2000_window_depths.txt'$BIN
depths <- as.data.frame(depths)
colnames(depths) <- col_names
col_names <- c('chr', 'bin', '103-83','1123-81','117-82','130-83','1464-85','1561-87','169-83','32-83','678-82','68-83','85-83','GEREcl7','L100','L100cl1','L127','L86','LEM2357','LEM2358cl3','LEM3464','LEM3469','LEM3469cl1','LEM3469cl5','LEM3469cl7','LEM3469cl8','LEM3469cl9','LEM3497','LEM3498','WANDERA','103-83','1123-81','117-82','130-83','1464-85','1561-87','169-83','32-83','678-82','68-83','85-83','GEREcl7','L100','L100cl1','L127','L86','LEM2357','LEM2358cl3','LEM3464','LEM3469','LEM3469cl1','LEM3469cl5','LEM3469cl7','LEM3469cl8','LEM3469cl9','LEM3497','LEM3498','WANDERA')


depths[which(depths[,1]==i),2]/1000
depths$chr <- as.factor(depths$chr)
levels(depths$chr) <- c(1:36)
#levels(depths$chr)

#Creating a numerig dataframe
depths$chr <- as.numeric(depths$chr)
depths$bin <- as.numeric(depths$bin)
for (i in 1:58)
{depths[,i] <- as.numeric(depths[,i])}

#summary(depths)


####################
### 2. Plotting  ###
####################
for (i in 1:36) {
  pdf(paste('N_ZERO_DEPTH', i,'.pdf',sep=''), useDingbats = F, width = 15)
  lab <- depths[which(depths[,1]==i),2]/1000
  lab2 <- as.character(lab)
  heatmap.2(t(depths[which(depths[,1]==i),c(names(depths)[3:30])]), 
            Colv = F, dendrogram= 'row', trace="none", symkey = F, 
            scale = 'none', cexRow = 0.5, cexCol = 0.5, margins = c(2,4), 
            labCol = lab2, labRow = names(depths)[3:30],  
            col=rich.colors(50), offsetCol=0, density.info = 'none')
  dev.off()
}


depths[,3:30]<- NULL
for (i in 1:36) {
  pdf(paste('hapwindowdepths', i,'.pdf',sep=''), useDingbats = F, width = 15)
  lab <- depths[which(depths[,1]==i),2]/1000
  lab2 <- as.character(lab)
  heatmap.2(t(depths[which(depths[,1]==i),c(names(depths)[31:58])]), 
            Colv = F, dendrogram= 'row', trace="none", symkey = F, 
            scale = 'none', cexRow = 0.5, cexCol = 0.5, margins = c(2,4), 
            labCol = lab2, labRow = names(depths)[31:58],  
            col=rich.colors(50), offsetCol=0, density.info = 'none')
  dev.off()
}


#########################################
#####     Somy heatmap per chr      #####
#####     Leishmania AETHIOPICA     #####
#########################################

#install.packages("R.utils")
library(R.utils)
library(stringr)
#install.packages("gplots")
library(gplots)
library(plyr)
library(heatmap3)

###########################
### 1. data import/prep ###
###########################
files <- list.files(path = 'depths_aeth', pattern = 'chr_depths.txt', full.names = T)
myfilelist <- lapply(files, read.table, header = T)
names(myfilelist) <- list.files(path = 'depths_aeth', pattern = 'chr_depths.txt', full.names=FALSE)
names <-  list.files(path = 'depths_aeth', pattern = 'chr_depths.txt', full.names=FALSE)
names

chromosomes <- myfilelist$`103-83.2000_chr_depths.txt`$CHR
chromosomes <- as.list(chromosomes)
#class(chromosomes)

somy <- list()
extracted_cols <- lapply(myfilelist, function(x) x[,2]) #note the added comma!
somy <- do.call(cbind, extracted_cols) #no quotes needed around cbind
somy <- as.data.frame(somy)
rownames(somy) <- chromosomes
col_names_somy <- c('103-83','1123-81','117-82','130-83','1464-85','1561-87','169-83','32-83','678-82','68-83','85-83','GEREcl7','L100','L100cl1','L127','L86','LEM2357','LEM2358cl3','LEM3464','LEM3469','LEM3469cl1','LEM3469cl5','LEM3469cl7','LEM3469cl8','LEM3469cl9','LEM3497','LEM3498','WANDERA')
#dim(somy)

somy2 <- somy
#Calculating the median SOMY per isolate
for (i in 1:29) {somy2[,i] <- somy[,i]/median[,i]}
colnames(somy2) <- col_names_somy
rownames(somy2) <- c(1:36)

####################
### 2. Plotting  ###
####################

pdf('SOMY.pdf', useDingbats = F, width = , height = )
heatmap.2(as.matrix(t(somy2)), trace = 'none', col=rich.colors(50), dendrogram = 'row', srtCol = 90,
          Colv = F, density.info = 'none', offsetCol=0, cexRow = 0.5, cexCol = 1, margins = c(2,4))
dev.off()
