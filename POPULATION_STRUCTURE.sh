#ADMIXTURE:
##Preparing data files for admixture (vcf to bed to bim file)
plink --vcf esembled_Aeth.18.NOL100cl1.NO678_82.vcf.gz --make-bed --out  esembled_Aeth.18.NOL100cl1.NO678_82 --allow-extra-chr
awk '{$1=0;print $0}'  esembled_Aeth.18.NOL100cl1.NO678_82.bim >  esembled_Aeth.18.NOL100cl1.NO678_82.bim.tmp
mv  esembled_Aeth.18.NOL100cl1.NO678_82.bim.tmp  esembled_Aeth.18.NOL100cl1.NO678_82.bim

##admixture for K1:10
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv  esembled_Aeth.18.NOL100cl1.NO678_82.bed $K | tee log${K}_NODUP.out; done

##txt file with the CV errors listed for K1:10
grep -h CV log*_NODUP.out > CV_NODUP.txt


#LDdecay:
##Preparing data file for V3(K=5)
bcftools view -a -s ^103-83,1123-81,117-82,130-83,1464-85,85-83,GEREcl7,L100,L127,LEM2358cl3,LEM3464,WANDERA esembled_Aeth.18.NOL100cl1.NO678_82.vcf.gz | bcftools view -e 'ALT=="."' - | bgzip > esembled_Aeth.pop3.vcf.gz   

##PopLDdecay
PopLDdecay -InVCF esembled_Aeth.pop3.vcf.gz -OutStat esembled_Aeth.pop3

##Fist plotting LDdecay
perl Plot_OnePop.pl -inFile esembled_Aeth.pop3.stat.gz -output Fig_pop3    
