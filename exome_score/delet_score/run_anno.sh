
#java -Xmx20G -jar ~/Programs/snpEff/snpEff.jar build -v -onlyReg GRCh37.75 #from release-75, although must move to the local Programs folder


chr=$1

#Set up exome
plink --bfile ~/athena/ukbiobank/exome/ukbb.exome.${chr} --keep-fam single_fam --recode vcf --out into.${chr}

cat into.${chr}.vcf | fgrep -v "#" | awk '{print "chr"$1 "\t" $2 "\t" ($2+1) "\t" $3}' > into.${chr}.bed

~/Programs/liftOver into.${chr}.bed ~/athena/refs/for_picard/hg38ToHg19.over.chain.gz out.${chr}.bed unMapped

rm unMapped into.${chr}.bed

head -100 into.${chr}.vcf | fgrep '#' > top
Rscript change_vcf.R ${chr}
cat top temp > out.${chr}.vcf #This is almost working, the snp id still has the wrong position

bcftools sort out.${chr}.vcf > temp; mv temp out.${chr}.vcf

rm into.${chr}.vcf top temp

#set up impute
plink2 --bgen ~/athena/ukbiobank/imputed/ukbb.${chr}.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${chr}.sample --keep-fam single_fam --make-bed --out temp
plink --bfile temp --recode vcf --out impute.${chr}
rm temp*

################# General Annotations #################
java -Xmx8g -jar ~/Programs/snpEff/snpEff.jar -v GRCh37.75 out.${chr}.vcf > anno_output/complete.chr${chr}.ann.vcf
gzip -f anno_output/complete.chr${chr}.ann.vcf

java -Xmx8g -jar ~/Programs/snpEff/snpEff.jar -v GRCh37.75 impute.${chr}.vcf > anno_output/impute.chr${chr}.ann.vcf
gzip -f anno_output/complete.chr${chr}.ann.vcf




#java -Xmx8g -jar ~/Programs/snpEff/SnpSift.jar phastCons ./phastCons out.${chr}.vcf > anno_output/complete.chr${chr}.phastcon.vcf
#gzip anno_output/complete.chr${chr}.phastcon.vcf


#################### Delet val annotations #######################
java -Xmx8g -jar ~/Programs/snpEff/SnpSift.jar dbnsfp -v -db ./dbNSFP4.1a.txt.gz out.${chr}.vcf > anno_output/complete.chr${chr}.dbnsfp.vcf
gzip -f anno_output/complete.chr${chr}.dbnsfp.vcf

java -Xmx8g -jar ~/Programs/snpEff/SnpSift.jar dbnsfp -v -db ./dbNSFP4.1a.txt.gz impute.${chr}.vcf > anno_output/impute.chr${chr}.dbnsfp.vcf
gzip -f anno_output/impute.chr${chr}.dbnsfp.vcf


#ls gene_set/*set | while read line;do
#  outname=`echo $line | cut -f2 -d"/" | cut -f1 -d'.'`
#  java -Xmx8g -jar ~/Programs/snpEff/SnpSift.jar geneSets -v ${line} anno_output/complete.chr${chr}.ann.vcf > anno_output/complete.chr${chr}.${outname}.vcf
#  gzip anno_output/complete.chr${chr}.${outname}.vcf
#done

#ls bb_files | while read line;do
#  bb_short=`echo ${line} | cut -f1 -d'.'`
#  java -Xmx8g -jar ~/Programs/snpEff/snpEff.jar -v -interval bb_files/${line} GRCh37.75 out.${chr}.vcf > anno_output/complete.chr${chr}.${bb_short}.vcf
#  gzip anno_output/complete.chr${chr}.${bb_short}.vcf
#done


rm out.${chr}.vcf
rm impute.${chr}.vcf

#ls ~/Programs/snpEff/data/GRCh37.75/ | fgrep regulation | fgrep gff | fgrep -v gz | cut -f2 -d'_' | cut -f1 -d'.' | while read line;do
#  java -Xmx8g -jar ~/Programs/snpEff/snpEff.jar -v -reg CD4 -reg GM06990 -reg GM12878 -reg H1ESC -reg HeLa-S3 -reg HepG2 -reg HMEC -reg HSMM -reg HUVEC -reg IMR90 -reg K562 -reg NH-A -reg NHEK GRCh37.75 out.22.vcf > anno_output/complete.chr22.reg.vcf
#done


