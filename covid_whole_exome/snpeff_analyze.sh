java -jar ~/bin/snpEff.jar download -v GRCh37.75
java -Xmx4g -jar ~/bin/snpEff.jar -v -stats ex2.html GRCh37.75 vcf_files/all_qc.vcf.gz > vcf_files/all_qc.ann.vcf
java -Xmx1g -jar ~/bin/SnpSift.jar caseControl -v -tfam vcf_files/pheno_for_gwas vcf_files/all_qc.ann.vcf > vcf_files/all_qc.cc.vcf


cat vcf_files/all_qc.cc.vcf | java -Xmx1g -jar ~/bin/SnpSift.jar filter "CC_ALL < 0.0001 & ((ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE'))" > candidate/snpeff.1.vcf

cat vcf_files/all_qc.cc.vcf | java -Xmx1g -jar ~/bin/SnpSift.jar filter "CC_DOM < 0.0001 & ((ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE'))" > candidate/snpeff.2.vcf

cat vcf_files/all_qc.cc.vcf | java -Xmx1g -jar ~/bin/SnpSift.jar filter "CC_REC < 0.0001 & ((ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE'))" > candidate/snpeff.3.vcf

cat vcf_files/all_qc.cc.vcf | java -Xmx1g -jar ~/bin/SnpSift.jar filter "CC_GENO < 0.0001 & ((ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE'))" > candidate/snpeff.4.vcf


