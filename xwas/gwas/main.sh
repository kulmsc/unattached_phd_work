
pheno_name=cd

Rscript prep_pheno.R $pheno_name

Rscript setup_match.R

for i in {1..22};do

  cat ~/athena/ukbiobank/qc/imputed/ukb_mfi_chr${i}_v3.txt | awk '$6 > 0.001 && $8 > 0.1 {print $2}' > snps

  plink2 --threads 10 --memory 16000 --bgen ~/athena/ukbiobank/imputed/ukbb.${i}.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${i}.sample --extract snps --keep-fam use_eid  --geno 0.1 --hwe 1e-50 midp --maf 0.001 --make-bpgen --out ukbb.${i}

  plink2 --threads 10 --memory 16000 --bpfile ukbb.${i} --pheno use_pheno --covar covars --logistic hide-covar --covar-variance-standardize --out results/${pheno_name}.${i}

  rm ukbb.${i}.*

done



