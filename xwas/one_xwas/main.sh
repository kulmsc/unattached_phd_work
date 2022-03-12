
pheno_name=cd

Rscript prep_pheno.R $pheno_name

Rscript setup_match.R

cat ~/athena/ukbiobank/qc/imputed/ukb_mfi_chrX_v3.txt | awk '$6 > 0.001 && $8 > 0.1 {print $2}' > snps

plink2 --threads 10 --memory 16000 --bgen ~/athena/ukbiobank/imputed/ukbb.23.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.23.sample --extract snps --keep-fam use_eid --geno 0.1 --hwe 1e-50 midp --maf 0.001 --make-bpgen --out ukbb.23

plink2 --bpfile ukbb.23 --keep-males --freq --out male
plink2 --bpfile ukbb.23 --keep-males --missing --out male
plink2 --bpfile ukbb.23 --keep-females --freq --out female
plink2 --bpfile ukbb.23 --keep-females --missing --out female

Rscript extra_qc.R

rm male*
rm female*

plink2 --threads 10 --memory 16000 --bpfile ukbb.23 --exclude bad_rsid --make-bed --out toxwas

rm ukbb.23*

Rscript fix_fam.R

cat toxwas.bim | cut -f2 > just_snps

i=1
split -l 50000 -d just_snps split_up_snps
rm just_snps
rm snps

ls split_up_snps* | while read line;do

  plink --bfile toxwas --extract $line --make-bed --out current_use

  runxwas --bfile current_use --xwas --strat-sex --fishers --covar no_sex_covars --ci 0.95 --out ${pheno_name}_collect_res/full.mod.${i}

  runxwas --bfile current_use --logistic --covar no_sex_covars --ci 0.95 --out ${pheno_name}_collect_res/logistic.1.mod.${i}

  runxwas --bfile current_use --xchr-model 2 --logistic --covar no_sex_covars --ci 0.95 --out ${pheno_name}_collect_res/logistic.2.mod.${i}

  let i=i+1
done

rm split_up_snps*
rm current_use*
rm toxwas*
