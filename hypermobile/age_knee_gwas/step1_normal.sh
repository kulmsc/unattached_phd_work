
#rm -f list_beds.txt
#for chr in {1..22}; do
#  echo "/home/kulmsc/athena/ukbiobank/calls/ukbb.${chr}" >> list_beds.txt
#done

#plink --merge-list list_beds.txt --make-bed --out ukb_cal_allChrs


taskset -c 30-38 regenie \
  --step 1 \
  --bed ../adhes_cap_gwas/ukb_cal_allChrs \
  --remove /home/kulmsc/athena/ukbiobank/custom_qc/impute_bad/array.all.fam \
  --exclude /home/kulmsc/athena/ukbiobank/custom_qc/impute_bad/array.all.snps \
  --phenoFile use_pheno_simple_${1} \
  --covarFile use_covar_simple_${1} \
  --bt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmpdir/regenie_local_preds \
  --out ukb_step1_agelocal_${1}

