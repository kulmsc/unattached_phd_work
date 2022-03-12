
#rm -f list_beds.txt
#for chr in {1..22}; do
#  echo "/home/kulmsc/athena/ukbiobank/calls/ukbb.${chr}" >> list_beds.txt
#done

#plink --merge-list list_beds.txt --make-bed --out ukb_cal_allChrs


regenie \
  --step 1 \
  --bed ukb_cal_allChrs \
  --remove /home/kulmsc/athena/ukbiobank/custom_qc/impute_bad/array.all.fam \
  --exclude /home/kulmsc/athena/ukbiobank/custom_qc/impute_bad/array.all.snps \
  --phenoFile new_diabetes_hypo_pheno \
  --covarFile use_diabetes_hypo_covar \
  --bt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmpdir/regenie_tmp_preds \
  --out ukb_step1_BT
