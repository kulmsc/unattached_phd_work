
for i in {1..22};do
regenie \
  --step 2 \
  --bgen /home/kulmsc/athena/ukbiobank/imputed/ukbb.${i}.bgen \
  --ref-first \
  --sample /home/kulmsc/athena/ukbiobank/imputed/ukbb.${i}.sample \
  --remove /home/kulmsc/athena/ukbiobank/custom_qc/impute_bad/impute.${i}.fam \
  --exclude /home/kulmsc/athena/ukbiobank/custom_qc/impute_bad/impute.${i}.snps \
  --phenoFile new_diabetes_hypo_pheno \
  --covarFile use_diabetes_hypo_covar \
  --bt \
  --firth 0.01 --approx \
  --pred ukb_step1_BT_pred.list \
  --bsize 400 \
  --split \
  --out ukb_step2_wdiabeteshypo_chr${i}
done

