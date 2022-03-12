#! /bin/bash -l
 
#SBATCH --partition=panda_physbio   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=rgwas
#SBATCH --time=72:00:00   # HH/MM/SS
#SBATCH --mem=14G   # memory requested, units available: K,M,G,T
 
source ~/.bashrc
 
 


#rm -f list_beds.txt
#for chr in {1..22}; do
#  echo "/home/kulmsc/athena/ukbiobank/calls/ukbb.${chr}" >> list_beds.txt
#done

#plink --merge-list list_beds.txt --make-bed --out ukb_cal_allChrs


#regenie \
#  --step 1 \
#  --bed ukb_cal_allChrs \
#  --remove /home/kulmsc/athena/ukbiobank/custom_qc/impute_bad/array.all.fam \
#  --exclude /home/kulmsc/athena/ukbiobank/custom_qc/impute_bad/array.all.snps \
#  --phenoFile use_pheno \
#  --covarFile use_covar \
#  --bt \
#  --bsize 1000 \
#  --lowmem \
#  --lowmem-prefix tmpdir/regenie_tmp_preds \
#  --out ukb_step1_BT

for i in {1..22};do
regenie \
  --step 2 \
  --bgen /home/kulmsc/athena/ukbiobank/imputed/ukbb.${i}.bgen \
  --ref-first \
  --sample /home/kulmsc/athena/ukbiobank/imputed/ukbb.${i}.sample \
  --phenoFile use_pheno \
  --covarFile use_covar \
  --bt \
  --firth 0.01 --approx \
  --pred ukb_step1_BT_pred.list \
  --bsize 400 \
  --split \
  --out ukb_step2_BT_chr${i}
done


