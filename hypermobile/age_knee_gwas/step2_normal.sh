#! /bin/bash -l

#SBATCH --partition=panda_physbio   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=rgwas
#SBATCH --time=72:00:00   # HH/MM/SS
#SBATCH --mem=16G   # memory requested, units available: K,M,G,T




for i in {1..22};do
#i=$1
taskset -c 10-18 regenie \
  --step 2 \
  --bgen /home/kulmsc/athena/ukbiobank/imputed/ukbb.${i}.bgen \
  --ref-first \
  --sample /home/kulmsc/athena/ukbiobank/imputed/ukbb.${i}.sample \
  --remove /home/kulmsc/athena/ukbiobank/custom_qc/impute_bad/impute.${i}.fam \
  --extract ../knee_rep_gwas/want_knee_ids \
  --phenoFile use_pheno_simple_${1} \
  --covarFile use_covar_simple_${1} \
  --bt \
  --firth 0.01 --approx \
  --pred ukb_step1_agelocal_${1}_pred.list \
  --bsize 400 \
  --split \
  --out res/ukb_step2_${1}_full_covar_chr${i}
  #sleep 300
done

