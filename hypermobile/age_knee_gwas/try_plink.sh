
for i in {1..22};do
  bgenix -g  /home/kulmsc/athena/ukbiobank/imputed/ukbb.${i}.bgen -incl-rsids ../knee_rep_gwas/want_knee_ids > temp.bgen
  plink2 --bgen temp.bgen ref-first --sample /home/kulmsc/athena/ukbiobank/imputed/ukbb.${i}.sample --make-bed --out small_files/use.${i}
  if [ -e small_files/use.${i}.bim ]; then
    cat small_files/use.${i}.bim | cut -f2 | sort | uniq -d > bad_ids
    len=`cat bad_ids | wc -l`
    if [ $len -gt 0 ];then
      plink --bfile small_files/use.${i} --exclude bad_ids --make-bed --out small_files/temp
      mv small_files/temp.bed small_files/use.${i}.bed
      mv small_files/temp.bim small_files/use.${i}.bim
      mv small_files/temp.fam small_files/use.${i}.fam
    fi
  fi
  rm bad_ids
done



ls small_files/*bed | cut -f1-2 -d'.' > merge_list
plink --merge-list merge_list --make-bed --out small_files/total


plink2 --bfile small_files/total --pheno use_pheno_simple --1 --covar use_covar_simple --glm interaction --parameters 1,2-15,26 --out small_files/res_age_inter

plink2 --bfile small_files/total --pheno use_pheno_simple --1 --covar use_covar_simple --glm interaction --parameters 1,2-15,28 --out small_files/res_bmi_inter

#Additive effect ('ADD')
#Dominance deviation ('DOMDEV')
#First covariate ('COVAR1' if not named in the file)
#Second covariate
#Dosage-first covariate interaction ('ADDxCOVAR1')
#Dominance deviation-first covariate interaction ('DOMDEVxCOVAR1')
#Dosage-second covariate interaction ('ADDxCOVAR2')
#Dominance deviation-second covariate interaction ('DOMDEVxCOVAR2')

#1,3-17,28
#14 total covars
