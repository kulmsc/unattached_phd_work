ls 2*vcf.gz | while read line;do
  index=`echo $line | cut -f1 -d'.'`
  zcat $line > temp.vcf
  bcftools view temp.vcf -Oz -o file.${index}.vcf.gz
  bcftools index file.${index}.vcf.gz
  rm temp.vcf  
done

bcftools merge fil*vcf.gz --missing-to-ref -Oz -o all_merged.vcf.gz

ls 2*vcf.gz | while read line;do
  id=`echo $line | cut -f1 -d'.'`
  vcftools --gzvcf $line --site-mean-depth --out ${id}
done

python faster_compile_depths.py

Rscript compile_depths.R

plink2 --pfile all_merged --extract bed1 bed_extract.txt --make-pgen --out all_qc --allow-extra-chr

plink2 --pfile all_qc --export vcf-4.2 --out all_qc

#plan is to convert get read depth for all individuals, where we have mean known read depth
#above 30 (or so for 3? or more people) keep the variant if not discard (may need to do this
#part within PLINK, then convert back to vcf and continute with SnpEff

#plink2 --vcf all_merged.vcf.gz --make-pgen --out all_merged --allow-extra-chr

#plink --vcf all_merged.vcf.gz --biallelic-only --make-bed --out all_merged --allow-extra-chr
