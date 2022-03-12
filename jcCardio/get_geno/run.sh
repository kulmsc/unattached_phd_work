
for i in {15..17};do
  cat ~/athena/ukbiobank/setup_exome/ukb_fe_exm_chrall_v1.bim | awk -v var="$i" '$1 == var {print $2}' > good_snp.$i

  plink --bfile ~/athena/ukbiobank/exome/ukbb.exome.${i} --keep-fam ../poster_jc/temp --extract good_snp.$i --recode vcf --out into.${i}

  #vep -i into.vcf --cache --dir_cache /home/kulmsc/athena/refs/vep_cache --force_overwrite --plugin LoF,loftee_path:/home/kulmsc/athena/refs/new_loftee,human_ancestor_fa:/home/kulmsc/athena/refs/for_loftee/human_ancestor.fa.gz,conservation_file:/home/kulmsc/athena/refs/for_loftee/phylocsf_gerp.sql,gerp_bigwig:/home/kulmsc/refs/for_loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw --refseq --symbol --use_given_ref -o anno.${i}.vcf

  rm anno.${i}.vcf.gz

  vep -i into.${i}.vcf --cache --dir_cache /home/kulmsc/athena/refs/vep_cache --force_overwrite --plugin LoF,loftee_path:/home/kulmsc/athena/refs/new_loftee/,human_ancestor_fa:/home/kulmsc/athena/refs/for_loftee/human_ancestor.fa.gz --refseq --symbol --use_given_ref -o anno.${i}.vcf

  gzip anno.${i}.vcf
  
  rm into.${i}.vcf
  rm good_snp.$i
  #vep -i into.vcf --format vcf --custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN --refseq --symbol --use_given_ref --force_overwrite -o clinvar.${i}.vcf

  #gzip -f clinvar.${i}.vcf

done
