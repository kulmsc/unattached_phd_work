#anno.22.vcf.gz comes from LOFTEE


for i in {1..22};do

  #determine the positions that align with an important transcript and HC variants
  zcat anno.${i}.vcf.gz | fgrep -w -f simple_transcripts.txt | fgrep HC > potent_data
  cat potent_data | cut -f14 | cut -f3 -d';' | cut -f2 -d'=' > gene.${i}
  cat potent_data | cut -f1 > id.${i}
  cat potent_data | cut -f5 > curr_transcript.${i}


  #get all data for these indiviudals and convert to raw
  plink --bfile ~/athena/ukbiobank/exome/ukbb.exome.${i} --extract id.${i} --recode A --out curr.${i}
  gzip -f curr.${i}.raw

done


#just getting the data from clinvar for positions called special by loftee
for i in {1..22};do

  #determine the positions that align with an important transcript and HC variants
  zcat clinvar.${i}.vcf.gz | fgrep -w -f id.${i} | cut -f14 | cut -d';' -f5 | cut -f2 -d'=' > clin_data.${i}
  zcat clinvar.${i}.vcf.gz | fgrep -w -f id.${i} | cut -f14 | cut -f3 -d';' > clin_cldn.${i}
done




#################################################################################



#getting files for all variants in special genes
#for i in {1..22};do
#  zcat anno.${i}.vcf.gz | fgrep -w -f simple_genes.txt > potent_data
#  cat potent_data | cut -f14 | tr ';' '\n' | fgrep 'SYMBOL=' | cut -f2 -d'=' > all_gene.${i}
#  cat potent_data | cut -f1 > all_id.${i}

#  plink --bfile ~/athena/ukbiobank/exome/ukbb.exome.${i} --extract all_id.${i} --recode A --out all_curr.${i}
#  gzip all_curr.${i}.raw

#done


#just getting the data from clinvar for specific transcripts
for i in {1..22};do

  #zcat anno.${i}.vcf.gz | fgrep -w -f simple_genes.txt | fgrep HC | cut -f1 > loftee_hc_all_data.${i}

  zcat anno.${i}.vcf.gz | fgrep -w -f simple_genes.txt | cut -f1 > temp
  zcat clinvar.${i}.vcf.gz | fgrep -w -f temp | cut -f1 > clin_id.${i}
  #zcat clinvar.${i}.vcf.gz | fgrep -w -f temp | cut -f14 | cut -d';' -f5 | cut -f2 -d'=' > clin_all_data.${i}
  #zcat clinvar.${i}.vcf.gz | fgrep -w -f temp | cut -f14 | cut -f3 -d';' > clin_all_cldn.${i}

done




#then in setup_mat read into
