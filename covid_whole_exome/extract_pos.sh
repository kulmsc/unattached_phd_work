rm all_pos

cat genes_to_investigate | cut -f1 | while read gene;do

  cat vcf_files/all_qc.ann.vcf | fgrep -w $gene | cut -f1,2 > temp
  len=`cat temp | wc -l`
  yes $gene | head -$len > temp2
  paste temp temp2 >> all_pos

done
rm temp temp2
