for eth in tot eas;do

#for chr in {1..22};do
for chr in 5;do
  

  #len=`zcat ready_files/${eth}.${chr}.annot.gz  | wc -l`
  #if [ $len == 1 ];then

  zcat baseline_annot/baseline.${chr}.annot.gz | fgrep -w -f hapmap_rsids > curr_annot
  cat curr_annot | cut -f3 > curr_rsids
  #cat curr_rsids | fgrep -v rs12186596 > x;mv x curr_rsids

  plink --bfile ~/athena/refs/1000genomes/${eth}.${chr} --extract curr_rsids --make-bed --out current

  plink --bfile current --freq --out freq_files/${eth}.${chr}

  cat current.bim | cut -f2 > temp
  zcat baseline_annot/baseline.${chr}.annot.gz | fgrep -w -f temp > curr_annot

  cat curr_annot | cut -f1 > temp
  cat curr_annot | cut -f3-100 > temp2
  cat current.bim | cut -f4 > temp3
  zcat baseline_annot/baseline.${chr}.annot.gz | head -1 > curr_annot
  paste temp temp3 temp2 >> curr_annot

  cat current.bim | cut -f2 > regression.snplist

  taskset -c 40-44 ldsc\
    --bfile current\
    --extract regression.snplist\
    --ld-wind-cm 1\
    --out weight_files/w.${eth}.${chr}


  taskset -c 40-44 ldsc\
   --l2\
   --bfile current\
   --ld-wind-cm 1\
   --annot curr_annot\
   --out baseline_scores/${eth}.${chr}


  mv curr_annot ${eth}.${chr}.annot
  gzip ${eth}.${chr}.annot

  rm temp temp2 temp3 curr_annot current.bim current.bed current.fam curr_rsids

  #fi

done

done
