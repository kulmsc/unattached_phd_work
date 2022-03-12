
ls many_mod_sets | while read line;do
  chromo=`echo $line | cut -f3 -d'.'`
  new=`echo $line | cut -f1,2,3 -d'.'`
  maybeset=`echo $line | cut -f4 -d'.'`
  if [ $maybeset == set ];then
    echo skip
  else
    new=`echo $line | cut -f1,2,3,4 -d'.'`
  fi
  if [ ! -e small_scores/${new}.profile.zst ];then
    plink  --memory 12000 --threads 12 --bfile ~/athena/ukbiobank/exome/ukbb.exome.${chromo} --keep-allele-order --score many_mod_sets/${line} 5 3 6 sum --out small_scores/${new}
    zstd --rm small_scores/${new}.profile
  fi
done
