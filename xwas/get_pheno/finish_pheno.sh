#ls raw_output/ | cut -f3 -d'.' | sort | uniq | while read author;do
author=lung
  rm pheno_defs/diag.${author}.txt.gz
  rm pheno_defs/time.${author}.txt.gz
  ls raw_output/diag.coding.${author}* | cut -f2 -d'/' | cut -f4 -d'.' | sort -n | while read num;do
    zcat raw_output/diag.coding.${author}.${num}.txt.gz >> pheno_defs/diag.${author}.txt
    zcat raw_output/diag.time.${author}.${num}.txt.gz >> pheno_defs/time.${author}.txt
  done
  gzip pheno_defs/diag.${author}.txt
  gzip pheno_defs/time.${author}.txt
#done
