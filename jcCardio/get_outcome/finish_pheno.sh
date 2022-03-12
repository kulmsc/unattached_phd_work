author=$1
rm pheno_defs/diag.${author}.txt.gz
rm pheno_defs/time.${author}.txt.gz
ls raw_output/diag.coding.${author}* | cut -f2 -d'/' | cut -f4 -d'.' | sort -n | while read num;do
  zcat raw_output/diag.coding.${author}.${num}.txt.gz >> pheno_defs/diag.${author}.txt
  zcat raw_output/diag.time.${author}.${num}.txt.gz >> pheno_defs/time.${author}.txt
done
gzip -f pheno_defs/diag.${author}.txt
gzip -f pheno_defs/time.${author}.txt
