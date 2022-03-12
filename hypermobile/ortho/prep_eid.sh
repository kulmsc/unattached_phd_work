i=1
cat icd_codes.txt | while read line;do
  cat ~/athena/ukbiobank/hesin/hesin_diag.txt | fgrep $line | cut -f1 | sort | uniq > eid_gwas_$i
  let i=i+1
done
