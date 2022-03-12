for i in 0 50000 100000 150000 200000 250000 300000 350000 400000 450000 500000 ;do
  cat outcome_icd_codes.txt  | tail -n+2 | cut -f1 | while read line;do
     y=${line,,}
     if [ ! -e raw_output/diag.coding.${y}.${i}.txt.gz ]; then
       echo diag.coding.${y}.${i}.txt.gz
       #python jc_pheno.py $i $line &
       #sleep 200
     fi
  done
done
