cat fresh_jc_codes.txt | tail -n +3 | while read line;do python jc_pheno.py 500000 $line; done

