#ac
#primary_icd_term=M750
#primary_opcs_term="XXX"
#other_terms="E02\|E03\|E05\|E10\|E11\|S42\|S43\|S46"

#knee
primary_icd_term="XXX"
primary_opcs_term="W401\|W411\|W422"
other_terms="M17\|M15\|M05\|M06\|M87\|M2556"

#hip
#primary_icd_term="M160\|M161\|M162\|M163\|M169"
#primary_opcs_term="W371\|W381\|W391\|W931\|W941\|W951"
#other_terms="M16\|M15\|M05\|M06\|M87\|M2555"

cat ~/athena/ukbiobank/hesin/hesin_diag.txt | grep $primary_icd_term | cut -f1 | sort | uniq > primary_eid
cat ~/athena/ukbiobank/hesin/hesin_oper.txt | grep $primary_opcs_term | cut -f1 | sort | uniq >> primary_eid

Rscript check_hesin.R


cat ~/athena/ukbiobank/hesin/hesin_diag.txt | grep $other_terms > subset_diag

cat subset_diag | fgrep -v -w -f primary_eid  > never_primary_diag_but_covar
cat subset_diag | fgrep -w -f primary_eid  > use_diag


cat kr_risk_icd | while read line;do
  code=`echo $line | tr '\t' ' ' | cut -f1 -d' '`
  name=`echo $line | tr '\t' ' ' | cut -f2 -d' '`
  code=`echo $code | sed -e 's/^"//' -e 's/"$//'`
  echo $code | tr '|' '\n' | while read x;do
    echo $x 
    echo $name
    cat never_primary_diag_but_covar | grep $x | cut -f1 | sort | uniq >> kr_covar.${name}.eid
    cat kr_covar.${name}.eid | wc -l
    echo
  done
done

ls kr_covar.* | while read line;do
  cat $line | sort | uniq > temp
  mv temp $line
done

rm never_primary_diag_but_covar
rm subset_diag
rm primary_eid
rm use_diag
