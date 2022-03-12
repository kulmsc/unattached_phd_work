  
filename=gene.5.vcf
cat $filename | fgrep ExcessHet | while read variant;do

  basic_info=`echo $variant | cut -f1-5 -d' '`
  case_info=`echo $variant | cut -f8 -d' ' | tr ';' '\n' | fgrep Cases | cut -f2 -d'='`
  control_info=`echo $variant | cut -f8 -d' ' | tr ';' '\n' | fgrep Controls | cut -f2 -d'='`
  all_info=`echo $variant | cut -f8 -d' ' | tr ';' '\n' | fgrep ANN | tr '|' '\n' | head -4 | tr '\n' '\t'`
  cc_all=`echo $variant | cut -f8 -d' ' | tr ';' '\n' | fgrep CC_ALL | cut -f2 -d'='`
  cc_geno=`echo $variant | cut -f8 -d' ' | tr ';' '\n' | fgrep CC_GENO | cut -f2 -d'='`
  cc_dom=`echo $variant | cut -f8 -d' ' | tr ';' '\n' | fgrep CC_DOM | cut -f2 -d'='`
  cc_rec=`echo $variant | cut -f8 -d' ' | tr ';' '\n' | fgrep CC_REC | cut -f2 -d'='`

  echo $basic_info $case_info $control_info $cc_all $cc_geno $cc_dom $cc_rec $all_info
done

