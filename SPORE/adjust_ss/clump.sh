chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}


i=1
cat all_specs/clump_param_specs | tail -n +2 | while read spec;do
  plim=`echo $spec | cut -f1 -d' '`
  r2lim=`echo $spec | cut -f2 -d' '`

  #if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.${chr}.clump.${i}.ss ]; then

    plink --memory 4000 --threads 1 --bfile geno_files/${low_author}.${chr} --clump temp_files/ss.${low_author}.${chr} --clump-snp-field RSID --clump-p1 $plim --clump-r2 $r2lim --out ${d}/out

    if [ -f ${d}/out.clumped ]; then
      sed -e 's/ [ ]*/\t/g' ${d}/out.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > ${d}/done_rsids
      fgrep -w -f ${d}/done_rsids temp_files/ss.${low_author}.${chr} > ~/athena/SPORE/mod_sets/${author}/${low_author}.${chr}.clump.${i}.ss
    fi
  #fi

  let i=i+1
done
