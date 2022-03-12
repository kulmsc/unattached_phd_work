chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}


i=1
cat all_specs/tweedie_param_specs | tail -n +2 | while read spec;do
  plim=`echo $spec | cut -f1 -d' '`

  if [ ! -e ~/athena/xwas/mod_sets/${author}/${low_author}.tweedie.${i}.ss ]; then

    plink --memory 4000 --threads 1 --bfile geno_files/${low_author}.${chr} --clump temp_files/ss.${low_author}.${chr} --clump-snp-field RSID --clump-p1 $plim --clump-r2 0.25 --out ${d}/out

    if [ -f ${d}/out.clumped ]; then

      plink --bfile geno_files/${low_author}.${chr} --freq --out ${d}/freq

      Rscript helper_scripts/add_tweedie.R $author $chr $d

      Rscript helper_scripts/tweedie.R $author $chr $d $i

    fi

  fi

  let i=i+3

done
