chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}

#https://github.com/YangLabHKUST/XPASS

cat all_specs/xpass_param_specs | while read spec;do
  aux_pop=`echo $spec | cut -f1 -d' '`
  index=`echo $spec | cut -f2 -d' '`  

  if [ $author != $aux_pop ];then
  #aux_pop=european

  zcat ../raw_ss/${aux_pop}/chr_ss/${aux_pop}_${chr}.ss.gz > temp_files/ss.aux.${chr}
  zcat ../raw_ss/${author}/chr_ss/${author}_${chr}.ss.gz > temp_files/ss.target.${chr}


  plink --memory 4000 --threads 1 --bfile geno_files/${author}.${chr} --clump temp_files/ss.target.${chr} --clump-snp-field RSID --clump-p1 0.000001 --clump-r2 0.1 --out ${d}/out_target
  len_target=`cat ${d}/out_target.log  | fgrep "No sig" | wc -l`

  plink --memory 4000 --threads 1 --bfile geno_files/${aux_pop}.${chr} --clump temp_files/ss.aux.${chr} --clump-snp-field RSID --clump-p1 0.000001 --clump-r2 0.1 --out ${d}/out_aux
  len_aux=`cat ${d}/out_aux.log  | fgrep "No sig" | wc -l`


  ess_aux=`cat temp_files/ss.aux.${chr} | head -2 | tail -n +2 | cut -f9`
  ess_target=`cat temp_files/ss.target.${chr} | head -2 | tail -n +2 | cut -f9`

  Rscript helper_scripts/prep_xpass_ss.R $chr $author $d $ess_target target
  Rscript helper_scripts/prep_xpass_ss.R $chr $aux_pop $d $ess_aux aux

  echo STARTING R SCRIPT
  echo chr $chr
  echo author $author
  echo aux_pop $aux_pop
  echo d $d
  echo len_target $len_target
  echo len_aux $len_aux
  echo index $index

  echo;echo;echo;echo;echo;echo;echo;echo;echo;echo

  Rscript helper_scripts/run_xpass.R $chr $author $aux_pop $d $len_target $len_aux $index

  echo;echo;echo;echo;echo;echo;echo;echo;echo;echo
    

  fi
done
