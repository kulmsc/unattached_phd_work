
ss_name=$1
prefix=$2
chr=$3
i=$4

echo $ss_name $prefix $chr $i

#Just logistics around controlling parallelism
#go_num=`head -1 temp_files/poss_go`
#grep -v -w $go_num temp_files/poss_go > temp_files/temp_poss
#mv temp_files/temp_poss temp_files/poss_go


if [ ! -e full_small_score_files/score.${prefix}.${ss_name}.${chr}.profile.zst ];then


  #The VCF ###
  plink --memory 12000 --threads 12 --bfile use_done --keep-allele-order --score ../mod_sets/${ss_name}.${chr}.ss 3 4 7 sum --allow-extra-chr --out small_score_files/score.${prefix}.full.${ss_name}.${chr}

  zstd --rm small_score_files/score.${prefix}.full.${ss_name}.${chr}.profile


#    ls ../impute/output_files/GEN*bed | cut -f4 -d'/' | cut -f1 -d'.' | sort | uniq | while read newfix;do
#    for subset in variable ref;do
#      cat ../impute/output_files/${newfix}.${chr}.${subset}.bim | cut -f2 > exrsid.${newfix}.${subset}
#      plink --memory 12000 --threads 12 --bfile use_done --extract exrsid.${newfix}.${subset} --keep-allele-order --score ../mod_sets/${ss_name}.${chr}.ss 3 4 7 sum --allow-extra-chr --out full_small_score_files/score.${newfix}.${subset}.${ss_name}.${chr}
#      zstd --rm full_small_score_files/score.${newfix}.${subset}.${ss_name}.${chr}.profile
#    done
#    done

else

  echo already exists

fi

