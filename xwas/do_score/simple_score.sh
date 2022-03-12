
ss_name=$1
author=$2
method=$3
chr=$4
i=$5

echo $ss_name $author $method $chr $i

#Just logistics around controlling parallelism
go_num=`head -1 temp_files/poss_go`
grep -v -w $go_num temp_files/poss_go > temp_files/temp_poss
mv temp_files/temp_poss temp_files/poss_go

low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
ver=`echo $ss_name | cut -f4 -d'.'`
use_chr=`echo $chr | cut -f1 -d'_'`

if [ ! -e small_score_files/score.${low_author}.${chr}.${ver}.${method}.profile.zst ];then

  cat ../mod_sets/${author}/${ss_name} | cut -f3  > temp_files/rsids.${i}

  bgenix -g ~/athena/ukbiobank/imputed/ukbb.${use_chr}.bgen -incl-rsids temp_files/rsids.${i} > temp_files/temp.${i}.bgen

  plink2_new --memory 12000 --threads 12 --bgen temp_files/temp.${i}.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${use_chr}.sample --keep-fam temp_files/brit_eid --make-bed --out temp_files/geno.${i}

  rm temp_files/temp.${i}.bgen

  cat temp_files/geno.${i}.bim | cut -f2 | sort | uniq -d > temp_files/dupids.${i}

  plink --bfile temp_files/geno.${i} --exclude temp_files/dupids.${i} --make-bed --out temp_files/geno2.${i}

  plink --memory 12000 --threads 12 --bfile temp_files/geno2.${i} --keep-allele-order --score ../mod_sets/${author}/${ss_name} 3 4 7 sum --out small_score_files/score.${low_author}.${chr}.${ver}.${method}

  zstd --rm small_score_files/score.${low_author}.${chr}.${ver}.${method}.profile

  rm temp_files/rsids.${i}
  rm temp_files/geno.${i}.*
  rm temp_files/geno2.${i}.*
  rm temp_files/dupids.${i}

else

  echo already exists

fi

#Just logistics around controlling parallelism
echo $go_num >> temp_files/poss_go

