
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


if [ ! -e small_score_files/score.${low_author}.${chr}.${ver}.${method}.profile.zst ];then

  cat ../mod_sets/${author}/${ss_name} | cut -f3 > temp_files/rsids.${i}

  cat temp_files/rsids.${i} | sort | uniq -d > temp_files/dupids.${i}

  cat ../mod_sets/${author}/${ss_name} | fgrep -w -v -f temp_files/dupids.${i} > temp_files/local_ss.${i}

  plink --memory 12000 --threads 12 --bfile temp_files/geno.${low_author}.${chr}  --exclude temp_files/dupids.${i} --keep-allele-order --score temp_files/local_ss.${i} 3 4 7 sum --out small_score_files/score.${low_author}.${chr}.${ver}.${method}

  zstd --rm small_score_files/score.${low_author}.${chr}.${ver}.${method}.profile


else

  echo already exists

fi

echo $go_num >> temp_files/poss_go
