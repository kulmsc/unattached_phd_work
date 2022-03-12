chr=$1
author=$2
dir=$3
d=comp_zone/dir${dir}
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`

cat temp_files/ss.${low_author}.${chr} |  awk '$6 > 0 {print $0}' > ${d}/newss
#cat ${d}/newss | cut -f3 | tail -n+2 | fgrep -w -f ~/athena/refs/hapmap_from_ldpred2 > ${d}/newrsids
cat ${d}/newss | cut -f3 | tail -n+2  > ${d}/newrsids

head -10000 geno_files/${low_author}.${chr}.fam | cut -f1 > ${d}/subset_inds
plink --bfile geno_files/${low_author}.${chr} --keep-fam ${d}/subset_inds --extract ${d}/newrsids --make-bed --out ${d}/for_ldpred2

Rscript helper_scripts/ldpred2.R $chr $author $d
