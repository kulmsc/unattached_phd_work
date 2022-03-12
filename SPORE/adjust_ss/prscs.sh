chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}

ess=`cat temp_files/ss.${low_author}.${chr} | head -2 | tail -n +2 | cut -f9`

cp helper_scripts/prscs_header ${d}/ss
cat temp_files/ss.${low_author}.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss

i=1
cat all_specs/prscs_param_specs | while read phi;do

  if [ ! -e ~/athena/SPORE/mod_sets/${author}/${low_author}.prs.${i}.ss ]; then

    python ~/Programs/PRScs/PRScs.py --ref_dir=/home/kulmsc/athena/refs/ldblk_1kg_eur --bim_prefix=geno_files/${low_author}.${chr} --sst_file=${d}/ss --phi=${phi} --n_gwas=${ess} --chrom=$chr --out_dir=${d}/output.$i

    Rscript helper_scripts/prscs_beta_switch.R $d $i $author $chr

  fi
  let i=i+1

done
