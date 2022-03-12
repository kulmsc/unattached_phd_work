chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}


if [ $author == "total" ];then


zcat ../raw_ss/european/chr_ss/european_${chr}.ss.gz > temp_files/ss.european.${chr}
zcat ../raw_ss/hispanic/chr_ss/hispanic_${chr}.ss.gz > temp_files/ss.hispanic.${chr}
zcat ../raw_ss/african/chr_ss/african_${chr}.ss.gz > temp_files/ss.african.${chr}
zcat ../raw_ss/eastasian/chr_ss/eastasian_${chr}.ss.gz > temp_files/ss.eastasian.${chr}
zcat ../raw_ss/total/chr_ss/total_${chr}.ss.gz > temp_files/ss.total.${chr}

ess_eur=`cat temp_files/ss.european.${chr} | head -2 | tail -n +2 | cut -f9`
ess_afr=`cat temp_files/ss.african.${chr} | head -2 | tail -n +2 | cut -f9`
ess_hsp=`cat temp_files/ss.hispanic.${chr} | head -2 | tail -n +2 | cut -f9`
ess_eas=`cat temp_files/ss.eastasian.${chr} | head -2 | tail -n +2 | cut -f9`
ess_tot=`cat temp_files/ss.total.${chr} | head -2 | tail -n +2 | cut -f9`


cp helper_scripts/prscs_header ${d}/ss.eur
cat temp_files/ss.european.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss.eur

cp helper_scripts/prscs_header ${d}/ss.afr
cat temp_files/ss.african.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss.afr

cp helper_scripts/prscs_header ${d}/ss.eas
cat temp_files/ss.eastasian.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss.eas

cp helper_scripts/prscs_header ${d}/ss.hsp
cat temp_files/ss.hispanic.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss.hsp

cp helper_scripts/prscs_header ${d}/ss.tot
cat temp_files/ss.total.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss.tot



i=1
cat all_specs/prscs_param_specs | while read phi;do
#phi=0.0001

#if [ ! -e ~/athena/doc_score/mod_sets/total/total.${chr}.prscsx.${i}.ss ]; then

#    python ~/Programs/PRScsx/PRScsx.py --ref_dir=/home/kulmsc/athena/refs --bim_prefix=geno_files/${low_author}.${chr} --sst_file=${d}/ss.afr,${d}/ss.eas,${d}/ss.hsp,${d}/ss.tot --phi=${phi} --n_gwas=${ess_afr},${ess_eas},${ess_hsp},${ess_tot} --pop=afr,eas,amr,eur --chrom=$chr --out_dir=${d} --out_name=x_afr,x_eas,x_hsp,x_tot

#    Rscript helper_scripts/prscsx_beta_switch.R $d afr african $chr $i
#    Rscript helper_scripts/prscsx_beta_switch.R $d eas eastasian $chr $i
#    Rscript helper_scripts/prscsx_beta_switch.R $d hsp hispanic $chr $i
#    Rscript helper_scripts/prscsx_beta_switch.R $d tot total $chr $i
#fi

#    cp ${d}/x* temp_files
#    rm ${d}/x*


    let i=i+1

#if [ ! -e ~/athena/doc_score/mod_sets/european/european.${chr}.prscsx.${i}.ss ]; then

    python ~/Programs/PRScsx/PRScsx.py --ref_dir=/home/kulmsc/athena/refs --bim_prefix=geno_files/${low_author}.${chr} --sst_file=${d}/ss.eur,${d}/ss.afr,${d}/ss.eas,${d}/ss.hsp --phi=${phi} --n_gwas=${ess_eur},${ess_afr},${ess_eas},${ess_hsp} --pop=eur,afr,eas,amr --chrom=$chr --out_dir=${d} --out_name=x_eur,x_afr,x_eas,x_hsp


    Rscript helper_scripts/prscsx_beta_switch.R $d afr african $chr $i $phi
    Rscript helper_scripts/prscsx_beta_switch.R $d eas eastasian $chr $i $phi
    Rscript helper_scripts/prscsx_beta_switch.R $d hsp hispanic $chr $i $phi
    Rscript helper_scripts/prscsx_beta_switch.R $d eur european $chr $i $phi

#fi

    let i=i+1

    #cp ${d}/x* temp_files
    #rm ${d}/x*

done


fi
