
author=$1
chr=$2
train_frac=0.2

lowauthor=`echo "$author" | tr '[:upper:]' '[:lower:]'`

if [ ! -e geno_files/${lowauthor}.${chr}.bed ];then
  echo no geno files within setup.sh
  zcat ../raw_ss/${author}/chr_ss/${lowauthor}_${chr}.ss.gz | cut -f3 > temp_files/rsids.${author}.${chr}
  zcat ../raw_ss/${author}/chr_ss/${lowauthor}_${chr}.ss.gz > temp_files/ss.${lowauthor}.${chr}

  bgenix -g ~/athena/ukbiobank/imputed/ukbb.${chr}.bgen -incl-rsids temp_files/rsids.${author}.${chr} > temp_files/${author}.${chr}.bgen

  plink2_new --memory 8000 --bgen temp_files/${author}.${chr}.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${chr}.sample --keep-fam ~/athena/doc_score/qc/cv_files/train_eid.${train_frac}.txt --maf 0.01 --geno 0.1 --hwe 1e-50 midp --make-bed --out geno_files/${lowauthor}.${chr}

  rm temp_files/${author}.${chr}.bgen
fi
