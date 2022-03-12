
author=$1
chr=$2
train_frac=0.2

lowauthor=`echo "$author" | tr '[:upper:]' '[:lower:]'`
realchr=`echo $chr | cut -f1 -d'_'`

if [ ! -e geno_files/${lowauthor}.${chr}.bed ];then
  echo no geno files within setup.sh
  zcat ../finished_ss/${author}/${lowauthor}.${chr}.ss.gz | head
  zcat ../finished_ss/${author}/${lowauthor}.${chr}.ss.gz | cut -f3 > temp_files/rsids.${author}.${chr}
  echo "CHR BP RSID A1 A2 SE BETA P" | tr ' ' '\t' > temp_files/ss.${lowauthor}.${chr}
  zcat ../finished_ss/${author}/${lowauthor}.${chr}.ss.gz >> temp_files/ss.${lowauthor}.${chr}

  bgenix -g ~/athena/ukbiobank/imputed/ukbb.${realchr}.bgen -incl-rsids temp_files/rsids.${author}.${chr} > temp_files/${author}.${chr}.bgen

  plink2_new --bgen temp_files/${author}.${chr}.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${realchr}.sample --keep-fam ~/athena/doc_score/qc/cv_files/train_eid.${train_frac}.txt --maf 0.01 --geno 0.1 --hwe 1e-50 midp --make-bed --out temp_files/almost.${lowauthor}.${chr}

  cat temp_files/almost.${lowauthor}.${chr}.bim | cut -f2 | sort | uniq -d > temp_files/dup_ids.${lowauthor}.${chr}

  plink --bfile temp_files/almost.${lowauthor}.${chr} --exclude temp_files/dup_ids.${lowauthor}.${chr} --make-bed --out geno_files/${lowauthor}.${chr}

  rm temp_files/dup_ids.${lowauthor}.${chr}
  rm temp_files/almost.${lowauthor}.${chr}*
  rm temp_files/${author}.${chr}.bgen
fi
