rm big_small_scores/*

#watch out for the grep !!!
for chr in {1..22};do
#for chr in 22;do

  #for cauthor in Daner;do
  cat all_specs/author_specs | while read cauthor;do
  #cauthor=Zeginni

    low_author=`echo "$cauthor" | tr '[:upper:]' '[:lower:]'`

    Rscript align_sumstats.R $chr $cauthor

    cat big_mod_set | cut -f1 | tail -n+2 > temp_files/all_rsid

    num_cols=`head -1 big_mod_set | cut -f3-300 | tr '\t' '\n' | wc -l`

    bgenix -g ~/athena/ukbiobank/imputed/ukbb.${chr}.bgen -incl-rsids temp_files/all_rsid > temp_files/temp.bgen

    plink2_new --memory 12000 --threads 12 --bgen temp_files/temp.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${chr}.sample --score big_mod_set 1 2 header-read cols=+scoresums --score-col-nums 3-${num_cols} --out big_small_scores/res.${chr}.${low_author}

    rm temp_files/temp.bgen temp_files/all_rsid big_mod_set

    gzip big_small_scores/res.${chr}.${low_author}.sscore

  done
done

