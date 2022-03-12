dir_name=otherredoortho

ls peak_snps.${dir_name}.* | cut -f3 -d'.' | while read yind;do

  cat peak_snps.${dir_name}.${yind} | while read snp;do
    chromo=`zcat ../${dir_name}/full_res_Y${yind}.txt.gz | fgrep -w $snp | cut -f1 -d' '`
    peak_ind=`zcat ../${dir_name}/full_res_Y${yind}.txt.gz | fgrep -w -n $snp | cut -f1 -d':'`
    let start_ind=peak_ind-1000
    let end_ind=peak_ind+1000

    zcat ../${dir_name}/full_res_Y${yind}.txt.gz | head -${end_ind} | tail -2000 |  cut -f3 -d' ' > spec_rsids

    bgenix -g ~/athena/ukbiobank/imputed/ukbb.${chromo}.bgen -incl-rsids spec_rsids -out > test.bgen

    plink2 --bgen test.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${chromo}.sample --make-bed --out test

    plink --bfile test --r2 triangle gz --out ld_info/prod.${dir_name}.${yind}.${snp}
    mv test.bim ld_info/bim.${dir_name}.${yind}.${snp}.txt
  done

done
