
#y_index=2
dir_name=otherredoortho
for y_index in {1..8};do

Rscript custom_clump.R 1 $y_index $dir_name

if [ -e peak_snps.${dir_name}.${y_index} ];then
  cat peak_snps.${dir_name}.${y_index} | while read snp;do
    chromo=`zcat ../${dir_name}/full_res_Y${y_index}.txt.gz | fgrep -w $snp | cut -f1 -d' '`
    peak_ind=`zcat ../${dir_name}/full_res_Y${y_index}.txt.gz | fgrep -w -n $snp | cut -f1 -d':'`
    let start_ind=peak_ind-25000
    let end_ind=peak_ind+25000

    zcat ../${dir_name}/full_res_Y${y_index}.txt.gz | head -${end_ind} | tail -50000 |  awk '{print $3,$4,$5,$6,$10,$11,10**-$13,$8}' > current.ma

    cat current.ma | cut -f1 -d' ' > temp
    bgenix -g ~/athena/ukbiobank/imputed/ukbb.${chromo}.bgen -incl-rsids temp -out > test.bgen

    plink2 --bgen test.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${chromo}.sample --make-bed --out test

    cat test.bim  | cut -f2 | sort | uniq -d > bad_snps

    plink --bfile test --exclude bad_snps --make-bed --out nowready
    cat nowready.bim | cut -f2 > plink_snps

    cat current.ma | fgrep -w -v -f bad_snps | fgrep -w -f plink_snps > nowready.ma

    gcta64  --bfile nowready --maf 0.01 --cojo-file nowready.ma --cojo-p 5e-7 --cojo-slct --out results/COJO.${y_index}.${snp}.${dir_name}

    rm plink_snps nowready* bad_snps test* temp current.ma
  done
fi

done
