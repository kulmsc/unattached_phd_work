spack load /g7nwnaf

dir_name=otherredoortho
for chrom in {1..22};do
  for y_index in {1..8};do

    echo "SNP P N" > current_pval_file
    zcat ../${dir_name}/full_res_Y${y_index}.txt.gz | tail -n+2 | awk -v var="$chrom" '$1 == var {print $0}' | awk '{print $3 " " 10^-($13) " " $8}' >> current_pval_file
    cat ~/athena/ukbiobank/calls/ukbb.${chrom}.fam | cut -f1 -d' ' | head -25000 > list_fam
    plink --bfile ~/athena/ukbiobank/calls/ukbb.${chrom} --keep-fam list_fam --make-bed --out current

    zcat ../ortho/full_res_Y1.txt.gz | tail -n+2 | awk -v var="$chrom" '$1 == var {print $0}' | awk '{print $3 " " $1 " " $2}' > curr_snploc

    magma --annotate window=5,1.5 --snp-loc curr_snploc --gene-loc NCBI37.3.gene.loc --out test
    #SNP location file should contain three columns: SNP ID, chromosome, and base pair position
    #The gene location file must contain at least four columns, in this order: gene ID, chromosome,start site, stop site

    magma --bfile current synonyms=dbsnp151.synonyms --gene-annot test.genes.annot --pval current_pval_file ncol=N
    #The p-value file must be a plain text data file with each row corresponding to a SNP. If MAGMA
    #detects a header in the file it will look for SNP IDs and p-values in the SNP and P column respectively.
    #If no header is found it will use the first column for SNP IDs and the second column for p-values

    rm test*
    rm curr_snploc
    rm current*

    Rscript finish.R $dir_name $y_index $chrom
  done
done
