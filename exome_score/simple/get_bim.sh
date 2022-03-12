chr=$1

cat /home/kulmsc/athena/ukbiobank/exome/ukbb.exome.${chr}.bim | awk '{print "chr"$1 "\t" $4 "\t" ($4+1) "\t" $2}' > into_bim

~/Programs/liftOver into_bim ~/athena/refs/for_picard/hg38ToHg19.over.chain.gz out_bim unMapped
