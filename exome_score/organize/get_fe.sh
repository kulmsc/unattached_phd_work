
cat ~/athena/ukbiobank/setup_exome/ukb_fe_exm_chrall_v1.bim | awk '{print "chr"$1 "\t" $4 "\t" ($4+1) "\t" $2}' > into_fe

~/Programs/liftOver into_fe ~/athena/refs/for_picard/hg38ToHg19.over.chain.gz out_fe unMapped
