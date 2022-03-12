
chromo=`cat to_investigate/chromo`
cat ~/athena/ukbiobank/qc/imputed/ukb_mfi_chr${chromo}_v3.txt | fgrep -w -f to_investigate/test_pos > to_investigate/small_qc
