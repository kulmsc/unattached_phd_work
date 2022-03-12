#ls res/*chr1_P1.reg* | while read line;do
#  num=`echo $line | cut -f3 -d"_"`
#  zcat $line > full_res_Y.${num}.txt
#done

cat res/ukb_step2_9_full_covar_chr6_P1.regenie > full_res_Y.9.txt

for chr in {7..22};do
  #ls res/*chr${chr}_P1.reg* | while read line;do
  ls res | fgrep regenie | fgrep _9_ | fgrep chr${chr}_P1 | while read line;do
    num=`echo $line | cut -f3 -d"_"`
    zcat res/$line | tail -n+2 >> full_res_Y.${num}.txt
  done
done
