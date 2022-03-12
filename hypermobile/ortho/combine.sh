num=diabeteshypocovar
ls *ukb_step2_wdiabeteshypo_chr1_* | while read line;do
  cat $line > full_res_Y${num}.txt
  let num=num+1
done

for chr in {2..22};do
#num=1
  ls *ukb_step2_wdiabeteshypo_chr${chr}_* | while read line;do
    cat $line | tail -n+2 >> full_res_Y${num}.txt
    #let num=num+1
  done
done
