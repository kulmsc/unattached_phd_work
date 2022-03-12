
rm -f list_beds.txt
for chr in {1..22}; do
  echo "/home/kulmsc/athena/ukbiobank/calls/ukbb.${chr}" >> list_beds.txt
done

plink --merge-list list_beds.txt --make-bed --out ukb_cal_allChrs

