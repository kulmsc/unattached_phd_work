
cat vcf_files/all_qc.ann.vcf | head -200 | grep '^#' > cand

cat genes_to_investigate | cut -f1 > temp
cat vcf_files/all_qc.cc.vcf | fgrep -w -f temp >> cand

plink --vcf cand --make-bed --out for_inter
plink --vcf cand --export A --out for_inter
