
plink2 --pfile all_qc --make-bed --out all_qc --max-alleles 2

plink --bfile all_qc --pheno pheno_for_gwas --assoc fisher-midp --out out

plink2 --pfile all_qc --pheno pheno_for_gwas --glm sex hide-covar --covar covar_for_gwas --covar-variance-standardize --out out

cat out.assoc.fisher | awk '$8 < 0.000001 {print $1 " " $3}' > ../candidate/plink.1.pos

cat out.assoc.logistic | awk '$13 < 0.005 {print $1 " " $2}' > ../candidate/plink.2.pos

cat all_qc.cc.vcf | head -200 | grep '^#' > cand
cat ../candidate/plink.1.pos | while read line;do
  takechr=`echo $line | cut -f1 -d' '`
  takepos=`echo $line | cut -f2 -d' '`
  cat all_qc.cc.vcf | awk '/^[^#]/{ print $0 }' | awk -v var="$takepos" '$2 == var {print $0}' | awk -v var="$takechr" '$1 == var {print $0}' >> cand
done
mv cand ../candidate/plink.1.vcf


cat all_qc.ann.vcf | head -200 | grep '^#' > cand
cat ../candidate/plink.2.pos | while read line;do
  takechr=`echo $line | cut -f1 -d' '`
  takepos=`echo $line | cut -f2 -d' '`
  cat all_qc.cc.vcf | awk '/^[^#]/{ print $0 }' | awk -v var="$takepos" '$2 == var {print $0}' | awk -v var="$takechr" '$1 == var {print $0}' >> cand
done
mv cand ../candidate/plink.2.vcf

