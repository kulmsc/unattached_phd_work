
cat genes_to_investigate | fgrep Casa | cut -f1 > gene_list
cat vcf_files/all_qc.cc.vcf | head -200 | grep '^#' > cand
cat vcf_files/all_qc.cc.vcf | awk '/^[^#]/{ print $0 }' | fgrep -f gene_list >> cand
mv cand candidate/gene.1.vcf


cat genes_to_investigate | fgrep Milani | cut -f1 > gene_list
cat vcf_files/all_qc.cc.vcf | head -200 | grep '^#' > cand
cat vcf_files/all_qc.cc.vcf | awk '/^[^#]/{ print $0 }' | fgrep -f gene_list >> cand
mv cand candidate/gene.2.vcf


cat genes_to_investigate | fgrep Gardner | cut -f1 > gene_list
cat vcf_files/all_qc.cc.vcf | head -200 | grep '^#' > cand
cat vcf_files/all_qc.cc.vcf | awk '/^[^#]/{ print $0 }' | fgrep -f gene_list >> cand
mv cand candidate/gene.3.vcf


cat genes_to_investigate | fgrep Custom | cut -f1 > gene_list
cat vcf_files/all_qc.cc.vcf | head -200 | grep '^#' > cand
cat vcf_files/all_qc.cc.vcf | awk '/^[^#]/{ print $0 }' | fgrep -f gene_list >> cand
mv cand candidate/gene.4.vcf

cp crispr_genes gene_list
cat vcf_files/all_qc.cc.vcf | head -200 | grep '^#' > cand
cat vcf_files/all_qc.cc.vcf | awk '/^[^#]/{ print $0 }' | fgrep -f gene_list >> cand
mv cand candidate/gene.5.vcf

