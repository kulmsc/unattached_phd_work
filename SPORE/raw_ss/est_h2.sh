#cat common_files/list_authors | while read author;do
author=european
ethnic=eur

  ldsc --h2 ${author}/${author}.munged.sumstats.gz --ref-ld-chr ~/athena/refs/custom_ldscores/$ethnic --w-ld-chr ~/athena/refs/custom_ldscores/$ethnic --out ${author}/${author}.h2

#done
