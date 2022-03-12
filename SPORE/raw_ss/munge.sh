cat common_files/list_authors | while read author;do
  low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
  #munge_sumstats --sumstats ${author}/clean_${low_author}.txt.gz --N-col ESS --out ${author}/${low_author}.munged
  munge_sumstats --sumstats ${author}/clean_${low_author}.txt.gz --N-col ESS --merge-alleles ~/athena/refs/eur_w_ld_chr/w_hm3.snplist --out ${author}/${low_author}.for_strat

done
