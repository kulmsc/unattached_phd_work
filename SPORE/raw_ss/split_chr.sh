cat common_files/list_authors | while read author;do
  low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`

  mkdir ${author}/chr_ss
  for chr in {1..23};do
    zcat ${author}/clean_${low_author}.txt.gz | head -1 > ${author}/chr_ss/${low_author}_${chr}.ss
    zcat ${author}/clean_${low_author}.txt.gz | awk -v var="$chr" '$1 == var {print $0}' >> ${author}/chr_ss/${low_author}_${chr}.ss
    gzip ${author}/chr_ss/${low_author}_${chr}.ss
  done

done
