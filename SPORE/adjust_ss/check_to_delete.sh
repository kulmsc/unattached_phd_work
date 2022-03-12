method_len=`cat all_specs/method_specs | wc -l`


cat all_specs/chr_specs | while read del_chr; do
  cat all_specs/author_specs | while read del_author; do
    del_lowauthor=`echo "$del_author" | tr '[:upper:]' '[:lower:]'`

    if [ -e done_check/${del_author}.${del_chr}.done ]; then
      done_len=`cat done_check/${del_author}.${del_chr}.done | wc -l`
      if [ $done_len -eq $method_len ];then
        rm geno_files/${del_lowauthor}.${del_chr}.*
      fi
    fi

  done
done



