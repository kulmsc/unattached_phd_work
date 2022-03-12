

ls raw_output/ | cut -f1,2 -d'.' | sort | uniq | while read data_type;do
  ls raw_output/ | cut -f4 -d'.' | sort | uniq | fgrep 0 | sort -n | while read num;do
    zcat raw_output/${data_type}.bentham.${num}.txt.gz >> combo_output/${data_type}.txt
    #zcat raw_output/${data_type}.bentham.${num}.txt.gz >> combo_output/${data_type}.txt
  done
  gzip -f combo_output/${data_type}.txt
done
