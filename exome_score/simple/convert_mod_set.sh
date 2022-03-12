
ls mod_sets | while read line;do
  len=`cat mod_sets/${line} | cut -f2 | head -1 | fgrep ':' | wc -l`
  if [ $len -gt 0 ];then
    cat mod_sets/${line} | cut -f2 | cut -d':' -f1 > chromo
    cat mod_sets/${line} | cut -f2 | cut -d':' -f2 > pos
    cat mod_sets/${line} | cut -f2 | cut -d':' -f3 > a1
    cat mod_sets/${line} | cut -f2 | cut -d':' -f4 > a2
    cat mod_sets/${line} | cut -f2 > ids
    cat mod_sets/${line} | cut -f3 > beta
    paste chromo pos a1 a2 ids beta > mod_sets/${line}
    rm chromo pos a1 a2 beta ids
  fi
done
