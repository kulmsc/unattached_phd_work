
ls many_mod_sets | while read line;do
  len=`cat many_mod_sets/${line} | cut -f2 | head -1 | fgrep ':' | wc -l`
  if [ $len -gt 0 ];then
    cat many_mod_sets/${line} | cut -f2 | cut -d':' -f1 > chromo
    cat many_mod_sets/${line} | cut -f2 | cut -d':' -f2 > pos
    cat many_mod_sets/${line} | cut -f2 | cut -d':' -f3 > a1
    cat many_mod_sets/${line} | cut -f2 | cut -d':' -f4 > a2
    cat many_mod_sets/${line} | cut -f2 > ids
    cat many_mod_sets/${line} | cut -f3 > beta
    paste chromo pos a1 a2 ids beta > many_mod_sets/${line}
    rm chromo pos a1 a2 beta ids
  fi
done
