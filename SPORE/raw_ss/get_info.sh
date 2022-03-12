

echo author,trait,sampe_size,cases,controls,ancestry,snps,ldsc_h2,ldsc_h2_se,hdl_h2,hdl_h2se > meta_stats

cat common_files/list_authors | while read author;do
  low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`

  trait=`cat ${author}/notes | head -2 | tail -1`
  samp_size=`cat ${author}/notes | head -3 | tail -1`
  cases=`cat ${author}/notes | head -4 | tail -1`
  controls=`cat ${author}/notes | head -5 | tail -1`
  ancestry=`cat ${author}/notes | head -6 | tail -1`
  snps=`zcat ${author}/clean_${low_author}.txt.gz | wc -l`
  ldsc_h2=`cat ${author}/${low_author}.h2.log | fgrep "scale h2" | cut -f2 -d':' | cut -f2 -d' '`
  ldsc_h2se=`cat ${author}/${low_author}.h2.log | fgrep "scale h2" | cut -f2 -d':' | cut -f3 -d' ' | cut -f2 -d'(' | cut -f1 -d')'`
  hdl_h2=`cat ${author}/${low_author}.HDL.h2.Rout | fgrep Heritability | cut -f6 -d' '`
  hdl_h2se=`cat ${author}/${low_author}.HDL.h2.Rout | fgrep Heritability | cut -f7 -d' ' | cut -f2 -d'(' | cut -f1 -d')'`

  echo $author,$trait,$samp_size,$cases,$controls,$ancestry,$snps,$ldsc_h2,$ldsc_h2se,$hdl_h2,$hdl_h2se >> meta_stats
done
