
eth=eur

for i in {1..22};do

#zcat baseline_scores/${eth}.${i}.l2.ldscore.gz | cut -f2 > base_rsid
#zcat impact_files/curr_impact.${eth}.${i}.l2.ldscore.gz | cut -f2 > impact_rsid

#zcat baseline_scores/${eth}.${i}.l2.ldscore.gz | head -1 > ready_files/${eth}.${i}.l2.ldscore
#zcat impact_files/curr_impact.${eth}.${i}.l2.ldscore.gz | head -1 > ready_files/curr_impact.${eth}.${i}.l2.ldscore

#zcat baseline_scores/${eth}.${i}.l2.ldscore.gz | fgrep -w -f impact_rsid > ready_files/${eth}.${i}.l2.ldscore
#zcat impact_files/curr_impact.${eth}.${i}.l2.ldscore.gz | fgrep -w -f base_rsid > ready_files/curr_impact.${eth}.${i}.l2.ldscore

#gzip ready_files/${eth}.${i}.l2.ldscore
#gzip ready_files/curr_impact.${eth}.${i}.l2.ldscore

#must cd first
#ln -s ../impact_files/curr_impact.${eth}.${i}.l2.M_5_50 curr_impact.${eth}.${i}.l2.M_5_50
#ln -s ../impact_files/curr_impact.${eth}.${i}.l2.M curr_impact.${eth}.${i}.l2.M
#ln -s ../baseline_scores/${eth}.${i}.l2.M_5_50 ${eth}.${i}.l2.M_5_50
#ln -s ../baseline_scores/${eth}.${i}.l2.M ${eth}.${i}.l2.M
#ln -s ../baseline_annot/baseline.${i}.annot.gz ${eth}.${i}.annot.gz


done
