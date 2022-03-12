chromo=$1

#GRCh38 ids
zcat ../delet_score/anno_output/complete.chr${chromo}.ann.vcf.gz | fgrep -v '#' | cut -f1,2,4,5 > exo_anno_1
zcat ../delet_score/anno_output/complete.chr${chromo}.ann.vcf.gz | fgrep -v '#' | cut -f8 | tr '|' '\t' | cut -f2,3,4,6,8,11 > exo_anno_2

#GRCh37
zcat ../delet_score/anno_output/impute.chr${chromo}.ann.vcf.gz | fgrep -v '#' | cut -f1,2,4,5 > impute_anno_1
zcat ../delet_score/anno_output/impute.chr${chromo}.ann.vcf.gz | fgrep -v '#' | cut -f8 | tr '|' '\t' | cut -f2,3,4,6,8,11 > impute_anno_2

#GRCh38 ids
zcat ../delet_score/anno_output/complete.chr${chromo}.dbnsfp.vcf.gz | fgrep -v '#' | fgrep dbNSFP | cut -f1,2,4,5 > exo_db_1
zcat ../delet_score/anno_output/complete.chr${chromo}.dbnsfp.vcf.gz | fgrep -v '#' | fgrep dbNSFP | cut -f8 > exo_db_2

#GRCh37
zcat ../delet_score/anno_output/impute.chr${chromo}.dbnsfp.vcf.gz | fgrep -v '#' | fgrep dbNSFP | cut -f1,2,4,5 > impute_db_1
zcat ../delet_score/anno_output/impute.chr${chromo}.dbnsfp.vcf.gz | fgrep -v '#' | fgrep dbNSFP | cut -f8 > impute_db_2

#head -1000 exo_db_2 | tr ';' '\n' | cut -f1 -d'=' | sort | uniq > poss_db
