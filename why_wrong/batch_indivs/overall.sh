
author=$1

#Rscript get_big_data.R $author train
#Rscript get_big_data.R $author test


Rscript by_actually_has_disease.R $author
Rscript by_extreme.R $author
Rscript by_residual.R $author
Rscript by_manual_feature.R $author

#Rscript temp.R $author 
#Rscript do_performance.R $author 
#sleep 100
