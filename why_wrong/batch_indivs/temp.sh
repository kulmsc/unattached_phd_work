
#Rscript get_author_big_data.R $1 train
#Rscript get_author_big_data.R $1 test

Rscript get_big_data.R $1 train 
mv use_big_data*RDS data
mv use*df*RDS data
Rscript get_big_data.R $1 test
mv use_big_data*RDS data
mv use*df*RDS data
