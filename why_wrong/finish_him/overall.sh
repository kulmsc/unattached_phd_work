author=$1
#author=bentham

#Rscript get_big_data.R $author train
#Rscript get_big_data.R $author test

for stat in all_around mean min count;do
  Rscript full_test.R $author $stat super_new
done
