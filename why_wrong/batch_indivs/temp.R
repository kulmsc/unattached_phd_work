x <- readRDS("temp.RDS")
train_group = x[[1]]
test_group = x[[2]]
spec_ind = x[[3]]
group_ind = x[[4]]
input_df = x[[5]]
train_df = x[[6]]
big_df = x[[7]]
#big_df <- readRDS("big_df.RDS")

library(doParallel)
library(glmnet)
registerDoParallel(4)
