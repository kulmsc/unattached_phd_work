library(glmnet)
library(Hmisc)

#author <- commandArgs(trailingOnly=TRUE)[1]
author <- "bentham"

#system(paste0("Rscript get_big_data.R ", author))
big_data <- readRDS(paste0("data/use_big_data.train.", author, ".RDS"))
train_df <- readRDS(paste0("data/use_train_df.", author, ".RDS"))
test_data <- readRDS(paste0("data/use_big_data.test.", author, ".RDS"))
valid_df <- readRDS(paste0("data/use_test_df.", author, ".RDS"))

big_data <- big_data[,colnames(big_data) %in% colnames(test_data)]
test_data <- test_data[,colnames(test_data) %in% colnames(big_data)]

#------------------------------------------------------
# Enet model over all features
#-----------------------------------------------------------

mod_func <- function(input_pheno, input_df, comp_by){

  cvmod <- cv.glmnet(y = input_pheno, x = as.matrix(input_df), alpha = 1)
  mod <- glmnet(y = input_pheno, x = as.matrix(input_df), alpha = 1, lambda = cvmod$lambda.min)


  enet_preds <- predict(mod, as.matrix(big_data))
  test_preds <- predict(mod, as.matrix(test_data))


  case_wrong_list <- list(train_df$eid[enet_preds > quantile(enet_preds, 0.99)],
                          train_df$eid[enet_preds > quantile(enet_preds, 0.995)],
                          train_df$eid[enet_preds > quantile(enet_preds, 0.9995)])

  case_test_wrong_list <- list(valid_df$eid[test_preds > quantile(enet_preds, 0.99)],
                               valid_df$eid[test_preds > quantile(enet_preds, 0.995)],
                               valid_df$eid[test_preds > quantile(enet_preds, 0.9995)])

  control_wrong_list <- list(train_df$eid[enet_preds < quantile(enet_preds, 0.01)],
                             train_df$eid[enet_preds < quantile(enet_preds, 0.005)],
                             train_df$eid[enet_preds < quantile(enet_preds, 0.0005)])

  control_test_wrong_list <- list(valid_df$eid[test_preds < quantile(enet_preds, 0.01)],
                                  valid_df$eid[test_preds < quantile(enet_preds, 0.005)],
                                  valid_df$eid[test_preds < quantile(enet_preds, 0.0005)])


  return(list(case_wrong_list, case_test_wrong_list, control_wrong_list, control_test_wrong_list))

  #saveRDS(mod, paste0("out_batch/extreme_", comp_by, ".model.", author, ".RDS"))
  #saveRDS(wrong_list, paste0("out_batch/extreme_groups_", comp_by, ".train.", author, ".RDS"))
  #saveRDS(test_wrong_list, paste0("out_batch/extreme_groups_", comp_by, ".test.", author, ".RDS"))

}

case_train <- train_df[train_df$pheno == 0,]
control_train <- train_df[train_df$pheno == 1,]

top_case <- case_train$eid[case_train$score > sort(case_train$score, decreasing = T)[100]]
top_control <- control_train$eid[control_train$score > sort(control_train$score, decreasing = T)[100]]

bottom_case <- case_train$eid[case_train$score < sort(case_train$score, decreasing = F)[100]]
bottom_control <- control_train$eid[control_train$score < sort(control_train$score, decreasing = F)[100]]

first = mod_func(train_df$pheno[train_df$eid %in% c(top_case, top_control)], big_data[train_df$eid %in% c(top_case, top_control),], "disease_in_hi_prs")
second = mod_func(train_df$pheno[train_df$eid %in% c(bottom_case, bottom_control)], big_data[train_df$eid %in% c(bottom_case, bottom_control),], "disease_in_lo_prs")
third = mod_func(train_df$pheno[train_df$eid %in% c(top_control, bottom_case)], big_data[train_df$eid %in% c(top_control, bottom_case),], "hi_lo_all")

#first second third - all list of lists, first set of elements are for case, control, train, test, - second set of elements are for cut offs

case_wrong_list <- c(first[[1]], second[[1]], third[[1]])
case_test_wrong_list <- c(first[[2]], second[[2]], third[[2]])

control_wrong_list <- c(first[[3]], second[[3]], third[[3]])
control_test_wrong_list <- c(first[[4]], second[[4]], third[[4]])

#each of these are 9 elements long 

saveRDS(case_wrong_list, paste0("out_batch/case_extreme_groups.train.", author, ".RDS"))
saveRDS(case_test_wrong_list, paste0("out_batch/case_extreme_groups.test.", author, ".RDS"))

saveRDS(control_wrong_list, paste0("out_batch/control_extreme_groups.train.", author, ".RDS"))
saveRDS(control_test_wrong_list, paste0("out_batch/control_extreme_groups.test.", author, ".RDS"))
