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



cvmod <- cv.glmnet(y = train_df$pheno, x = as.matrix(big_data), alpha = 1, family = "binomial")
#mod <- glmnet(y = train_df$pheno, x = as.matrix(big_data), alpha = 1, lambda = cvmod$lambda.min, family = "binomial")
exit()
mod <- glmnet(y = train_df$pheno, x = as.matrix(big_data), alpha = 1, lambda = 0.0001, family = "binomial")

enet_preds <- predict(mod, as.matrix(big_data))
test_preds <- predict(mod, as.matrix(test_data))




wrong_list <- list(train_df$eid[enet_preds > quantile(enet_preds, 0.99)],
                   train_df$eid[enet_preds > quantile(enet_preds, 0.995)],
                   train_df$eid[enet_preds > quantile(enet_preds, 0.9995)],
                   train_df$eid[enet_preds > quantile(enet_preds, 0.9999)])

test_wrong_list <- list(valid_df$eid[test_preds > quantile(enet_preds, 0.99)],
                        valid_df$eid[test_preds > quantile(enet_preds, 0.995)],
                        valid_df$eid[test_preds > quantile(enet_preds, 0.9995)],
                        valid_df$eid[test_preds > quantile(enet_preds, 0.9999)])

control_wrong_list <- list(train_df$eid[enet_preds < quantile(enet_preds, 0.01)],
                           train_df$eid[enet_preds < quantile(enet_preds, 0.005)],
                           train_df$eid[enet_preds < quantile(enet_preds, 0.0005)],
                           train_df$eid[enet_preds < quantile(enet_preds, 0.0001)])

control_test_wrong_list <- list(valid_df$eid[test_preds < quantile(enet_preds, 0.01)],
                                valid_df$eid[test_preds < quantile(enet_preds, 0.005)],
                                valid_df$eid[test_preds < quantile(enet_preds, 0.0005)],
                                valid_df$eid[test_preds < quantile(enet_preds, 0.0001)])


saveRDS(mod, paste0("out_batch/wrong.model.", author, ".RDS"))

saveRDS(wrong_list, paste0("out_batch/case_groups.train.", author, ".RDS"))
saveRDS(test_wrong_list, paste0("out_batch/case_groups.test.", author, ".RDS"))

saveRDS(control_wrong_list, paste0("out_batch/control_groups.train.", author, ".RDS"))
saveRDS(control_test_wrong_list, paste0("out_batch/control_groups.test.", author, ".RDS"))
