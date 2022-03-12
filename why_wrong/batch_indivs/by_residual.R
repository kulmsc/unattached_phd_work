library(glmnet)
library(Hmisc)

#author <- commandArgs(trailingOnly=TRUE)[1]
author <- "rheenen"

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



#now do I pick based on the coef of the best model, or do I force a select number of nonzero coefs?
#there is an argument just to keep everything and move onto the cross validation (but since there are infinite ways to remove indivs from a continous feature this does not seem possible)
if("sex" %in% colnames(train_df)){
  base_mod <- glm(pheno ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, train_df, family = "binomial")
} else {
  base_mod <- glm(pheno ~ age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, train_df, family = "binomial")
}
resids <- base_mod$residuals

cvmod <- cv.glmnet(y = resids, x = as.matrix(big_data), alpha = 1)
mod <- glmnet(y = resids, x = as.matrix(big_data), alpha = 1, lambda = cvmod$lambda.min)



enet_preds <- predict(mod, as.matrix(big_data))
test_preds <- predict(mod, as.matrix(test_data))


wrong_list <- list(train_df$eid[enet_preds > quantile(enet_preds, 0.99)],
                   train_df$eid[enet_preds > quantile(enet_preds, 0.995)],
                   train_df$eid[enet_preds > quantile(enet_preds, 0.9995)],
                   train_df$eid[enet_preds > quantile(enet_preds, 0.9999)])

test_wrong_list <- list(valid_df$eid[test_preds > quantile(test_preds, 0.99)],
                        valid_df$eid[test_preds > quantile(test_preds, 0.995)],
                        valid_df$eid[test_preds > quantile(test_preds, 0.9995)],
                        valid_df$eid[test_preds > quantile(test_preds, 0.9999)])

saveRDS(mod, paste0("out_batch/resid.model.", author, ".RDS"))
saveRDS(wrong_list, paste0("out_batch/resid_groups.train.", author, ".RDS"))
saveRDS(test_wrong_list, paste0("out_batch/resid_groups.test.", author, ".RDS"))
