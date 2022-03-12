library(pROC)

author <- commandArgs(trailingOnly=TRUE)[1]
#author <- "christophersen"

train_df <- readRDS("use_train_df.RDS")
test_df <- readRDS("use_test_df.RDS")

resid_groups <- readRDS(paste0("out_batch/resid_groups.train.", author, ".RDS"))
extreme_groups <- readRDS(paste0("out_batch/extreme_groups.train.", author, ".RDS"))
wrong_groups <- readRDS(paste0("out_batch/wrong_groups.train.", author, ".RDS"))
time_groups <- readRDS("out_batch/timeline_groups.train.RDS")

test_resid_groups <- readRDS(paste0("out_batch/resid_groups.test.", author, ".RDS"))
test_extreme_groups <- readRDS(paste0("out_batch/extreme_groups.test.", author, ".RDS"))
test_wrong_groups <- readRDS(paste0("out_batch/wrong_groups.test.", author, ".RDS"))
test_time_groups <- readRDS("out_batch/timeline_groups.test.RDS")


all_groups <- list(resid_groups, extreme_groups, wrong_groups, time_groups)

#------------------------------------------------------------------------
#  AUC GETTING FUNCTION
#----------------------------------------------------------------------------------

#need to write function to do cross validation
#Need to set up repeated cross validation
run_cv <- function(input_df, remove_eids, input_folds = 3, input_repeats = 3){
  eid <- input_df$eid
  all_auc <- rep(0.5, 9)
  size_group <- floor(nrow(train_df)/input_folds)
  start_spots <- floor(seq(1, size_group, length.out=input_repeats+1))
  repeat_list <- list()
  k <- 1
  for(i in 1:input_repeats){
    folds_list <- list()
    for(j in 1:input_folds){
      print(c(i,j))
      test_group <- eid[(start_spots[i]+((j-1)*size_group)):(start_spots[i]+(j*size_group))]
      train_group <- eid[!(eid %in% test_group)]
      #print("a")
      train_group <- train_group[!(train_group %in% remove_eids)]
      test_group <- test_group[!(test_group %in% remove_eids)]
      #print("b")

      mod <- glm(pheno ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + score, data = input_df[input_df$eid %in% train_group,], family = "binomial")
      #print("mod")
      preds <- predict(mod, input_df[input_df$eid %in% test_group,])
      #print("preds")
      auc <- roc(input_df[input_df$eid %in% test_group,]$pheno ~ preds)$auc
      #print("auc")

      #print(auc)
      all_auc[k] <- auc
      k <- k + 1
    }
  }
  return(all_auc)
}



#----------------------------------------------------------------
# Do IT
#-----------------------------------------------------------------------------

base_auc <- run_cv(train_df, 1)

saveRDS(base_auc, paste0("res/base_auc.", author, ".RDS"))




