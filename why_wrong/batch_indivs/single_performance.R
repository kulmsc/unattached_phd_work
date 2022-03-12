library(pROC)

#author <- commandArgs(trailingOnly=TRUE)[1]
group_by <- "resid"
author <- "christophersen"

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
run_cv <- function(input_df, remove_eids, input_folds = 3, input_repeats = 30){

  #THIS IS A NEW LINE THAT MIGHT HELP !!!!!!!!!!!!!
  #remove_eids <- intersect(remove_eids, input_df$df[input_df$pheno ==0])

  eid <- input_df$eid
  all_auc <- rep(0.5, input_folds * input_repeats)
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

base_auc <- mean(run_cv(train_df, 1))

step_through <- function(all_groups, start_remove, group_by){

  auc_min <- list(rep(NA, length(all_groups[[1]])),
                   rep(NA, length(all_groups[[2]])),
                   rep(NA, length(all_groups[[3]])),
                   rep(NA, length(all_groups[[4]])))
  auc_mean <- list(rep(NA, length(all_groups[[1]])),
                   rep(NA, length(all_groups[[2]])),
                   rep(NA, length(all_groups[[3]])),
                   rep(NA, length(all_groups[[4]])))

  for(i in 1:length(all_groups)){
    for(j in 1:length(all_groups[[i]])){
      cvres <- run_cv(train_df, c(all_groups[[i]][[j]], start_remove ))
      auc_mean[[i]][j] <- mean(cvres) - base_auc
      auc_min[[i]][j] <- min(cvres) - base_auc
    }
  }

  auc_perf <- auc_mean

  if(group_by == "resid"){
    perm_remove <- c(all_groups[[1]][which.max(auc_perf[[1]])], start_remove)
    best_auc <- max(auc_perf[[1]])
    all_groups[[1]][which.max(auc_perf[[1]])] <- 123
    remove_group <- 1
    remove_index <- which.max(auc_perf[[1]] )

  } else if(group_by == "extreme") {
    perm_remove <- c(all_groups[[2]][which.max(auc_perf[[2]] )], start_remove)
    best_auc <- max(auc_perf[[2]])
    all_groups[[2]][which.max(auc_perf[[2]] )] <- 123
    remove_group <- 2
    remove_index <- which.max(auc_perf[[2]] )

  } else if(group_by == "wrong") {
    perm_remove <- c(all_groups[[3]][which(auc_perf[[3]] == max(unlist(auc_perf)))], start_remove)
    best_auc <- max(auc_perf[[3]])
    all_groups[[3]][which(auc_perf[[3]] == max(unlist(auc_perf)))] <- 123
    remove_group <- 3
    remove_index <- which(auc_perf[[3]] == max(unlist(auc_perf)))

  } else if(group_by == "time") {
    perm_remove <- c(all_groups[[4]][which.max(auc_perf[[4]] )], start_remove)
    best_auc <- max(auc_perf[[4]])
    all_groups[[4]][which.max(auc_perf[[4]] )] <- 123
    remove_group <- 4
    remove_index <- which.max(auc_perf[[4]] )


  } else {
    print("butts")
  }

  return(list(best_auc , all_groups, perm_remove, remove_group, remove_index))
}

#need each of the all_groups to have names

print(1)
back_list <- step_through(all_groups, c(123456789), "resid")
best_auc <- back_list[[1]]
remove_group <- back_list[[3]]
remove_index <- back_list[[5]]

if(group_by == "resid"){
  best_test_remove <- test_resid_groups[[remove_index]]
} else if(group_by == "exteme"){
  best_test_remove <- test_extreme_groups[[remove_index]]
} else if(group_by == "wrong"){
  best_test_remove <- test_wrong_groups[[remove_index]]
} else if(group_by == "time"){
  best_test_remove <- test_time_groups[[remove_index]]
}

saveRDS(list(remove_group, best_test_remove, best_auc, NA, NA), paste0("res/perf.", author, ".", group_by, ".RDS"))



