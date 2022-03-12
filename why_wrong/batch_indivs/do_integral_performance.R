library(pROC)
library(jaccard)
library(glmnet)

author <- commandArgs(trailingOnly=TRUE)[1]
#author <- "bentham"

source("integral_funcs.R")


big_data <- readRDS(paste0("data/use_big_data.train.", author, ".RDS"))
train_df <- readRDS(paste0("data/use_train_df.", author, ".RDS"))
test_data <- readRDS(paste0("data/use_big_data.test.", author, ".RDS"))
test_df <- readRDS(paste0("data/use_test_df.", author, ".RDS"))

uval <- apply(big_data, 2, function(x) length(unique(x)))
limited_bins <- colnames(big_data)[uval==2][colSums(big_data[,uval==2]) < 750]
#zval <- apply(big_data[,uval>2], 2, function(x) sum(x==0))
#many_zero <- colnames(big_data)[uval>2][zval > 10000]
big_data <- big_data[,!(colnames(big_data) %in% limited_bins)]

uval <- apply(big_data, 2, function(x) length(unique(x)))
bad_col <- rep(0, ncol(big_data))
inds <- (1:ncol(big_data))[uval == 2]
for(i in 1:length(inds)){
  for(j in i:length(inds)){
    if(i != j & bad_col[inds[i]] == 0 & bad_col[inds[j]] == 0){
      if(jaccard(big_data[,i], big_data[,j]) > 0.9){
        bad_col[inds[j]] <- 1
      }
    }
  }
}
cor_col <- rep(0, ncol(big_data))
inds <- (1:ncol(big_data))[uval > 2]
for(i in 1:length(inds)){
  for(j in i:length(inds)){
    if(i != j & cor_col[inds[i]] == 0 & cor_col[inds[j]] == 0){
      if(cor(big_data[,i], big_data[,j]) > 0.9){
        cor_col[inds[j]] <- 1
      }
    }
  }
}


big_data <- big_data[,bad_col + cor_col == 0]

big_data <- big_data[,colnames(big_data) %in% colnames(test_data)]
test_data <- test_data[,colnames(test_data) %in% colnames(big_data)]


big_df <- big_data


resid_groups <- readRDS(paste0("out_batch/resid_groups.train.", author, ".RDS"))
#extreme_groups <- readRDS(paste0("out_batch/extreme_groups.train.", author, ".RDS"))
#wrong_groups <- readRDS(paste0("out_batch/wrong_groups.train.", author, ".RDS"))
time_groups <- readRDS("out_batch/timeline_groups.train.RDS")

all_case_groups <- readRDS(paste0("out_batch/case_groups.train.", author, ".RDS"))
all_control_groups <- readRDS(paste0("out_batch/control_groups.train.", author, ".RDS"))
extreme_case_groups <- readRDS(paste0("out_batch/case_extreme_groups.train.", author, ".RDS"))
extreme_control_groups <- readRDS(paste0("out_batch/control_extreme_groups.train.", author, ".RDS"))

test_resid_groups <- readRDS(paste0("out_batch/resid_groups.test.", author, ".RDS"))
#test_extreme_groups <- readRDS(paste0("out_batch/extreme_groups.test.", author, ".RDS"))
#test_wrong_groups <- readRDS(paste0("out_batch/wrong_groups.test.", author, ".RDS"))
test_time_groups <- readRDS("out_batch/timeline_groups.test.RDS")

test_all_case_groups <- readRDS(paste0("out_batch/case_groups.test.", author, ".RDS"))
test_all_control_groups <- readRDS(paste0("out_batch/control_groups.test.", author, ".RDS"))
test_extreme_case_groups <- readRDS(paste0("out_batch/case_extreme_groups.test.", author, ".RDS"))
test_extreme_control_groups <- readRDS(paste0("out_batch/control_extreme_groups.test.", author, ".RDS"))

all_groups <- list(resid_groups, time_groups, all_case_groups, all_control_groups, extreme_case_groups, extreme_control_groups)
saveRDS(all_groups, "all_groups.RDS")
all_actions <- c("remove", "remove", "case", "control", "case", "control")
big_switch <- "all_remove"
#----------------------------------------------------------------







#------------------------------------------------------------------------
#  AUC GETTING FUNCTION
#----------------------------------------------------------------------------------

#need to write function to do cross validation
#Need to set up repeated cross validation
run_cv <- function(input_df, remove_eids, the_action, group_ind, spec_ind, input_folds = 3, input_repeats = 3){

  eid <- input_df$eid
  constant_pheno <- input_df$pheno

  all_auc <- rep(0.5, input_folds * input_repeats)
  size_group <- floor(nrow(train_df)/input_folds)
  start_spots <- floor(seq(1, size_group, length.out=input_repeats+1))
  repeat_list <- list()
  k <- 1
  for(i in 1:input_repeats){
    folds_list <- list()
    for(j in 1:input_folds){

      print(paste("CV:", i, j, group_ind))
      test_group <- eid[(start_spots[i]+((j-1)*size_group)):(start_spots[i]+(j*size_group))]
      train_group <- eid[!(eid %in% test_group)]


      if(group_ind  == 1){
        saveRDS(list(train_group, test_group, spec_ind, group_ind, input_df, train_df, big_df), "temp.RDS")
        print("done")
        remove_eids <- predict_residual(train_group, test_group, spec_ind, group_ind, input_df)
      } else if(group_ind == 3){
        remove_eids <- predict_pheno_all(train_group, test_group, spec_ind, group_ind, input_df)
      } else if(group_ind == 4){
        remove_eids <- predict_pheno_all(train_group, test_group, spec_ind, group_ind, input_df)
      } else if(group_ind == 5){
        remove_eids <- predict_pheno_extreme(train_group, test_group, spec_ind, group_ind, input_df)
      } else if(group_ind == 6){
        remove_eids <- predict_pheno_extreme(train_group, test_group, spec_ind, group_ind, input_df)
      }

      print("done getting remove eids")

      if(the_action == "remove" | big_switch == "all_remove"){
        train_group <- train_group[!(train_group %in% remove_eids)]
        test_group <- test_group[!(test_group %in% remove_eids)]

      } else if(the_action == "case"){
        if(big_switch == "explicit"){
          input_df$pheno[input_df$eid %in% remove_eids] <- 1
        } else if(big_switch == "check_pheno"){
          train_group <- train_group[!(train_group %in% remove_eids & input_df$pheno[input_df$eid %in% train_group] == 0)]
          test_group <- test_group[!(test_group %in% remove_eids & input_df$pheno[input_df$eid %in% test_group] == 0)]
        }

      } else if(the_action == "control"){
        if(big_switch == "explicit"){
          input_df$pheno[input_df$eid %in% remove_eids] <- 0
        } else if(big_switch == "check_pheno") {
          train_group <- train_group[!(train_group %in% remove_eids & input_df$pheno[input_df$eid %in% train_group] == 1)]
          test_group <- test_group[!(test_group %in% remove_eids & input_df$pheno[input_df$eid %in% test_group] == 1)]
        }
      }

      if(length(unique(input_df$pheno[input_df$eid %in% test_group])) != 2 |
         length(unique(input_df$pheno[input_df$eid %in% train_group])) != 2){
         all_auc[k] <- 0.5
      } else {
        if("sex" %in% colnames(input_df)){
          mod <- glm(pheno ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + score, data = input_df[input_df$eid %in% train_group,], family = "binomial")
        } else {
          mod <- glm(pheno ~ age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + score, data = input_df[input_df$eid %in% train_group,], family = "binomial")
        }

        preds <- predict(mod, input_df[input_df$eid %in% test_group,])
        auc <- roc(input_df[input_df$eid %in% test_group,]$pheno ~ preds, quiet=T)$auc
        input_df$pheno <- constant_pheno

        #print(auc)
        all_auc[k] <- auc
      }

      k <- k + 1
    }
  }
  return(all_auc)

}



#----------------------------------------------------------------
# Do IT
#-----------------------------------------------------------------------------

base_auc <- mean(run_cv(train_df, 1, "nothing", 0, 0))
saveRDS(base_auc, paste0("res/base_auc.", author, ".RDS"))


step_through <- function(all_groups, start_remove){

  auc_min <- lapply(all_groups, function(x) rep(NA, length(x)))
  auc_mean <- lapply(all_groups, function(x) rep(NA, length(x)))
  auc_best <- lapply(all_groups, function(x) rep(NA, length(x)))


  for(i in 1:length(all_groups)){
    for(j in 1:length(all_groups[[i]])){
      print(paste("big group:", i, "  spec group:", j))

      cvres <- run_cv(train_df, c(all_groups[[i]][[j]], start_remove ), all_actions[i], i, j)
      auc_mean[[i]][j] <- mean(cvres) - base_auc
      auc_min[[i]][j] <- min(cvres) - base_auc
      auc_best[[i]][j] <- sum(cvres > base_auc)
    }
  }

  return_list <- list()
  ii <- 1
  for(auc_perf in list(auc_mean, auc_min, auc_best)){
    #saveRDS(auc_perf, "auc_perf.RDS")
    #saveRDS(all_groups, "all_groups.RDS")

    for(qq in 1:length(auc_perf)){
      if(max(unlist(auc_perf)) %in% auc_perf[[qq]]){
        perm_remove <- c(all_groups[[qq]][which(auc_perf[[qq]] == max(unlist(auc_perf)))], start_remove)
        all_groups[[qq]][which(auc_perf[[qq]] == max(unlist(auc_perf)))] <- 123
        remove_group <- qq
        remove_index <- which(auc_perf[[qq]] == max(unlist(auc_perf)))
      }
    }


    return_list[[ii]] <- list(unlist(auc_perf)[which.max(unlist(auc_perf))] , all_groups, perm_remove, remove_group, remove_index)
    ii <- ii + 1
  }

  return(return_list)
}

#need each of the all_groups to have names

print("start")
back_list_big <- step_through(all_groups, c(123456789))
best_auc_1 <- back_list_big[[1]][[1]]
all_groups_1 <- back_list_big[[1]][[2]]
perm_remove_1 <- back_list_big[[1]][[3]]
remove_group_1 <- back_list_big[[1]][[4]]
remove_index_1 <- back_list_big[[1]][[5]]

print(2)
back_list <- step_through(all_groups_1, perm_remove_1)
best_auc_2 <- back_list[[1]][[1]]
all_groups_2 <- back_list[[1]][[2]]
perm_remove_2 <- back_list[[1]][[3]]
remove_group_2 <- back_list[[1]][[4]]
remove_index_2 <- back_list[[1]][[5]]

print(3)
back_list <- step_through(all_groups_2, perm_remove_2)
best_auc_3 <- back_list[[1]][[1]]
all_groups_3 <- back_list[[1]][[2]]
perm_remove_3 <- back_list[[1]][[3]]
remove_group_3 <- back_list[[1]][[4]]
remove_index_3 <- back_list[[1]][[5]]

#print(4)
#back_list <- step_through(all_groups_3, perm_remove_3)
#best_auc_4 <- back_list[[1]]
#all_groups_4 <- back_list[[2]]
#perm_remove_4 <- back_list[[3]]
#remove_group_4 <- back_list[[4]]
#remove_index_4 <- back_list[[5]]

#print(5)
#back_list <- step_through(all_groups_4, perm_remove_4)
#best_auc_5 <- back_list[[1]]
#all_groups_5 <- back_list[[2]]
#perm_remove_5 <- back_list[[3]]
#remove_group_5 <- back_list[[4]]
#remove_index_5 <- back_list[[5]]



all_remove <- list(perm_remove_1, perm_remove_2, perm_remove_3)
all_auc <- c(best_auc_1, best_auc_2, best_auc_3)
all_remove_group <- c(remove_group_1, remove_group_2, remove_group_3)
all_remove_index <- c(remove_index_1, remove_index_2, remove_index_3)
best_ind <- which.max(all_auc)

all_groups <- readRDS("all_groups.RDS")
all_test_groups <- list(test_resid_groups, test_time_groups, test_all_case_groups, test_all_control_groups, test_extreme_case_groups, test_extreme_control_groups)

all_remove_group <- all_remove_group[1:best_ind]
all_remove_index <- all_remove_index[1:best_ind]
best_test_remove <- list()
best_train_remove <- list()
for(i in 1:length(all_remove_group)){
  best_test_remove[[i]] <- all_test_groups[[all_remove_group[i]]][all_remove_index[i]]
  best_train_remove[[i]] <- all_groups[[all_remove_group[i]]][all_remove_index[i]]
}
#best_test_remove <- unique(unlist(best_test_remove))
#best_train_remove <- unique(unlist(best_train_remove))


###
mean_train_remove <- all_groups[[back_list_big[[1]][[4]]]][back_list_big[[1]][[5]]]
mean_test_remove <- all_test_groups[[back_list_big[[1]][[4]]]][back_list_big[[1]][[5]]]
mean_auc <- back_list_big[[1]][[1]]

min_train_remove <- all_groups[[back_list_big[[2]][[4]]]][back_list_big[[2]][[5]]]
min_test_remove <- all_test_groups[[back_list_big[[2]][[4]]]][back_list_big[[2]][[5]]]
min_auc <- back_list_big[[2]][[1]]

count_train_remove <- all_groups[[back_list_big[[3]][[4]]]][back_list_big[[3]][[5]]]
count_test_remove <- all_test_groups[[back_list_big[[3]][[4]]]][back_list_big[[3]][[5]]]
count_auc <- back_list_big[[3]][[1]]


saveRDS(list(best_train_remove, best_test_remove, all_auc, all_remove_group, all_remove_index), paste0("res/intge.all_around.", author, ".", big_switch, ".RDS"))

saveRDS(list(mean_train_remove, mean_test_remove, mean_auc, back_list_big[[1]][[4]], back_list_big[[1]][[5]]), paste0("res/intge.mean.", author, ".", big_switch, ".RDS"))

saveRDS(list(min_train_remove, min_test_remove, min_auc, back_list_big[[2]][[4]], back_list_big[[2]][[5]]), paste0("res/intge.min.", author, ".", big_switch, ".RDS"))

saveRDS(list(count_train_remove, count_test_remove, count_auc, back_list_big[[3]][[4]], back_list_big[[3]][[5]]), paste0("res/intge.count.", author, ".", big_switch, ".RDS"))
#iterate over all options keeping the best always removed

