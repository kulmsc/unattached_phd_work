library(survival)
library(PRROC)
library(pROC)
library(epitools)

args <- commandArgs(trailingOnly=TRUE)
#author <- "christophersen"
author <- args[1]
batch_type <- "mean"
#batch_type <- args[2]
case_type <- "all_remove"
#case_type <- args[3]

res <- readRDS(paste0("../batch_indivs/res/togo.", batch_type,".", author, ".", case_type,".RDS"))

train_df <- readRDS(paste0("../batch_indivs/data/use_train_df.", author, ".RDS"))
test_df <- readRDS(paste0("../batch_indivs/data/use_test_df.", author, ".RDS"))

better_train_df <- data.frame(train_df)
better_test_df <- data.frame(test_df)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
all_actions <- c("remove", "remove", "case", "control", "case", "control")
big_switch <- case_type

train_remove_eids <- unlist(res[[1]])
test_remove_eids <- unlist(res[[2]])
the_action <- all_actions[res[[4]]]

if(the_action == "remove" | big_switch == "all_remove"){
  better_train_df <- better_train_df[!(train_df$eid %in% train_remove_eids),]
  better_test_df <- better_test_df[!(test_df$eid %in% test_remove_eids),]

} else if(the_action == "case"){
  if(big_switch == "explicit"){
    better_train_df$pheno[better_train_df$eid %in% train_remove_eids] <- 1
    better_test_df$pheno[better_test_df$eid %in% test_remove_eids] <- 1
  } else if(big_switch == "check_pheno"){
    better_train_df <- better_train_df[!(better_train_df$eid %in% train_remove_eids & better_train_df$pheno == 0),]
    better_test_df <- better_test_df[!(better_test_df$eid %in% test_remove_eids & better_test_df$pheno == 0),]
  }

} else if(the_action == "control"){
  if(big_switch == "explicit"){
    better_train_df$pheno[better_train_df$eid %in% train_remove_eids] <- 1
    better_test_df$pheno[better_test_df$eid %in% test_remove_eids] <- 1
  } else if(big_switch == "check_pheno") {
    better_train_df <- better_train_df[!(better_train_df$eid %in% train_remove_eids & better_train_df$pheno == 1),]
    better_test_df <- better_test_df[!(better_test_df$eid %in% test_remove_eids & better_test_df$pheno == 1),]
  }
}




#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


poss_covars <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_covar", stringsAsFactors=F, sep = "\t")
poss_covars[,1] <- tolower(poss_covars[,1])
poss_author <- gsub("-", "", author)
extra_covars <- strsplit(poss_covars[poss_covars[,1] == poss_author,2], ",")[[1]]
extra_covar <- gsub(" ", "_", extra_covars)
extra_covar <- colnames(train_df)[colnames(train_df) %in% extra_covar & colnames(train_df) != "age"]

#saveRDS(list(best_train_remove, best_test_remove, all_auc, all_remove_group, all_remove_index), paste0("res/perf.", author, ".RDS"))


#-------------------------------------------------------
# THE FUNCTION
#------------------------------------------------------------------
if("sex" %in% colnames(train_df)){
  base_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
} else {
  base_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
}

few_names <- colnames(train_df)[!(colnames(train_df) %in% c("pheno", "score", "eid", "age", "sex",  "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))]
if(sum(!grepl("PGS", few_names)) > 0){
  run_extra_covar <- TRUE
} else {
  run_extra_covar <- FALSE
}

if(any(grepl("PGS", colnames(train_df)))){
  run_other_scores <- TRUE
} else {
  run_other_scores <- FALSE
}

other_scores <- train_df[,grepl("PGS", colnames(train_df)),drop = F]


get_over_here <- function(df_train, df_test, df_type){
    base_mod <- glm(as.formula(paste0("pheno ~ ", base_covars)), data = df_train, family = "binomial")
    base_pred <- predict(base_mod, df_test)
    base_roc <- roc(df_test$pheno ~ base_pred, quiet=T)
    base_auc <- as.numeric(ci.auc(base_roc))
    base_pr <- pr.curve(scores.class0 = base_pred, weights.class0 = df_test$pheno, curve = T)

    # NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
    # if there are additional covariates then I should create a model that includes these covariates
    # then with this model form a complementary prediction from which ROC and PR curves can be drawn
    # this is the same process that happens in the odds ratio for loop after if(run_extra_covar) ...
    # NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
    if(run_extra_covar){
      df_train_extra <- df_train[complete.cases(df_train),]
      df_test_extra <- df_test[complete.cases(df_test),]
      extra_base_mod <- glm(as.formula(paste0("pheno ~ ", base_covars, "+", paste(extra_covar, collapse = "+"))), data = df_train_extra, family = "binomial")
      extra_base_pred <- predict(extra_base_mod, df_test_extra)
      extra_base_roc <- roc(df_test_extra$pheno ~ extra_base_pred, quiet=T)
      extra_base_auc <- as.numeric(ci.auc(extra_base_roc))
      extra_base_pr <- pr.curve(scores.class0 = extra_base_pred, weights.class0 = df_test_extra$pheno, curve = T)
    }


    all_base_odds_ratio <- matrix(0, nrow = 6, ncol = 3)
    all_extra_base_odds_ratio <- matrix(0, nrow = 6, ncol = 3)
    k <- 1
    for(cut_off in c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995)){
      base_group <- rep(1, length(base_pred))
      base_group[base_pred < quantile(base_pred, 0.2)] <- 0
      base_group[base_pred > quantile(base_pred, cut_off)] <- 2
      base_odds_table <- matrix(c(sum(df_test$pheno == 1 & base_group == 2), sum(df_test$pheno == 0 & base_group == 2),
                         sum(df_test$pheno == 1 & base_group == 0), sum(df_test$pheno == 0 & base_group == 0)), nrow = 2)
      base_odds_ratio <- oddsratio.wald(base_odds_table)
      all_base_odds_ratio[k,] <- base_odds_ratio$measure[2,c(2,1,3)]

      if(run_extra_covar){
        base_group <- rep(1, length(extra_base_pred))
        base_group[extra_base_pred < quantile(extra_base_pred, 0.2)] <- 0
        base_group[extra_base_pred > quantile(extra_base_pred, cut_off)] <- 2

        base_odds_table <- matrix(c(sum(df_test_extra$pheno == 1 & base_group == 2), sum(df_test_extra$pheno == 0 & base_group == 2),
                         sum(df_test_extra$pheno == 1 & base_group == 0), sum(df_test_extra$pheno == 0 & base_group == 0)), nrow = 2)
        base_odds_ratio <- oddsratio.wald(base_odds_table)
        all_extra_base_odds_ratio[k,] <- base_odds_ratio$measure[2,c(2,1,3)]
      }
      k <- k + 1
    }

      #Make the model
      score_mod <- glm(as.formula(paste0("pheno ~", base_covars, " + score")), data = df_train, family = "binomial")
      score_pred <- predict(score_mod, df_test)
      # NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
      # just here as above we create models similar to the basic score model, although now we need additional covariates
      # the first if statements will add extra covariates to this mode, the second statement adds the extra scores
      # This series if statements will come up again several times after each type of statistic is prepared (AUC, OR, PR)
      # NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
      if(run_extra_covar){
        extra_score_mod <- glm(as.formula(paste0("pheno ~ ", base_covars, " + score +", paste(extra_covar, collapse = "+"))), data = df_train_extra, family = "binomial")
        extra_score_pred <- predict(extra_score_mod, df_test_extra)
      }
      if(run_other_scores){
        all_other_score_mod <- list()
        all_other_score_pred <- list()
        all_other_score_roc <- list()
        all_other_score_auc <- list()
        all_other_score_pr <- list()
        for(k in 1:ncol(other_scores)){
          all_other_score_mod[[k]] <- glm(as.formula(paste0("pheno ~ ", base_covars, " + ", colnames(other_scores)[k])), data = df_train, family = "binomial")
          all_other_score_pred[[k]] <- predict(all_other_score_mod[[k]], df_test)
        }
      }

      # AUC #####################################
      score_roc <- roc(df_test$pheno ~ score_pred, quiet=T)
      score_auc <- as.numeric(ci.auc(score_roc))
      if(run_extra_covar){
        extra_score_roc <- roc(df_test_extra$pheno ~ extra_score_pred, quiet=T)
        extra_score_auc <- as.numeric(ci.auc(extra_score_roc))
      }
      if(run_other_scores){
        for(k in 1:ncol(other_scores)){
          all_other_score_roc[[k]] <- roc(df_test$pheno ~ all_other_score_pred[[k]], quiet=T)
          all_other_score_auc[[k]] <- as.numeric(ci.auc(all_other_score_roc[[k]]))
        }
      }

      # PR #####################################
      score_pr <- pr.curve(scores.class0 = score_pred, weights.class0 = df_test$pheno, curve = T)
      if(run_extra_covar){
        extra_score_pr <- pr.curve(scores.class0 = extra_score_pred, weights.class0 = df_test_extra$pheno, curve = T)
      }
      if(run_other_scores){
        for(k in 1:ncol(other_scores)){
          all_other_score_pr[[k]] <- pr.curve(scores.class0 = all_other_score_pred[[k]], weights.class0 = df_test$pheno, curve = T)
        }
      }


      # ODDS RATIO #############################
      all_other_score_odds_ratio <- rep(list(matrix(0, nrow = 6, ncol = 3)), ncol(other_scores))
      all_extra_score_odds_ratio <- matrix(0, nrow = 6, ncol = 3)
      all_score_odds_ratio <- matrix(0, nrow = 6, ncol = 3)
      k <- 1
      for(cut_off in c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995)){
        score_group <- rep(1, length(score_pred))
        score_group[score_pred < quantile(score_pred, 0.2)] <- 0
        score_group[score_pred > quantile(score_pred, cut_off)] <- 2

        score_odds_table <- matrix(c(sum(df_test$pheno == 1 & score_group == 2), sum(df_test$pheno == 0 & score_group == 2),
                         sum(df_test$pheno == 1 & score_group == 0), sum(df_test$pheno == 0 & score_group == 0)), nrow = 2)
        score_odds_ratio <- oddsratio.wald(score_odds_table)
        all_score_odds_ratio[k,] <- score_odds_ratio$measure[2,c(2,1,3)]

        if(run_extra_covar){
          score_group <- rep(1, length(extra_score_pred))
          score_group[extra_score_pred < quantile(extra_score_pred, 0.2)] <- 0
          score_group[extra_score_pred > quantile(extra_score_pred, cut_off)] <- 2

          score_odds_table <- matrix(c(sum(df_test_extra$pheno == 1 & score_group == 2), sum(df_test_extra$pheno == 0 & score_group == 2),
                         sum(df_test_extra$pheno == 1 & score_group == 0), sum(df_test_extra$pheno == 0 & score_group == 0)), nrow = 2)
          score_odds_ratio <- oddsratio.wald(score_odds_table)
          all_extra_score_odds_ratio[k,] <- score_odds_ratio$measure[2,c(2,1,3)]
        }

        if(run_other_scores){
          for(l in 1:ncol(other_scores)){
            score_group <- rep(1, length(all_other_score_pred[[l]]))
            score_group[all_other_score_pred[[l]] < quantile(all_other_score_pred[[l]], 0.2)] <- 0
            score_group[all_other_score_pred[[l]] > quantile(all_other_score_pred[[l]], cut_off)] <- 2

            score_odds_table <- matrix(c(sum(df_test$pheno == 1 & score_group == 2), sum(df_test$pheno == 0 & score_group == 2),
                           sum(df_test$pheno == 1 & score_group == 0), sum(df_test$pheno == 0 & score_group == 0)), nrow = 2)
            score_odds_ratio <- oddsratio.wald(score_odds_table)
            all_other_score_odds_ratio[[l]][k,] <- score_odds_ratio$measure[2,c(2,1,3)]
          }
        }

     k <- k + 1
  }

  all_base_holder <- list("auc" = list(base_roc, base_auc), "or" = all_base_odds_ratio, "pr" = base_pr)
  all_score_holder <- list("auc" = list(score_roc, score_auc), "or" = all_score_odds_ratio, "pr" = score_pr)

  if(run_extra_covar){
    all_extra_holder <- list("base_auc" = list(extra_base_roc, extra_base_auc), "base_pr" = extra_base_pr, "base_or" = all_extra_base_odds_ratio, "auc" = list(extra_score_roc, extra_score_auc), "pr" = extra_score_pr, "or" = all_extra_score_odds_ratio)
  } else {
    all_extra_holder <- list(NULL)
  }

  if(run_other_scores){
    all_other_holder <- list("auc" = list(all_other_score_roc, all_other_score_auc), "pr" = all_other_score_pr, "or" = all_other_score_odds_ratio)
  } else {
    all_other_holder <- list(NULL)
  }

  final_holder <- list("base" = all_base_holder, "score" = all_score_holder, "extra" = all_extra_holder, "other" = all_other_holder)
  #return(final_holder)
  saveRDS(final_holder, paste0("acc_results/", df_type, ".", batch_type, ".", case_type, ".", author, ".RDS"))
}

get_over_here(train_df, test_df, "all_indiv")
get_over_here(better_train_df, better_test_df, "better_indiv")
