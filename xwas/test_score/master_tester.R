library(survival)
library(PRROC)
library(pROC)
library(epitools)

#author <- "Shah" #Liu-2, Malik, Nikpay, Okada, Onengut, Phelan, Rheenen
author <- "lung"
#score_add_on <- "22"
goon = "custom"
score_add_on <- commandArgs(trailingOnly=TRUE)
phen_method <- "icd_selfrep"
score_method <- "auc_best_name"
subrate_style <- "slow"
#train_frac <- 0.6
#test_frac <- 1 - train_frac


full_test_eid <- read.table("~/athena/doc_score/qc/cv_files/test_eid.0.5.txt", stringsAsFactors = F)
full_train_eid <- read.table("~/athena/doc_score/qc/cv_files/train_eid.0.5.txt", stringsAsFactors = F)
tune_eids <- read.table("../tune_score/tune_eids", stringsAsFactors=F)

test_eid <- full_test_eid[!(full_test_eid[,1] %in% tune_eids[,1]),1]
train_eid <- c(full_train_eid[,1], tune_eids[,1])

#THERE MAY BE PROBLEMS WITH SURV_DF

#Read in the PGSs
#Right at the top we read in the train and test eids explicitly
all_scores <-  readRDS(paste0("../do_score/final_scores/all_score_", score_add_on, ".", tolower(author), ".RDS"))
all_eid <- read.table("../do_score/final_eid", header = T, stringsAsFactors=F)
eid <- all_eid[,1]
all_scores <- all_scores[eid %in% train_eid | eid %in% test_eid,]
eid <- eid[eid %in% train_eid | eid %in% test_eid]

#Normalize the scores
best_score <- read.table(paste0("../tune_score/tune_results/", tolower(author), ".best.ss"), stringsAsFactors=F)
best_score <- best_score[best_score[,1] == score_method,2]
best_score <- paste0(author, ".6.clump") #6 is best for cad

scores <- all_scores[,grepl(tolower(author), colnames(all_scores))]
#scores <- apply(scores, 2, function(x) (x-mean(x)) / (max(abs((x-mean(x)))) * 0.01) )
scores <- apply(scores, 2, function(x) (x-min(x)) / ((max(x) - min(x)) * 0.01) )
keep_scores <- scores
scores <- scores[,colnames(scores) == best_score, drop=F]



#Read in the phenotypes, order is: cancer sr, noncancer sr, icd9, icd10, oper, meds
#selfasses date: 2018-11-22; hesin date: 21/01/2000
pheno <- read.table(paste0("../get_pheno/pheno_defs/diag.", tolower(author), ".txt.gz"), stringsAsFactors=F)
dates <- read.table(paste0("../get_pheno/pheno_defs/time.", tolower(author), ".txt.gz"), stringsAsFactors=F)
for(i in 1:ncol(dates)){
  if(i %in% c(1,2,6)){
    dates[dates[,i] == "__________", i] <- "2020-12-31"
    dates[,i] <- as.Date(dates[,i], "%Y-%m-%d")
  } else {
    dates[dates[,i] == "__________", i] <- "31/12/2020"
    dates[,i] <- as.Date(dates[,i], "%d/%m/%Y")
  }
}

if(phen_method == "icd"){
  pheno <- pheno[,3:4]
  dates <- dates[,3:4]
} else if(phen_method == "selfrep"){
  pheno <- pheno[,1:2]
  dates <- dates[,1:2]
} else if(phen_method == "icd_selfrep"){
  pheno <- pheno[,1:4]
  dates <- dates[,1:4]
} else if(phen_method == "all" | phen_method == "double"){
  print("doing nothing")
}

dates <- apply(dates, 1, min)
dates[dates == as.Date("2020-12-31")] <- NA
if(phen_method == "double"){
  pheno <- rowSums(pheno)
  pheno[pheno == 1] <- 0
  pheno[pheno > 1] <- 1

  dates[pheno == 0] <- NA
} else {
  pheno <- rowSums(pheno)
  pheno[pheno > 1] <- 1
}


#Read in the eids used that are the same order as the pheno and dates, then subset the pheno and dates accordingly
pheno_eids <- read.table("../get_pheno/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]
scores <- scores[eid %in% pheno_eids,,drop=F]
eid <- eid[eid %in% pheno_eids]
dates <- dates[pheno_eids %in% eid]
pheno <- pheno[pheno_eids %in% eid]
pheno_eids <- pheno_eids[pheno_eids %in% eid]
dates <- dates[order(pheno_eids)[rank(eid)]]
pheno <- pheno[order(pheno_eids)[rank(eid)]]


#Read in the base covars
covars <- readRDS("~/athena/doc_score/analyze_score/get_covars/base_covars.RDS")
covars <- covars[covars[,1] %in% eid,]
covars <- covars[order(covars[,1])[rank(eid)],]

#Set up survival analysis data frame
#Artifically decide start date is 1999, that way all are even, if date is prior then remove it
#The current maximum date possible is 31 May 2020
death <- read.table("~/athena/ukbiobank/hesin/death.txt", stringsAsFactors=F, header = T)
death[,5] <- unlist(lapply(death[,5], function(x) paste0(strsplit(x, "/")[[1]][3], "-", strsplit(x, "/")[[1]][2], "-", strsplit(x, "/")[[1]][1])))
death <- death[!duplicated(death[,1]),]
death <- death[death[,1] %in% eid,]
add_on <- death[1,]
add_on[5] <- ""
add_eid <- eid[!(eid %in% death[,1])]
add_on <- add_on[rep(1, length(add_eid)),]
add_on$eid <- add_eid
death <- rbind(death, add_on)
death <- death[order(death[,1])[rank(eid)],]

censor <- read.table("~/athena/doc_score/analyze_score/get_covars/covar_data/censor_covars", stringsAsFactors=F, header = T, sep = ",")
if(sum(!(eid %in% censor[,1])) > 0){
  add_on <- matrix(0, nrow = sum(!(eid %in% censor[,1])), ncol = 3)
  add_on[,1] <- eid[!(eid %in% censor[,1])]
  colnames(add_on) <- colnames(censor)
  censor <- rbind(censor, add_on)
}
censor <- censor[censor[,1] %in% eid,]
censor <- censor[order(censor[,1])[rank(eid)],]


start_date <- rep("1999-01-01", nrow(scores))
end_date <- rep("2020-05-31", nrow(scores))
end_date[censor[,2] != ""] <- censor[censor[,2] != "", 2]
end_date[death[,5] != ""] <- death[death[,5] != "", 5]
end_date[!is.na(dates)] <- dates[!is.na(dates)]

is_death_date <- rep(0, nrow(scores))
is_death_date[death[,5] != ""] <- 1
surv_df <- data.frame(time = as.numeric(as.Date(end_date) - as.Date(start_date)), pheno, is_death_date, covars, score = scores[,1])
df <- data.frame(pheno, covars[,-1], score = scores[,1])
#df <- data.frame(pheno, covars[,-1], score = scores)


#need to remove people that had diagnosis before accepted start of the study
eid <- eid[surv_df$time > 0]
df <- df[surv_df$time > 0,]
covars <- covars[surv_df$time > 0,]
pheno <- pheno[surv_df$time > 0]
surv_df <- surv_df[surv_df$time > 0,]




#Subset the sex
base_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")



#Set up Fine and Gray
print("finegray")
event_type <- rep("censor", nrow(surv_df))
event_type[surv_df$pheno == 1] <- "diagnosis"
event_type[surv_df$is_death_date == 1] <- "death"
surv_df$event_type <- as.factor(event_type)
surv_df$eid <- eid
fg_diag <- finegray(Surv(time, event_type) ~ ., data = surv_df[,-which(colnames(surv_df) %in% c("is_death_date", "pheno"))], etype="diagnosis")
fg_death <- finegray(Surv(time, event_type) ~ ., data = surv_df[,-which(colnames(surv_df) %in% c("is_death_date", "pheno"))], etype="death")
fg_scores_diag <- fg_diag[,grepl(tolower(author), colnames(fg_diag))]
fg_scores_death <- fg_death[,grepl(tolower(author), colnames(fg_death))]
fg_diag <- fg_diag[,!grepl("ss", colnames(fg_diag))]
fg_death <- fg_death[,!grepl("ss", colnames(fg_death))]


go_small <- TRUE

df_train <- df[eid %in% train_eid,]
df_test <- df[eid %in% test_eid,]
if(go_small){
  sub_test_eid <- eid[eid %in% test_eid]
  case_eid <- sub_test_eid[df_test$pheno == 1]
  control_eid <- sub_test_eid[df_test$pheno == 0]
  new_test_eid <- c(sample(case_eid, round(length(case_eid)/2)), sample(control_eid, round(length(control_eid)/2)))
  new_train_eid <- c(case_eid[!(case_eid %in% new_test_eid)], control_eid[!(control_eid %in% new_test_eid)])
  df_train <- df_test[sub_test_eid %in% new_train_eid,]
  df_test <- df_test[sub_test_eid %in% new_test_eid,]

  surv_df_train <- surv_df[eid %in% new_train_eid,]
  surv_df_test <- surv_df[eid %in% new_test_eid,]
  fg_diag_train <- fg_diag[fg_diag$eid %in% new_train_eid,]
  fg_diag_test <- fg_diag[fg_diag$eid %in% new_test_eid,]
  fg_death_train <- fg_death[fg_death$eid %in% new_train_eid,]
  fg_death_test <- fg_death[fg_death$eid %in% new_test_eid,]

} else {

  surv_df_train <- surv_df[eid %in% train_eid,]
  surv_df_test <- surv_df[eid %in% test_eid,]
  fg_diag_train <- fg_diag[fg_diag$eid %in% train_eid,]
  fg_diag_test <- fg_diag[fg_diag$eid %in% test_eid,]
  fg_death_train <- fg_death[fg_death$eid %in% train_eid,]
  fg_death_test <- fg_death[fg_death$eid %in% test_eid,]

}

#exit()


    
    ###########################################################
    #                    SURVIVAL ANALYSES                    #
    ###########################################################
print("base surv")
    
    #BASE ###############################
    base_mod <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", base_covars)), data = fg_diag_train, weight = fgwt)
    
    base_conc_obj <- survConcordance(Surv(fgstart, fgstop, fgstatus) ~ predict(base_mod, fg_diag_test), data = fg_diag_test, weight = fgwt)
    base_conc <- c(base_conc_obj$conc, base_conc_obj$std.err)
    

    if(goon == "survfit"){
      base_survfit <- survfit(base_mod, fg_diag_test, se.fit = F)
      base_full_cumhaz <- base_survfit$cumhaz[,!duplicated(fg_diag_test$eid)]
      final_cumhaz <- base_full_cumhaz[nrow(base_survfit$cumhaz),]
      the_time <- base_survfit$time
    
    
    } else {
      H0 <- basehaz(base_mod, centered = FALSE)
      coef <- base_mod$coefficients
      good_order <- colnames(surv_df_test)[colnames(surv_df_test) %in% names(coef)]
      coef <- coef[order(names(coef))[rank(good_order)]]
      prop_haz <- rowSums(t(t(surv_df_test[,colnames(surv_df_test) %in% names(coef)]) * coef))

      base_full_cumhaz <- matrix(0, nrow = length(H0$time), ncol = length(prop_haz))
      for(k in 1:nrow(base_full_cumhaz)){
        base_full_cumhaz[k,] <- H0$hazard[k] * exp(prop_haz)
      }

      final_cumhaz <- base_full_cumhaz[nrow(base_full_cumhaz),]
      the_time <- H0$time
    }
    
    group_factor <- rep(1, length(final_cumhaz))
    group_factor[final_cumhaz < quantile(final_cumhaz, 0.2)] <- 0
    group_factor[final_cumhaz > quantile(final_cumhaz, 0.8)] <- 2
    
    base_avg_cumhaz <- data.frame(time = the_time,
                                  mean_lo = apply(base_full_cumhaz[,group_factor==0], 1, mean),
                                  mean_mid = apply(base_full_cumhaz[,group_factor==1], 1, mean),
                                  mean_hi = apply(base_full_cumhaz[,group_factor==2], 1, mean),
                                  sd_lo = apply(base_full_cumhaz[,group_factor==0], 1, sd),
                                  sd_mid = apply(base_full_cumhaz[,group_factor==1], 1, sd),
                                  sd_hi = apply(base_full_cumhaz[,group_factor==2], 1, sd))





      #WITH SCORE ###########################
print("score surv")

      #Make the model
      score_mod <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", base_covars, " + score")), data = fg_diag_train, weight = fgwt)
      
      #CONCORDANCE ###################################
      score_conc_obj <- survConcordance(Surv(fgstart, fgstop, fgstatus) ~ predict(score_mod, fg_diag_test), data = fg_diag_test, weight = fgwt)
      score_conc <- c(score_conc_obj$conc, score_conc_obj$std.err)

     
      #SURVFIT ##########################################
      #SLOW
      print("slow")
      if(subrate_style == "slow"){

	if(goon == "survfit"){
          all_cumhaz <- list()
          all_time <- list()
          start_vals <- seq(1, nrow(surv_df_test), 1000)
          end_vals <- c(start_vals[-1]-1, nrow(surv_df_test))
          for(k in 1:length(start_vals)){
            g1 <- surv_df_test$eid[start_vals[k]:end_vals[k]]
            score_survfit <- survfit(score_mod, fg_diag_test[fg_diag_test$eid %in% g1,], se.fit = F)
            all_cumhaz[[k]] <- score_survfit$cumhaz[,!duplicated(fg_diag_test$eid[fg_diag_test$eid %in% g1])]
            all_time[[k]] <- score_survfit$time
          }
          score_full_cumhaz <- do.call("cbind", all_cumhaz)

          final_cumhaz <- score_full_cumhaz[nrow(score_survfit$cumhaz),]
          the_time <- score_survfit$time

        } else {
          H0 <- basehaz(score_mod, centered = FALSE)
          coef <- score_mod$coefficients
          good_order <- colnames(surv_df_test)[colnames(surv_df_test) %in% names(coef)]
          coef <- coef[order(names(coef))[rank(good_order)]]
          prop_haz <- rowSums(t(t(surv_df_test[,colnames(surv_df_test) %in% names(coef)]) * coef))
          
          score_full_cumhaz <- matrix(0, nrow = length(H0$time), ncol = length(prop_haz))
          for(k in 1:nrow(score_full_cumhaz)){
            score_full_cumhaz[k,] <- H0$hazard[k] * exp(prop_haz)
          }

          final_cumhaz <- score_full_cumhaz[nrow(score_full_cumhaz),]
          the_time <- H0$time
        }
        
        group_factor <- rep(1, length(final_cumhaz))
        group_factor[final_cumhaz < quantile(final_cumhaz, 0.2)] <- 0
        group_factor[final_cumhaz > quantile(final_cumhaz, 0.8)] <- 2
        
        score_avg_cumhaz <- data.frame(time = the_time,
                                        mean_lo = apply(score_full_cumhaz[,group_factor==0], 1, mean),
                                        mean_mid = apply(score_full_cumhaz[,group_factor==1], 1, mean),
                                        mean_hi = apply(score_full_cumhaz[,group_factor==2], 1, mean),
                                        sd_lo = apply(score_full_cumhaz[,group_factor==0], 1, sd),
                                        sd_mid = apply(score_full_cumhaz[,group_factor==1], 1, sd),
                                        sd_hi = apply(score_full_cumhaz[,group_factor==2], 1, sd))
        score_avg_cumhaz <- score_avg_cumhaz[!duplicated(score_avg_cumhaz$mean_mid),]



      } else if(subrate_style == "fast"){ 
        #FAST
        group_list <- list()
        for(k in 1:length(names(score_mod$coef))){
         group_list[[k]] <-  as.numeric(quantile(surv_df_train[[names(score_mod$coef)[k]]], c(0.1, 0.5, 0.9)))
         if(sign(score_mod$coef[k] == -1)){
          group_list[[k]] <- rev(group_list[[k]])
         }
        }
        mean_df <- data.frame(do.call("cbind", group_list))
        colnames(mean_df) <- names(score_mod$coef)
        if(sign(score_mod$coef[2]) == -1){
         mean_df$sex <- c(0.9, 0.5, 0.1)
        } else {
         mean_df$sex <- c(0.1, 0.5, 0.9)
        }
        
        score_pred <- survfit(score_mod, newdata = mean_df)
        score_full_cumhaz <- NULL
        
        score_avg_cumhaz <- data.frame(score_pred$time, score_pred$cumhaz, score_pred$std.err)
        colnames(score_avg_cumhaz) <- c("time", "mean_lo", "mean_mid", "mean_hi", "sd_lo", "sd_mid", "sd_hi")
        score_avg_cumhaz <- score_avg_cumhaz[!duplicated(score_avg_cumhaz$mean_mid),]
      }
      
      


    # Set up the normal model data
    print("normal")
  
    ###########################################################
    #                    NORMAL MODELS                        #
    ###########################################################

    base_mod <- glm(as.formula(paste0("pheno ~ ", base_covars)), data = df_train, family = "binomial")
    base_pred <- predict(base_mod, df_test)
    base_roc <- roc(df_test$pheno ~ base_pred, quiet=T)
    base_auc <- as.numeric(ci.auc(base_roc))
    base_pr <- pr.curve(scores.class0 = base_pred, weights.class0 = df_test$pheno, curve = T)



    all_base_odds_ratio <- matrix(0, nrow = 6, ncol = 3)
    k <- 1
    for(cut_off in c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995)){
      base_group <- rep(1, length(base_pred))
      base_group[base_pred < quantile(base_pred, 0.2)] <- 0
      base_group[base_pred > quantile(base_pred, cut_off)] <- 2
      base_odds_table <- matrix(c(sum(df_test$pheno == 1 & base_group == 2), sum(df_test$pheno == 0 & base_group == 2),
                         sum(df_test$pheno == 1 & base_group == 0), sum(df_test$pheno == 0 & base_group == 0)), nrow = 2)
      base_odds_ratio <- oddsratio.wald(base_odds_table)
      all_base_odds_ratio[k,] <- base_odds_ratio$measure[2,c(2,1,3)]

      k <- k + 1
    }
 

#exit()
      #Make the model
      score_mod <- glm(as.formula(paste0("pheno ~", base_covars, " + score")), data = df_train, family = "binomial")
      score_pred <- predict(score_mod, df_test)

      # AUC #####################################
      score_roc <- roc(df_test$pheno ~ score_pred, quiet=T)
      score_auc <- as.numeric(ci.auc(score_roc))

      # PR #####################################
      score_pr <- pr.curve(scores.class0 = score_pred, weights.class0 = df$pheno, curve = T)


      # ODDS RATIO #############################
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


        k <- k + 1
      }
    

#Get answers together ##################################################
#previously also held the base_full_cumhaz but memory gets to be alot
all_base_holder <- list("conc" = base_conc, "survfit" = base_avg_cumhaz, "auc" = list(base_roc, base_auc), "or" = all_base_odds_ratio, "pr" = base_pr)
all_extra_holder <- list(NULL)
all_other_holder <- list(NULL)

misc_info <- list("phen_method" = phen_method, "subrate_style" = subrate_style, "score_method" = score_method)
final_obj <- list("conc" = score_conc, "survfit" = score_avg_cumhaz, "auc" = list(score_roc, score_auc), "or" = all_score_odds_ratio, "pr" = score_pr, "base" = all_base_holder, "score_names" = best_score, "misc" = misc_info)
saveRDS(final_obj, paste0("final_stats/", author, ".", score_add_on, ".RDS"))


