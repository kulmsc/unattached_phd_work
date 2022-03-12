library(vroom)
library(Hmisc)
library(glmnet)
library(data.table)

#author <- "michailidou"
author <- commandArgs(trailingOnly=TRUE)[1]

############################# READ IN AND SET UP #########################################

train_df <- readRDS(paste0("../get_data/df_output/df_train.", author, ".RDS"))
test_df <- readRDS(paste0("../get_data/df_output/df_test.", author, ".RDS"))
surv_train_df <- readRDS(paste0("../get_data/df_output/survdf_train.", author, ".RDS"))
surv_test_df <- readRDS(paste0("../get_data/df_output/survdf_test.", author, ".RDS"))

big_data <- readRDS("../get_data/big_data.RDS")
colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))] <- paste0("X", colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))])

big_data <- big_data[big_data$eid %in% train_df$eid,]
train_df <- train_df[train_df$eid %in% big_data$eid,]
big_data <- big_data[order(big_data$eid)[rank(train_df$eid)],]


if(mean(train_df$score[train_df$pheno == 1]) < mean(train_df$score[train_df$pheno == 0])){
  train_df$score <- train_df$score * -1
}


#need to remove phenotypes that occured before enrollment
#note that the big_data only has ICDs, etc. that occured before enrollment
#can compare phenotypes such that where a big data col matches with the phenotype call we remove the person from the analysis
#then also remove all relevant cols just to make sure



pheno_decoder <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_defs", stringsAsFactors=F, header=T)
useful_decoder <- as.character(pheno_decoder[tolower(pheno_decoder[,1]) == author,4:9][1,])
real_names <- c("CANCER", "NONCANCER", "ICD9", "ICD10", "OPCS", "MEDS")
bad_names <- c()
bad_people <- list()
j <- 1

#remove individuals who were already diagnosed
for(i in 1:length(useful_decoder)){
  if(!is.na(useful_decoder[i])){
    if(grepl("|", useful_decoder[i])){
      for(subcode in strsplit(useful_decoder[i], "|", fixed = T)[[1]]){
        bad_names <- c(bad_names, paste0(real_names[i], "_", subcode))   
        if(real_names[i] %in% c("ICD9", "ICD10", "OPCS")){
          bad_people[[j]] <- big_data$eid[big_data[,which(colnames(big_data) == paste0(real_names[i], "_", subcode))] == 1]
          j <- j + 1
        }
      }

    } else {
      bad_names <- c(bad_names, paste0(real_names[i], "_", useful_decoder[i]))
      bad_people[[j]] <- big_data$eid[big_data[,which(colnames(big_data) == paste0(real_names[i], "_", useful_decoder[i]))] == 1]
      j <- j + 1
    }
    
  }
}

bad_people <- unique(unlist(bad_people))
big_data <- big_data[!(big_data$eid %in% bad_people),]
train_df <- train_df[!(train_df$eid %in% bad_people),]

big_data <- big_data[,-which(colnames(big_data) %in% bad_names)]
minor_counts <- readRDS("minor_counts.RDS")
big_data <- big_data[,!(colnames(big_data) %in% names(minor_counts)[minor_counts < 250])]

na_counts <- apply(big_data, 2, function(x) sum(is.na(x)))
big_data <- big_data[,na_counts < nrow(big_data)*0.9]


### IMPUTATION ###
na_count <- apply(big_data, 2, function(x) sum(is.na(x)))
#what common variable should be included: age, sex, 
inc_vars <- c("SURVEY_X31.0.0", #sex
              "SURVEY_X34.0.0", #year of birth
              "SURVEY_X50.0.0", #standing height
              "SURVEY_X135.0.0", #number illnesses
              "SURVEY_X137.0.0", #number of medication
              "SURVEY_X1647.0.0", #born in UK
              "SURVEY_X2178.0.0", #health ratio
              "SURVEY_X20161.0.0", #pack years smoking
              "SURVEY_X20458.0.0", #general happiness
              "SURVEY_X21002.0.0", #weight
              "SURVEY_X22599.0.0", #jobs held
              "SURVEY_X23104.0.0", #BMI
              "SURVEY_X30020.0.0") #hemoglobin

if(length(unique(big_data$SURVEY_X31.0.0)) == 1){
  big_data <- big_data[,-which(colnames(big_data) == "SURVEY_X31.0.0")]
  inc_vars <- inc_vars[inc_vars != "SURVEY_X31.0.0"]
}



the_formula <- as.formula(paste0("~ ", paste(inc_vars, collapse = " + ")))
f <- aregImpute(the_formula, nk=0, tlinear=FALSE, data = big_data, B = 50, n.impute = 1)
for(i in 1:length(inc_vars)){
  if(!is.null(f$imputed[[i]]) & any(is.na(big_data[inc_vars[i]]))){
    big_data[inc_vars[i]][is.na(big_data[inc_vars[i]])] <- unname(f$imputed[[inc_vars[i]]][,1])
  }
}

bad_names <- names(na_count)[na_count > 0]
for(i in 1:length(bad_names)){
  the_formula <- as.formula(paste0("~ ", paste(c(inc_vars, bad_names[i]), collapse = " + ")))
  f <- aregImpute(the_formula, nk=0, tlinear=FALSE, data = big_data, B = 50, n.impute = 1)
  big_data[bad_names[i]][is.na(big_data[bad_names[i]])] <- unname(f$imputed[[bad_names[i]]][,1])
}

u_counts <- apply(big_data, 2, function(x) length(unique(x)))
big_data <- big_data[,u_counts > 1]

for(i in 1:ncol(big_data)){
  big_data[,i] <- (big_data[,i] - min(big_data[,i]))/(max(big_data[,i]) - min(big_data[,i]))
}


#print("REMEMBER TO SAVE BIG_DATA")
#exit()
################################### RESIDUALS  ###########################################
if("sex" %in% colnames(train_df)){
  base_mod <- glm(pheno ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, train_df, family = "binomial")
} else {
  base_mod <- glm(pheno ~ age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, train_df, family = "binomial")
}
resids <- base_mod$residuals

count_zero_resid <- function(in_lambda){
  mod <- glmnet(x = as.matrix(big_data), y = resids, alpha = 1, lambda = in_lambda)
  zero_count <- sum(as.numeric(mod$beta) != 0) - 10
  return(zero_count)
}

root_res <- uniroot(count_zero_resid, lower = 0.00001, upper = 15)
mod <- glmnet(x = as.matrix(big_data), y = resids, alpha = 1, lambda = root_res$root)

keeper_names_resid <- rownames(as.matrix(mod$beta))[as.numeric(mod$beta) != 0]

#############

count_zero_pheno <- function(in_lambda){
  mod <- glmnet(x = as.matrix(big_data), y = train_df$pheno, alpha = 1, lambda = in_lambda)
  zero_count <- sum(as.numeric(mod$beta) != 0) - 10
  return(zero_count)
}

root_res <- uniroot(count_zero_pheno, lower = 0.00001, upper = 15)
mod <- glmnet(x = as.matrix(big_data), y = train_df$pheno, alpha = 1, lambda = root_res$root)

keeper_names_pheno <- rownames(as.matrix(mod$beta))[as.numeric(mod$beta) != 0]





#exit()

########
if(length(keeper_names_resid) < 10){
  keeper_names_resid <- c(keeper_names_resid, replicate(10-length(keeper_names_resid), paste0(sample(letters, 5), collapse="")))
}
if(length(keeper_names_pheno) < 10){
  keeper_names_pheno <- c(keeper_names_pheno, replicate(10-length(keeper_names_pheno), paste0(sample(letters, 5), collapse="")))
}
if(length(keeper_names_resid) > 10){
  keeper_names_resid <- keeper_names_resid[1:10]
}
if(length(keeper_names_pheno) > 10){
  keeper_names_pheno <- keeper_names_pheno[1:10]
}

save_coef <- data.frame("resid" = keeper_names_resid, "pheno" = keeper_names_pheno)


saveRDS(save_coef, paste0("resid_res/", author, ".new.RDS"))


