library(survival)
library(PRROC)
library(pROC)
library(epitools)

#author <- "Shah" #Liu-2, Malik, Nikpay, Okada, Onengut, Phelan, Rheenen
#author <- "Michailidou"
author <- commandArgs(trailingOnly=TRUE)
phen_method <- "icd_selfrep"
score_method <- "auc_best_name"
subrate_style <- "slow"
train_frac <- 0.6
test_frac <- 1 - train_frac




#THERE MAY BE PROBLEMS WITH SURV_DF

#Read in the PGSs
#Right at the top we read in the train and test eids explicitly
all_scores <-  readRDS(paste0("~/athena/doc_score/do_score/final_scores/all_score.", tolower(author), ".RDS"))
all_eid <- read.table("~/athena/doc_score/do_score/all_specs/for_eid.fam", stringsAsFactors=F)
train_eid <- read.table(paste0("~/athena/doc_score/qc/cv_files/train_eid.", train_frac, ".txt"), stringsAsFactors=F)
test_eid <- read.table(paste0("~/athena/doc_score/qc/cv_files/test_eid.", test_frac, ".txt"), stringsAsFactors=F)
eid <- all_eid[,1]
all_scores <- all_scores[eid %in% train_eid[,1] | eid %in% test_eid[,1],]
eid <- eid[eid %in% train_eid[,1] | eid %in% test_eid[,1]]

#Normalize the scores
best_score <- read.table(paste0("~/athena/doc_score/analyze_score/tune_score/tune_results/", tolower(author), ".best.ss"), stringsAsFactors=F)
best_score <- best_score[best_score[,1] == score_method,2]
scores <- all_scores[,grepl(tolower(author), colnames(all_scores))]
scores <- scores[,apply(scores, 2, function(x) length(unique(x)) > 3)]
scores <- apply(scores, 2, function(x) (x-mean(x)) / (max(abs((x-mean(x)))) * 0.01) )
scores <- scores[,colnames(scores) == best_score,drop=F]

# NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
# NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
#                             ADDED COVARIATES                                      #
# NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
# NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW


#get other scores - scores that I did not compile directly if they are available
#First check to see if there are any additional scores for the disease currently under analysis
#Then if so (in the for loop) we read in the other scores and sort them appropriately
other_scores <- readRDS("~/athena/doc_score/analyze_score/other_scores/final_scores/all_score.1.RDS")
other_defs <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_scores", stringsAsFactors=F, header=T)
other_defs <- other_defs[other_defs[,1] == author,]
if(any(colnames(other_scores) %in% other_defs$PGS_Catalog_name)){
  run_other_scores <- TRUE
  other_scores <- other_scores[,colnames(other_scores) %in% other_defs$PGS_Catalog_name,drop=F]
  other_eid <- read.table("~/athena/doc_score/analyze_score/other_scores/final_scores/eid", stringsAsFactors=F)
  other_scores <- data.frame(eid = other_eid, other_scores)
} else {
  run_other_scores <- FALSE
}

# END NEW #########################################################################################
# END NEW #########################################################################################



#Read in the phenotypes, order is: cancer sr, noncancer sr, icd9, icd10, oper, meds
#selfasses date: 2018-11-22; hesin date: 21/01/2000
pheno <- read.table(paste0("~/athena/doc_score/analyze_score/construct_defs/pheno_defs/diag.", tolower(author), ".txt.gz"), stringsAsFactors=F)
dates <- read.table(paste0("~/athena/doc_score/analyze_score/construct_defs/pheno_defs/time.", tolower(author), ".txt.gz"), stringsAsFactors=F)
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
pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]
scores <- scores[eid %in% pheno_eids, , drop = F]
eid <- eid[eid %in% pheno_eids]
train_eid <- train_eid[train_eid[,1] %in% pheno_eids,]
test_eid <- test_eid[test_eid[,1] %in% pheno_eids,]

dates <- dates[pheno_eids %in% eid]
pheno <- pheno[pheno_eids %in% eid]
pheno_eids <- pheno_eids[pheno_eids %in% eid]
scores <- scores[eid %in% pheno_eids, , drop = F]
eid <- eid[eid %in% pheno_eids]
dates <- dates[order(pheno_eids)[rank(eid)]]
pheno <- pheno[order(pheno_eids)[rank(eid)]]
#eid is the correct order


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


#need to remove people that had diagnosis before accepted start of the study
eid <- eid[surv_df$time > 0]
df <- df[surv_df$time > 0,]
covars <- covars[surv_df$time > 0,]
pheno <- pheno[surv_df$time > 0]
surv_df <- surv_df[surv_df$time > 0,]

# NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
# NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
#                             ADDED COVARIATES                                      #
# NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW
# NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW

#get the other covariates
#Specifically, if there are covariates other than age and sex that other resources believe are risk factors for this disease then we read them in and order them such that they align with the already existing covariates
#We can pull these additional covariates from two different place however, from the ICD or non ICD records
#If from the non-ICD we first get the names from extra covars, and then subset those names out of the single_col_covars
poss_covars <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_covar", stringsAsFactors=F, sep = "\t")
poss_author <- gsub("-", "", author)
extra_covars <- strsplit(poss_covars[poss_covars[,1] == poss_author,2], ",")[[1]]
extra_covars <- gsub(" ", "_", extra_covars)

single_col_covars <- read.table("~/athena/doc_score/analyze_score/get_covars/covar_data/single_col_covars", stringsAsFactors=F, header=T)
single_col_covars <- single_col_covars[,colnames(single_col_covars) %in% c("eid", extra_covars), drop = F]
single_col_covars <- single_col_covars[,!colnames(single_col_covars) %in% c("age", "sex"),drop=F]
if(ncol(single_col_covars) > 1){
  single_col_covars <- single_col_covars[single_col_covars$eid %in% eid,,drop=F]
  single_col_covars <- single_col_covars[order(single_col_covars$eid)[rank(eid)],] #nrow(single_col_covars) may be < length(eid)
  single_col_covars <- single_col_covars[,-1,drop=F]
} else {
  single_col_covars <- NULL
}

#If from the ICD record we read in the specific covariate file made from the ICDs that goes with the disease
#Similar with non-ICD we first have to check if there are any non-ICD covariates relevant, and then if so we go on and read them in and proceed with sorting
#With non-ICD or ICD covariates we leave the covariate file NULL if there is nothing relevant

hesin_decode <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_to_covar_hesin", stringsAsFactors=F)
hesin_decode <- strsplit(hesin_decode[hesin_decode[,1] == poss_author,2], ",")[[1]]
sort_covar <- read.table("~/athena/doc_score/analyze_score/descript_defs/covar_defs_hesin", stringsAsFactors=F)
hesin_decode <- sort_covar[sort_covar[,1] %in% hesin_decode,1]
if(any(grepl(tolower(author), list.files("~/athena/doc_score/analyze_score/get_covars/hesin_covars/")))){
  hesin_covar <- read.table(paste0("~/athena/doc_score/analyze_score/get_covars/hesin_covars/diag.coding.", tolower(author), ".txt.gz"), stringsAsFactors=F, header=F)
  hesin_eid <- read.table(paste0("~/athena/doc_score/analyze_score/get_covars/hesin_covars/eid.txt.gz"), stringsAsFactors=F, header=F)
  colnames(hesin_covar) <- hesin_decode
  hesin_covar <- hesin_covar[hesin_eid[,1] %in% eid,,drop=F]
  hesin_covar <- hesin_covar[order(hesin_eid[,1])[rank(eid)],,drop=F]
  hesin_covar <- hesin_covar[,colSums(hesin_covar) > 0,drop=F]
} else {
  hesin_covar <- NULL
}

#We finish this extra covariate process by combining the ICD and non-ICD files into a singe extra_covar 
#Also we set a variable indicating whether or not we should try to run models with extra covariates
run_extra_covar <- TRUE
if(!is.null(single_col_covars) & !is.null(hesin_covar)){
  extra_covar <- cbind(single_col_covars, hesin_covar)
} else if(!is.null(single_col_covars)){
  extra_covar <- single_col_covars
} else if(!is.null(hesin_covar)){
  extra_covar <- hesin_covar
} else {
  run_extra_covar <- FALSE
}

# END NEW ####################################################################################
##############################################################################################

#add in the other interesting covariates
if(run_other_scores){
  other_scores <- other_scores[other_scores[,1] %in% eid,]
  other_scores <- other_scores[order(other_scores[,1])[rank(eid)],]
  other_scores <- other_scores[,-1,drop=F]
}

if(run_other_scores & run_extra_covar){
  df <- cbind(df, extra_covar, other_scores)
  surv_df <- cbind(surv_df, extra_covar, other_scores)
} else if(run_other_scores){
  df <- cbind(df, other_scores)
  surv_df <- cbind(surv_df, other_scores)
} else if(run_extra_covar){
  df <- cbind(df, extra_covar)
  surv_df <- cbind(surv_df, extra_covar)
}

#Subset the sex
author_defs <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_defs", stringsAsFactors=F, header=T)
subset_sex <- author_defs$sex[author_defs$author == author]
if(subset_sex == "F"){
  eid <- eid[df$sex == 0]
  df <- df[df$sex == 0,]
  surv_df <- surv_df[surv_df$sex == 0,]
  df <- df[,-which(colnames(df) == "sex")]
  surv_df <- surv_df[,-which(colnames(surv_df) == "sex")]
  base_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
}else if(subset_sex == "M"){
  eid <- eid[df$sex == 1]
  df <- df[df$sex == 1,]
  surv_df <- surv_df[surv_df$sex == 1,]
  df <- df[,-which(colnames(df) == "sex")]
  surv_df <- surv_df[,-which(colnames(surv_df) == "sex")]
  base_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
} else {
  base_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
}



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



df_train <- df[eid %in% train_eid[,1],]
df_test <- df[eid %in% test_eid[,1],]
surv_df_train <- surv_df[eid %in% train_eid[,1],]
surv_df_test <- surv_df[eid %in% test_eid[,1],]
fg_diag_train <- fg_diag[fg_diag$eid %in% train_eid[,1],]
fg_diag_test <- fg_diag[fg_diag$eid %in% test_eid[,1],]
fg_death_train <- fg_death[fg_death$eid %in% train_eid[,1],]
fg_death_test <- fg_death[fg_death$eid %in% test_eid[,1],]

df_train$eid <- eid[eid %in% train_eid[,1]]
df_test$eid <- eid[eid %in% test_eid[,1]]
surv_df_train$eid <- eid[eid %in% train_eid[,1]]
surv_df_test$eid <- eid[eid %in% test_eid[,1]]

saveRDS(df_train, paste0("df_output/df_train.", tolower(author), ".RDS"))
saveRDS(df_test, paste0("df_output/df_test.", tolower(author), ".RDS"))
saveRDS(surv_df_train, paste0("df_output/survdf_train.", tolower(author), ".RDS"))
saveRDS(surv_df_test, paste0("df_output/survdf_test.", tolower(author), ".RDS"))



