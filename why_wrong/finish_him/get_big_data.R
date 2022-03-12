library(vroom)
library(data.table)

#author <- "christophersen"
#tort <- "train"
args <- commandArgs(trailingOnly=TRUE)
author <- args[1]
tort <- args[2]

############################# READ IN AND SET UP #########################################

train_df <- readRDS(paste0("../get_data/df_output/df_train.", author, ".RDS"))
test_df <- readRDS(paste0("../get_data/df_output/df_test.", author, ".RDS"))
surv_train_df <- readRDS(paste0("../get_data/df_output/survdf_train.", author, ".RDS"))
surv_test_df <- readRDS(paste0("../get_data/df_output/survdf_test.", author, ".RDS"))

if(tort == "train"){
  print("skip")
} else {
  train_df <- test_df
}

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
minor_counts <- readRDS("../comp_outcomes/minor_counts.RDS")
big_data <- big_data[,!(colnames(big_data) %in% names(minor_counts)[minor_counts < 250])]
#eid <- train_df$eid

#add extra data
#poss_covars <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_covar", stringsAsFactors=F, sep = "\t")
#poss_covars[,1] <- tolower(poss_covars[,1])
#poss_author <- gsub("-", "", author)
#extra_covars <- strsplit(poss_covars[poss_covars[,1] == poss_author,2], ",")[[1]]
#extra_covars <- gsub(" ", "_", extra_covars)

#single_col_covars <- read.table("~/athena/doc_score/analyze_score/get_covars/covar_data/single_col_covars", stringsAsFactors=F, header=T)
#single_col_covars <- single_col_covars[,colnames(single_col_covars) %in% c("eid", extra_covars), drop = F]
#single_col_covars <- single_col_covars[,!colnames(single_col_covars) %in% c("age", "sex"),drop=F]
#if(ncol(single_col_covars) > 1){
#  single_col_covars <- single_col_covars[single_col_covars$eid %in% eid,,drop=F]
#  single_col_covars <- single_col_covars[order(single_col_covars$eid)[rank(eid)],] #nrow(single_col_covars) may be < length(eid)
#  single_col_covars <- single_col_covars[,-1,drop=F]
#} else {
#  single_col_covars <- NULL
#}

#If from the ICD record we read in the specific covariate file made from the ICDs that goes with the disease
#Similar with non-ICD we first have to check if there are any non-ICD covariates relevant, and then if so we go on and read them in and proceed with sorting
#With non-ICD or ICD covariates we leave the covariate file NULL if there is nothing relevant

#hesin_decode <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_to_covar_hesin", stringsAsFactors=F)
#hesin_decode[,1] <- tolower(hesin_decode[,1])
#hesin_decode <- strsplit(hesin_decode[hesin_decode[,1] == poss_author,2], ",")[[1]]
#sort_covar <- read.table("~/athena/doc_score/analyze_score/descript_defs/covar_defs_hesin", stringsAsFactors=F)
#hesin_decode <- sort_covar[sort_covar[,1] %in% hesin_decode,1]
#if(any(grepl(tolower(author), list.files("~/athena/doc_score/analyze_score/get_covars/hesin_covars/")))){
#  hesin_covar <- read.table(paste0("~/athena/doc_score/analyze_score/get_covars/hesin_covars/diag.coding.", tolower(author), ".txt.gz"), stringsAsFactors=F, header=F)
#  hesin_eid <- read.table(paste0("~/athena/doc_score/analyze_score/get_covars/hesin_covars/eid.txt.gz"), stringsAsFactors=F, header=F)
#  colnames(hesin_covar) <- hesin_decode
#  hesin_covar <- hesin_covar[hesin_eid[,1] %in% eid,,drop=F]
#  hesin_covar <- hesin_covar[order(hesin_eid[,1])[rank(eid)],,drop=F]
#  hesin_covar <- hesin_covar[,colSums(hesin_covar) > 0,drop=F]
#} else {
#  hesin_covar <- NULL
#}

#if(!is.null(single_col_covars) & !is.null(hesin_covar)){
#  extra_covar <- cbind(single_col_covars, hesin_covar)
#} else if(!is.null(single_col_covars)){
#  extra_covar <- single_col_covars
#} else if(!is.null(hesin_covar)){
#  extra_covar <- hesin_covar
#} else {
#  run_extra_covar <- FALSE
#}





saveRDS(big_data, paste0("use_big_data.", tort, ".RDS"))
saveRDS(train_df, paste0("use_", tort, "_df.RDS"))




