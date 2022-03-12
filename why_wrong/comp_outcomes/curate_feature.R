library(vroom)
library(data.table)

author <- "mahajan"
#author <- commandArgs(trailingOnly=TRUE)[1]

############################# READ IN AND SET UP #########################################

train_df <- readRDS(paste0("../get_data/df_output/df_train.", author, ".RDS"))
test_df <- readRDS(paste0("../get_data/df_output/df_test.", author, ".RDS"))
surv_train_df <- readRDS(paste0("../get_data/df_output/survdf_train.", author, ".RDS"))
surv_test_df <- readRDS(paste0("../get_data/df_output/survdf_test.", author, ".RDS"))

big_data <- readRDS("../get_data/big_data.RDS")
colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))] <- paste0("X", colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))])

night <- read.csv("~/athena/ukbiobank/setup_sixthPhenos/ukb46898.csv.gz", stringsAsFactors=F, header=T)
eid_night <- night$eid
night <- night[,14:348]
night <- night[,grepl(".0.0", colnames(night))]
eid_night <- eid_night[!is.na(night[,1])]
night <- night[!is.na(night[,1]),]


big_data <- big_data[big_data$eid %in% surv_train_df$eid,]
surv_train_df <- surv_train_df[surv_train_df$eid %in% big_data$eid,]
big_data <- big_data[order(big_data$eid)[rank(surv_train_df$eid)],]


if(mean(train_df$score[train_df$pheno == 1]) < mean(train_df$score[train_df$pheno == 0])){
  surv_train_df$score <- surv_train_df$score * -1
  train_df$score <- train_df$score * -1
}


hi_prs_no_disease_eid <- train_df$eid[train_df$pheno == 0][train_df$score[train_df$pheno == 0] >= sort(train_df$score[train_df$pheno == 0], decreasing=T)[100]]
hi_prs_yes_disease_eid <- train_df$eid[train_df$pheno == 1][train_df$score[train_df$pheno == 1] >= sort(train_df$score[train_df$pheno == 1], decreasing=T)[100]]
lo_prs_no_disease_eid <- train_df$eid[train_df$pheno == 0][train_df$score[train_df$pheno == 0] <= sort(train_df$score[train_df$pheno == 0], decreasing=F)[100]]
lo_prs_yes_disease_eid <- train_df$eid[train_df$pheno == 1][train_df$score[train_df$pheno == 1] <= sort(train_df$score[train_df$pheno == 1], decreasing=F)[100]]



#will have to do with the time, only looking at individuals who were diagnosed at baseline (or near then)
author_feats <- read.csv("author_feats", stringsAsFactors=F, header=T)

if(!is.na(author_feats$code[author_feats$author == author])){
  feat_name <- author_feats$name[author_feats$author == author]
  feat_code <- author_feats$code[author_feats$author == author]
  feat_feature <- author_feats$feature[author_feats$author == author]
  if(feat_code == "top"){
    feat_code <- quantile(big_data[feat_feature][,1], 0.9, na.rm=T)
    comper <- "greater"
  } else {
    comper <- "x"
  }

  date_ass <- read.csv("date_assessed.csv", stringsAsFactors=F, header=T)
  date_ass <- date_ass[date_ass$eid %in% surv_train_df$eid,]
  date_ass <- date_ass[order(date_ass$eid)[rank(surv_train_df$eid)],]
  date_ass$X53.0.0 <- as.Date(date_ass$X53.0.0)

  surv_train_df$date_diag <- as.Date("1999-01-01") + surv_train_df$time

  surv_train_df$use_pheno <- surv_train_df$pheno
  surv_train_df$use_pheno[surv_train_df$date_diag > date_ass$X53.0.0 + 100] <- 0

  if(comper == "greater"){
    two_tabs <- c(table(surv_train_df$use_pheno[big_data[feat_feature][,1] > feat_code & surv_train_df$score > quantile(surv_train_df$score, 0.9) ]),
                  table(surv_train_df$use_pheno[big_data[feat_feature][,1] > feat_code]))
  } else {
    two_tabs <- c(table(surv_train_df$use_pheno[big_data[feat_feature][,1] == feat_code & surv_train_df$score > quantile(surv_train_df$score, 0.9) ]),
                  table(surv_train_df$use_pheno[big_data[feat_feature][,1] == feat_code]))
  }
  two_tabs <- c(feat_name, feat_code, two_tabs)
  write.table(two_tabs, paste0("manual_tabs/", author, ".tabs.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
}





