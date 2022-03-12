
#author <- "bentham"
#author <- commandArgs(trailingOnly=TRUE)[1]

############################# READ IN AND SET UP #########################################
author <- "jc.pvc"
comp <- "pvc"

use_eid <- readRDS(paste0("~/athena/jcCardio/setup_mats/for_why_not.all.", comp, ".RDS"))
ppl_1 <- use_eid[[1]]
ppl_2 <- use_eid[[2]]

#train_df <- readRDS(paste0("../get_data/df_output/df_train.", author, ".RDS"))
#test_df <- readRDS(paste0("../get_data/df_output/df_test.", author, ".RDS"))
#surv_train_df <- readRDS(paste0("../get_data/df_output/survdf_train.", author, ".RDS"))
#surv_test_df <- readRDS(paste0("../get_data/df_output/survdf_test.", author, ".RDS"))

big_data <- readRDS("../get_data/big_data.RDS")
colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))] <- paste0("X", colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))])

#big_data <- big_data[big_data$eid %in% train_df$eid,]
#train_df <- train_df[train_df$eid %in% big_data$eid,]
#big_data <- big_data[order(big_data$eid)[rank(train_df$eid)],]


#if(mean(train_df$score[train_df$pheno == 1]) < mean(train_df$score[train_df$pheno == 0])){
#  train_df$score <- train_df$score * -1
#}

#hi_prs_no_disease_eid <- train_df$eid[train_df$pheno == 0][train_df$score[train_df$pheno == 0] >= sort(train_df$score[train_df$pheno == 0], decreasing=T)[100]]
#hi_prs_yes_disease_eid <- train_df$eid[train_df$pheno == 1][train_df$score[train_df$pheno == 1] >= sort(train_df$score[train_df$pheno == 1], decreasing=T)[100]]
#lo_prs_no_disease_eid <- train_df$eid[train_df$pheno == 0][train_df$score[train_df$pheno == 0] <= sort(train_df$score[train_df$pheno == 0], decreasing=F)[100]]
#lo_prs_yes_disease_eid <- train_df$eid[train_df$pheno == 1][train_df$score[train_df$pheno == 1] <= sort(train_df$score[train_df$pheno == 1], decreasing=F)[100]]

###########################################################

get_time <- function(eid_list, save_name){

  data_collector <- data.frame(matrix(0, nrow = length(eid_list), ncol = 10))
  colnames(data_collector) <- c("last_icd_diag", "last_oper_diag", "last_diag", "total_visit", "total_icd", "total_oper", "is_dead", "death_date", "is_dropped", "age")


  write.table(eid_list, "search_eid.txt", row.names = F, col.names = F)
  system("./get_hesin.sh")

  hesin <- read.table("small_hesin", stringsAsFactors = F, sep = "\t")
  hesin_diag <- read.table("small_hesin_diag", stringsAsFactors = F, sep = "\t")
  hesin_oper <- read.table("small_hesin_oper", stringsAsFactors = F, sep = "\t")

  hesin[,5] <- as.Date(hesin[,5], "%d/%m/%Y")

  for(i in 1:length(eid_list)){
    curr_eid <- eid_list[i]
    #time of last diagnosis
    data_collector$last_icd_diag[i] <- as.character(max(hesin[hesin[,1] == curr_eid & hesin[,2] %in% hesin_diag[hesin_diag[,1] == curr_eid, 2], 5], na.rm = T))
    data_collector$last_oper_diag[i] <- as.character(max(hesin[hesin[,1] == curr_eid & hesin[,2] %in% hesin_diag[hesin_oper[,1] == curr_eid, 2], 5], na.rm = T))
    data_collector$last_diag[i] <- as.character(max(hesin[hesin[,1] == curr_eid, 5], na.rm = T))

    #get total doctors visits
    data_collector$total_visit[i] <- sum(hesin[,1] == curr_eid)
    data_collector$total_icd[i] <- sum(hesin_diag[,1] == curr_eid)
    data_collector$total_oper[i] <- sum(hesin_oper[,1] == curr_eid)

    #get death date?
    death <- read.table("small_death", stringsAsFactors=F)
    death[,5] <- as.Date(death[,5], "%d/%m/%Y")
    data_collector$is_dead[i] <- curr_eid %in% death[,1]
    if(data_collector$is_dead[i]){
      data_collector$death_date[i] <- as.character(death[death[,1] %in% curr_eid,5])
    } else {
      data_collector$death_date[i] <- NA
    }

    drop_outs <- read.table("~/athena/ukbiobank/hesin/drop_out.txt", stringsAsFactors=F)
    data_collector$is_dropped[i] <- curr_eid %in% drop_outs[,1]


    #get age
    data_collector$age[i] <- 2021 - big_data["SURVEY_X34.0.0"][big_data$eid == curr_eid,1]
  }

  saveRDS(data_collector, paste0("timeline_res/", author, ".", save_name, ".time.RDS"))
}

get_time(ppl_1, paste0(author, ".healthy"))
get_time(ppl_2, paste0(author, ".disease"))
#get_time(hi_prs_yes_disease_eid, "hi_prs_yes_disease_eid")
#get_time(lo_prs_yes_disease_eid, "lo_prs_yes_disease_eid")
#get a sense if they dropped off the map
#not sure how to empirically measure this
#can instead simply plot a graph of the timelines 
