
author <- "bentham"
#author <- commandArgs(trailingOnly=TRUE)[1]

############################# READ IN AND SET UP #########################################

#train_df <- readRDS(paste0("../get_data/df_output/df_train.", author, ".RDS"))
#test_df <- readRDS(paste0("../get_data/df_output/df_test.", author, ".RDS"))
#surv_train_df <- readRDS(paste0("../get_data/df_output/survdf_train.", author, ".RDS"))
#surv_test_df <- readRDS(paste0("../get_data/df_output/survdf_test.", author, ".RDS"))

big_data <- readRDS("../get_data/big_data.RDS")
colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))] <- paste0("X", colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))])

#train_big_data <- big_data[big_data$eid %in% train_df$eid,]
#train_df <- train_df[train_df$eid %in% train_big_data$eid,]
#train_big_data <- train_big_data[order(train_big_data$eid)[rank(train_df$eid)],]

#test_big_data <- big_data[big_data$eid %in% test_df$eid,]
#test_df <- test_df[test_df$eid %in% test_big_data$eid,]
#test_big_data <- test_big_data[order(test_big_data$eid)[rank(test_df$eid)],]

#if(mean(train_df$score[train_df$pheno == 1]) < mean(train_df$score[train_df$pheno == 0])){
#  train_df$score <- train_df$score * -1
#  test_df$score <- test_df$score * -1
#}

#yes_disease_eid <- train_df$eid[train_df$pheno == 0]
#no_disease_eid <- train_df$eid[train_df$pheno == 1]


###########################################################
drop_outs <- read.table("~/athena/ukbiobank/hesin/drop_out.txt", stringsAsFactors=F)

get_time <- function(eid_list, save_name, big_data){

  data_collector <- data.frame(matrix(0, nrow = length(eid_list), ncol = 18))
  colnames(data_collector) <- c("last_icd_diag", "last_oper_diag", "last_diag", "first_icd_diag", "first_oper_diag", "first_diag", "total_time", "long_vists", "num_eme_admit", "num_psy_admit", "total_visit", "total_icd", "total_oper", "is_dead", "death_date", "age_at_death", "is_dropped", "age")

  data_collector$death_date <- NA
  data_collector$age_at_death <- NA
  #epidur, total time in hospital
  #admisorc_uni
  #classpat_uni

  write.table(eid_list, "search_eid.txt", row.names = F, col.names = F)
  system("./get_hesin.sh")

  hesin <- read.table("small_hesin", stringsAsFactors = F, sep = "\t")
  hesin$admi <- substr(hesin[,23], 1, 2)
  hesin_diag <- read.table("small_hesin_diag", stringsAsFactors = F, sep = "\t")
  hesin_oper <- read.table("small_hesin_oper", stringsAsFactors = F, sep = "\t")
  death <- read.table("small_death", stringsAsFactors=F)

  hesin[,5] <- as.Date(hesin[,5], "%d/%m/%Y")
  death[,5] <- as.Date(death[,5], "%d/%m/%Y")
  death$char_date <- as.character(death[,5])

  for(i in 1:length(eid_list)){
    curr_eid <- eid_list[i]
    print(i)

    #time of last diagnosis
    data_collector$last_icd_diag[i] <- max(hesin[hesin[,1] == curr_eid & hesin[,2] %in% hesin_diag[hesin_diag[,1] == curr_eid, 2], 5], na.rm = T)
    data_collector$last_oper_diag[i] <- max(hesin[hesin[,1] == curr_eid & hesin[,2] %in% hesin_diag[hesin_oper[,1] == curr_eid, 2], 5], na.rm = T)
    data_collector$last_diag[i] <- max(hesin[hesin[,1] == curr_eid, 5], na.rm = T)

    #time of last diagnosis
    data_collector$first_icd_diag[i] <- min(hesin[hesin[,1] == curr_eid & hesin[,2] %in% hesin_diag[hesin_diag[,1] == curr_eid, 2], 5], na.rm = T)
    data_collector$first_oper_diag[i] <- min(hesin[hesin[,1] == curr_eid & hesin[,2] %in% hesin_diag[hesin_oper[,1] == curr_eid, 2], 5], na.rm = T)
    data_collector$first_diag[i] <- min(hesin[hesin[,1] == curr_eid, 5], na.rm = T)

    #data_collector$total_time[i] <- sum(hesin[hesin[,1] == curr_eid,7], na.rm=T) #total time
    #data_collector$long_visits[i] <- sum(hesin[hesin[,1] == curr_eid,7] > 3, na.rm=T)
    #data_collector$num_eme_admit[i] <- sum(hesin[hesin[,1] == curr_eid,43] == 20) #emergency admissions
    #data_collector$num_psy_admit[i] <- sum(hesin[hesin[,1] == curr_eid,43] == 40) #psy admissions

    #get total doctors visits
    #data_collector$total_visit[i] <- sum(hesin[,1] == curr_eid)
    data_collector$total_icd[i] <- sum(hesin_diag[,1] == curr_eid)
    data_collector$total_oper[i] <- sum(hesin_oper[,1] == curr_eid)

    #get death date?
    #data_collector$is_dead[i] <- curr_eid %in% death[,1]
    if(data_collector$is_dead[i]){
      data_collector$death_date[i] <- death[death[,1] %in% curr_eid,5]
      data_collector$age_at_death[i] <- as.numeric(strsplit(death[death[,1] %in% curr_eid,6], "-")[[1]][1]) - big_data["SURVEY_X34.0.0"][big_data$eid == curr_eid,1]
    }

    #data_collector$is_dropped[i] <- curr_eid %in% drop_outs[,1]

    #get age
    #data_collector$age[i] <- 2021 - big_data["SURVEY_X34.0.0"][big_data$eid == curr_eid,1]
  }

  split_hesin <- split(hesin, hesin$V1)

  total_time <- unlist(lapply(split_hesin, function(x) sum(x[,7], na.rm=T)))
  long_visits <- unlist(lapply(split_hesin, function(x) sum(x[,7] > 3, na.rm=T)))
  num_eme_admit<- unlist(lapply(split_hesin, function(x) sum(x[,43] == 20, na.rm=T)))
  num_psy_admit <- unlist(lapply(split_hesin, function(x) sum(x[,43] == 40, na.rm=T)))
  total_visit <- unlist(lapply(split_hesin, nrow))

  insert_df <- data.frame("eid" = names(split_hesin), total_time, long_visits, num_eme_admit, num_psy_admit, total_visit)
  temp_df <- data.frame(matrix(0, nrow = nrow(data_collector) - nrow(insert_df), ncol = ncol(insert_df)))
  colnames(temp_df) <- colnames(insert_df)
  temp_df$eid <- eid_list[!(eid_list %in% insert_df$eid)]
  insert_df <- rbind(temp_df, insert_df)
  insert_df$eid <- as.numeric(insert_df$eid)
  insert_df <- insert_df[order(insert_df$eid)[rank(eid_list)],]
  
  data_collector$total_time <- insert_df$total_time
  data_collector$long_visits <- insert_df$long_visits
  data_collector$num_eme_admit <- insert_df$num_eme_admit
  data_collector$num_psy_admit <- insert_df$num_psy_admit
  data_collector$total_visit <- insert_df$total_visit

  data_collector$is_dropped[eid_list %in% drop_outs[,1]] <- 1
  data_collector$age <- 2021 - big_data["SURVEY_X34.0.0"][,1] 
  data_collector$is_dead[eid_list %in% death[,1]] <- 1

  saveRDS(data_collector, paste0("timeline_res/", author, ".", save_name, ".time.RDS"))
}

#get_time(no_disease_eid, "all_no_disease_eid", train_big_data)
#get_time(yes_disease_eid, "all_yes_disease_eid", train_big_data)
#get_time(test_df$eid, "test_eid", test_big_data)
#get_time(train_df$eid, "train_eid", train_big_data)
#get_time(big_data$eid, "new_all_eid", big_data)

#get a sense if they dropped off the map
#not sure how to empirically measure this
#can instead simply plot a graph of the timelines 
