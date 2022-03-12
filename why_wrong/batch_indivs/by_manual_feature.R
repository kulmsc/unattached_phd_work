author <- "bentham"

#system(paste0("Rscript get_big_data.R ", author))
big_data <- readRDS(paste0("use_big_data.train.", author, ".RDS"))
train_df <- readRDS(paste0("use_train_df.", author, ".RDS"))
test_data <- readRDS(paste0("use_big_data.test.", author, ".RDS"))
valid_df <- readRDS(paste0("use_test_df.", author, ".RDS"))

#34-0.0 - year of birth
train_age_1 <- train_df$eid[big_data$SURVEY_X34.0.0 < quantile(big_data$SURVEY_X34.0.0, 0.05)]
train_age_2 <- train_df$eid[big_data$SURVEY_X34.0.0 < quantile(big_data$SURVEY_X34.0.0, 0.01)]

test_age_1 <- valid_df$eid[test_data$SURVEY_X34.0.0 < quantile(test_data$SURVEY_X34.0.0, 0.05)]
test_age_2 <- valid_df$eid[test_data$SURVEY_X34.0.0 < quantile(test_data$SURVEY_X34.0.0, 0.01)]


#2178-0.0 - overall health rating
train_health_1 <- train_df$eid[big_data$SURVEY_X2178.0.0 == 0]

test_health_1 <- valid_df$eid[test_data$SURVEY_X2178.0.0 == 0]


#Count ICD and OPCS
train_num_icd <- rowSums(big_data[,grep("ICD", colnames(big_data))])
train_num_ehr <- rowSums(big_data[,c(grep("ICD", colnames(big_data), grep("OPCS", colnames(big_data))))])

test_num_icd <- rowSums(test_data[,grep("ICD", colnames(test_data))])
test_num_ehr <- rowSums(test_data[,c(grep("ICD", colnames(test_data), grep("OPCS", colnames(test_data))))])

train_icd_1 <- train_df$eid[train_num_icd == 0]
test_icd_1 <- valid_df$eid[test_num_icd == 0]

train_ehr_1 <- train_df$eid[train_num_ehr == 0]
test_ehr_1 <- valid_df$eid[test_num_ehr == 0]


####################################################################

train_list <- list(train_age_1, train_age_2, train_health_1, train_icd_1, train_ehr_1)
test_list <- list(test_age_1, test_age_2, test_health_1, test_icd_1, test_ehr_1)

saveRDS(train_list, paste0("out_batch/manual_stuff.train.", author, ".RDS"))
saveRDS(test_list, paste0("out_batch/manual_stuff.test.", author, ".RDS"))
