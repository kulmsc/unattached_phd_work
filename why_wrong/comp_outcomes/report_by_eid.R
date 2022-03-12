library(vroom)

#author <- "jc"
author <- commandArgs(trailingOnly=TRUE)[1]

############################# READ IN AND SET UP #########################################

train_df <- readRDS(paste0("../get_data/df_output/df_train.", author, ".RDS"))
test_df <- readRDS(paste0("../get_data/df_output/df_test.", author, ".RDS"))
surv_train_df <- readRDS(paste0("../get_data/df_output/survdf_train.", author, ".RDS"))
surv_test_df <- readRDS(paste0("../get_data/df_output/survdf_test.", author, ".RDS"))

big_data <- readRDS("../get_data/big_data.RDS")
colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))] <- paste0("X", colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))])

big_data <- big_data[big_data$eid %in% test_df$eid,]
test_df <- test_df[test_df$eid %in% big_data$eid,]
big_data <- big_data[order(big_data$eid)[rank(test_df$eid)],]


if(mean(test_df$score[test_df$pheno == 1]) < mean(test_df$score[test_df$pheno == 0])){
  test_df$score <- test_df$score * -1
}

test_df <- test_df[order(test_df$score, decreasing = T),]
hi_prs_no_disease_eid <- test_df$eid[test_df$pheno == 0][test_df$score[test_df$pheno == 0] >= sort(test_df$score[test_df$pheno == 0], decreasing=T)[100]]
hi_prs_yes_disease_eid <- test_df$eid[test_df$pheno == 1][test_df$score[test_df$pheno == 1] >= sort(test_df$score[test_df$pheno == 1], decreasing=T)[100]]
test_df <- test_df[order(test_df$score, decreasing = F),]
lo_prs_no_disease_eid <- test_df$eid[test_df$pheno == 0][test_df$score[test_df$pheno == 0] <= sort(test_df$score[test_df$pheno == 0], decreasing=F)[100]]
lo_prs_yes_disease_eid <- test_df$eid[test_df$pheno == 1][test_df$score[test_df$pheno == 1] <= sort(test_df$score[test_df$pheno == 1], decreasing=F)[100]]

###########################################################################
make_report <- function(curr_eid, where_from, rank_eid){
#curr_eid = hi_prs_no_disease_eid[1] 
#where_from = "hi_prs_no_disease_eid"

base_info <- c()
#What do I want in my report
#age, sex, bmi, height, weight, BP
base_info["age"] <- 2021 - big_data["SURVEY_X34.0.0"][big_data$eid == curr_eid,1]
sex <- big_data["SURVEY_X31.0.0"][big_data$eid == curr_eid,1]
if(sex == 0){
  base_info["sex"] <- "Female"
} else {
  base_info["sex"] <- "Male"
}
base_info["weight"] <- big_data["SURVEY_X21002.0.0"][big_data$eid == curr_eid,1]
base_info["height"] <- big_data["SURVEY_X50.0.0"][big_data$eid == curr_eid,1]
base_info["bmi"] <- round(as.numeric(big_data["SURVEY_X21002.0.0"][big_data$eid == curr_eid,1])/((as.numeric(big_data["SURVEY_X50.0.0"][big_data$eid == curr_eid,1])/100)^2), 4)
base_info["bp_diastolic"] <- big_data["SURVEY_X4079.0.0"][big_data$eid == curr_eid,1]
base_info["bp_systolic"] <- big_data["SURVEY_X4080.0.0"][big_data$eid == curr_eid,1]


#alcohol, smoking
addict_pos <- c("No", "Yes")
base_info["ever_addicted"] <- addict_pos[as.numeric(big_data["SURVEY_X20401.0.0"][big_data$eid == curr_eid,1]) + 1]
alc_pos <- c("daily", "multi-weekly", "weekly", "monthy", "occasion", "never")
base_info["alc_val"] <- alc_pos[which(big_data[1,grep("1558.0", colnames(big_data))] == 1)]
smoke_pos <- c("never", "previous", "current")
base_info["smoke_val"] <- smoke_pos[which(big_data[1,grep("20116.0", colnames(big_data))] == 1)]
base_info["pack_years"] <- big_data["SURVEY_X20161.0.0"][big_data$eid == curr_eid,1]

get_one_hot <- function(type_name, coding_file){
  all_data <- big_data[big_data$eid == curr_eid, grep(type_name, colnames(big_data))]
  data_codes <- names(all_data)[all_data != 0]
  data_codes <- unlist(lapply(strsplit(data_codes, "_"), function(x) x[2]))
  if(!is.null(data_codes)){
    coding_dict <- as.data.frame(vroom(paste0("~/athena/prep_uk_data/", coding_file), delim="\t"))
    data_vals <- coding_dict[coding_dict[,1] %in% data_codes,2]
    if(length(data_vals) == 0){
      data_vals <- "None"
    }
  } else {
    data_vals <- "None"
  }
  return(data_vals)
}

#ICD codes
icd10_vals <- get_one_hot("ICD10", "coding19.tsv")

icd9_vals <- get_one_hot("ICD9", "coding87.tsv")

noncancer_vals <- get_one_hot("NONCANCER", "coding6.tsv")

cancer_vals <- get_one_hot("CANCER", "coding3.tsv")

meds_vals <- get_one_hot("MEDS", "coding4.tsv")

opcs_vals <- get_one_hot("OPCS", "coding240.tsv")

#other survey answers (exercise, diet, happiness, job)
other_data <- c()
poss_happiness <- c("extremely happy", "very happy", "moderately happy", "moderately unhappy", "very unhappy", "exteremely unhappy")
other_data["general_happiness"] <- poss_happiness[big_data["SURVEY_X20458.0.0"][big_data$eid == curr_eid,1]][1]
other_data["happiness_own_health"] <- poss_happiness[big_data["SURVEY_X20459.0.0"][big_data$eid == curr_eid,1]][1]
other_data["sleep_duration"] <- paste0(big_data["SURVEY_X1160.0.0"][big_data$eid == curr_eid,1], " hours")
poss_income <- c("less than 18k", "18k to 31k", "31k to 52k", "52k to 100k", "> 100k")
other_data["income_before_tax"] <- poss_income[big_data["SURVEY_X738.0.0"][big_data$eid == curr_eid,1]][1]
other_data["age_complete_full_time"] <- paste0(big_data["SURVEY_X845.0.0"][big_data$eid == curr_eid,1], " years")
other_data["time_at_current_address"] <- paste0(big_data["SURVEY_X699.0.0"][big_data$eid == curr_eid,1], " years")
poss_sun_uv <- c("never/rarely", "sometimes", "most of the time", "always", "do not go in sun")
other_data["use_sun_uv_protection"] <- poss_sun_uv[big_data["SURVEY_X2267.0.0"][big_data$eid == curr_eid,1]][1]
other_data["fresh_fruit_intake"] <- paste0(big_data["SURVEY_X1309.0.0"][big_data$eid == curr_eid,1], " servings")
other_data["water_intake"] <- paste0(big_data["SURVEY_X1528.0.0"][big_data$eid == curr_eid,1], " glasses")
other_data["duration_of_walks"] <- paste0(big_data["SURVEY_X874.0.0"][big_data$eid == curr_eid,1], " minutes")
poss_pace <- c("slow", "average", "brisk")
other_data["walking_pace"] <- poss_pace[big_data["SURVEY_X924.0.0"][big_data$eid == curr_eid,1]][1]
poss_visit <- c("almost daily", "2-4 per week", "1 per week", "1 per month", "1 few months", "never", "none")
other_data["freq_friends_visit"] <- poss_visit[big_data["SURVEY_X1031.0.0"][big_data$eid == curr_eid,1]][1]
other_data["time_outdoors_summer"] <- paste0(big_data["SURVEY_X1050.0.0"][big_data$eid == curr_eid,1], " hours")
other_data["time_outdoors_winter"] <- paste0(big_data["SURVEY_X1060.0.0"][big_data$eid == curr_eid,1], " hours")
other_data["time_watch_tv"] <- paste0(big_data["SURVEY_X1070.0.0"][big_data$eid == curr_eid,1], " hours")
other_data["time_use_computer"] <- paste0(big_data["SURVEY_X1080.0.0"][big_data$eid == curr_eid,1], " hours")
other_data["time_driving"] <- paste0(big_data["SURVEY_X1090.0.0"][big_data$eid == curr_eid,1], " hours")
poss_test <- c("No", "Yes")
other_data["ever_had_psa"] <- poss_test[big_data["SURVEY_X2365.0.0"][big_data$eid == curr_eid,1] + 1]
other_data["ever_had_mammogram"] <- poss_test[big_data["SURVEY_X2674.0.0"][big_data$eid == curr_eid,1] + 1]


#Jobs
all_jobs <- readRDS("../get_data/all_jobs.RDS")
curr_jobs <- all_jobs[all_jobs$eid == curr_eid,-1]
curr_jobs <- curr_jobs[!is.na(curr_jobs)]

job_decoder <- as.data.frame(vroom("~/athena/prep_uk_data/coding2.tsv", delim = "\t"))
all_jobs <- job_decoder$meaning[job_decoder$coding %in% curr_jobs]
if(length(all_jobs) == 0){
  all_jobs <- "None"
}


#Make final report ###################################
max_val = max(unlist(lapply(list(icd10_vals, icd9_vals, noncancer_vals, cancer_vals, meds_vals, opcs_vals), length)))

report <- data.frame(matrix("", nrow = 37, ncol = max_val + 1), stringsAsFactors = F)
report[1:11,2] <- base_info
report[1:11,1] <- names(base_info)
report[12, 2:(length(icd10_vals)+1)] <- icd10_vals
report[13, 2:(length(icd9_vals)+1)] <- icd9_vals
report[14, 2:(length(noncancer_vals)+1)] <- noncancer_vals
report[15, 2:(length(cancer_vals)+1)] <- cancer_vals
report[16, 2:(length(meds_vals)+1)] <- meds_vals
report[17, 2:(length(opcs_vals)+1)] <- opcs_vals
report[12:17,1] <- c("icd10", "icd9", "noncancer", "cancer", "meds", "opcs")
report[18:36,2] <- other_data
report[18:36,1] <- names(other_data)
report[37,2:(length(all_jobs)+1)] <- all_jobs
report[37,1] <- "jobs"

write.table(report, paste0("report_res/", author, ".", where_from, ".", curr_eid, ".", rank_eid, ".txt"), row.names = F, col.names = F, quote = F, sep = "\t")
}


#use_eid <- readRDS("~/athena/jcCardio/setup_mats/for_why_not.af.RDS")
#for(eid in use_eid[[1]]){
#  make_report(eid, "healthy_lpv", 1)
#}

make_report(hi_prs_no_disease_eid[1], "hi_prs_no_disease_eid", 1)
make_report(hi_prs_yes_disease_eid[1], "hi_prs_yes_disease_eid", 1)
make_report(lo_prs_no_disease_eid[1], "lo_prs_no_disease_eid", 1)
make_report(lo_prs_yes_disease_eid[1], "lo_prs_yes_disease_eid", 1)

make_report(hi_prs_no_disease_eid[2], "hi_prs_no_disease_eid", 1)
make_report(hi_prs_yes_disease_eid[2], "hi_prs_yes_disease_eid", 1)
make_report(lo_prs_no_disease_eid[2], "lo_prs_no_disease_eid", 1)
make_report(lo_prs_yes_disease_eid[2], "lo_prs_yes_disease_eid", 1)

make_report(hi_prs_no_disease_eid[3], "hi_prs_no_disease_eid", 1)
make_report(hi_prs_yes_disease_eid[3], "hi_prs_yes_disease_eid", 1)
make_report(lo_prs_no_disease_eid[3], "lo_prs_no_disease_eid", 1)
make_report(lo_prs_yes_disease_eid[3], "lo_prs_yes_disease_eid", 1)
