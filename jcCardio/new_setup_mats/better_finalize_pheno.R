pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]

exome_fam <- read.table("~/athena/ukbiobank/exome/ukbb.exome.22.fam", stringsAsFactors=F)
use_eids <- pheno_eids[pheno_eids %in% exome_fam[,1]]

#############################3
#READ IN DATA
#################################
#read in all of the possible definitions

all_files <- list.files("../get_pheno/pheno_defs/", "time")
pheno_list <- list()
print(length(all_files))
for(i in 1:length(all_files)){
  print(i)
  pheno_list[[i]] <- read.table(paste0("../get_pheno/pheno_defs/", all_files[i]), stringsAsFactors = F)
  pheno_list[[i]] <- pheno_list[[i]][pheno_eids %in% exome_fam[,1],]

  for(j in 1:5){
    if(any(grepl("/", unique(pheno_list[[i]][pheno_list[[i]][,j] != "__________", j])))){
      pheno_list[[i]][,j] <- as.Date(pheno_list[[i]][,j], "%d/%m/%Y")
    } else {

      pheno_list[[i]][,j] <- as.Date(pheno_list[[i]][,j], "%Y-%m-%d")
    }
  }

  temp_df <- data.frame(apply(pheno_list[[i]], 1, function(x) any(!is.na(x))) *1,
                        apply(pheno_list[[i]], 1, function(x) min(x, na.rm=T)), stringsAsFactors = F)
  temp_df[,2] <- as.Date(temp_df[,2])
  colnames(temp_df) <- c("status", "date")
  pheno_list[[i]] <- temp_df

}
names(pheno_list) <- unlist(lapply(strsplit(all_files, ".", fixed=T), function(x) x[2]))


#seperate out the data
#make list for primary phenotypes, comorbidities, and negative outcome
af <- pheno_list[["incident_af"]]
prev_af <- pheno_list[["prev_af"]]
vt <- pheno_list[["va"]]
va <- pheno_list[["va"]]

update_pheno <- function(primary, new){
  for(i in 1:nrow(primary)){
    if(new[i,1] == 1){
      if(primary[i,1] == 0){
        primary[i,1] <- 1
        primary[i,2] <- new[i,2]
      } else if(primary[i,2] > new[i,2]){
        primary[i,2] <- new[i,2]
      }
    }
  }
  return(primary)
}

va <- update_pheno(va, pheno_list[["pvc"]])
va <- update_pheno(va, pheno_list[["cardiac_arrest"]])

primary_list <- list(af, vt, va, pheno_list[["pvc"]], pheno_list[["cardiac_arrest"]])
names(primary_list) <- c("af", "vt", "va", "pvc", "cardiac_arrest")

comor_names <- c("htn", "hld", "dm", "cad_in", "asthma", "copd", "osa", "hypo", "hyper", "ckd", "pvd", "depress")
comor_list <- list()
for(i in 1:length(comor_names)){
  comor_list[[i]] <- pheno_list[[comor_names[i]]]
}
names(comor_list) <- comor_names




#################################################
# Survey Stuff
###############################################
#include race, age, sex, tobacco use, frequent alcohol use

system("zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -d',' -f1,5,6,28,31,41,217,1254,1291,1444 > survey.txt")
survey <- read.table("survey.txt", sep = ",", header = T, stringsAsFactors = F)
colnames(survey) <- c("eid", "sex", "yob", "waist", "hip_ratio", "date_attend", "alcohol", "smoking", "BMI", "tobacco")
survey$hip_ratio <- survey$waist/survey$hip_ratio
survey$alcohol[survey$alcohol == -3] <- NA
survey$smoking[survey$smoking == -3] <- NA
survey$tobacco[survey$tobacco < 0] <- NA
survey$tobacco[survey$tobacco > 0] <- 1
survey$age <- 2021 - survey$yob
survey <- survey[,-which(colnames(survey) == "yob")]
survey$freq_etoh <- 0
survey$freq_etoh[survey$alcohol == 1] <- 1
survey$date_attend <- as.Date(survey$date_attend)
survey <- survey[survey$eid %in% use_eids,]

brit_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/brit_fam", stringsAsFactors=F)
asian_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/asian_fam", stringsAsFactors=F)
african_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/african_fam", stringsAsFactors=F)
euro_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/euro_fam", stringsAsFactors=F)

survey$brit <- 0
survey$asian <- 0
survey$african <- 0
survey$euro <- 0
survey$brit[survey$eid %in% brit_fam[,1]] <- 1
survey$asian[survey$eid %in% asian_fam[,1]] <- 1
survey$african[survey$eid %in% african_fam[,1]] <- 1
survey$euro[survey$eid %in% euro_fam[,1]] <- 1

death_date <- read.table("~/athena/ukbiobank/hesin/death.txt", stringsAsFactors=F, header=T)

survey$all_death <- 0
survey$all_date_death <- NA

for(i in 1:nrow(death_date)){
  survey$all_death[survey$eid == death_date$eid[i]] <- 1
  survey$all_date_death[survey$eid == death_date$eid[i]] <- death_date$date_of_death[i]
}
survey$all_date_death <- as.Date(survey$all_date_death, "%d/%m/%Y")


survey <- survey[order(survey$eid)[rank(use_eids)],]

#avg_time_studied <- as.numeric(mean(AF$date[-which(AF$date < survey$date_attend)] - survey$date_attend[-which(AF$date < survey$date_attend)], na.rm = T))




######################################################
# OUTCOMES
###########################################################

all_death <- data.frame(0, survey$all_date_death)
all_death[!is.na(all_death[,2]), 1] <- 1

outcome_names <- c("cardiomyopathy", "ischemic_stroke", "hf")

remove_after <- function(done_df, OUT){
  temp_bool <- OUT[,2] > done_df[,2]
  temp_bool[is.na(temp_bool)] <- FALSE

  done_df[temp_bool, 1] <- 0
  done_df[temp_bool, 2] <- NA
  return(done_df)
}

outcome_corrected <- list()
for(i in 1:length(primary_list)){
  outcome_corrected[[i]] <- list()

  for(j in 1:length(outcome_names)){
    outcome_corrected[[i]][[j]] <- remove_after(pheno_list[[outcome_names[j]]], primary_list[[i]])
  }

  outcome_corrected[[i]][[(length(outcome_names)+1)]] <- all_death
  names(outcome_corrected[[i]]) <- c(outcome_names, "death")
}
names(outcome_corrected) <- names(primary_list)


######################################################

comor_df <- do.call("cbind", comor_list)
colnames(comor_df) <- paste0("comor_", gsub(".", "_", colnames(comor_df), fixed=T))
survey <- survey[,1:(ncol(survey)-2)]


for(i in 1:length(primary_list)){
  total_outcome <- do.call("cbind", outcome_corrected[[i]])
  
  outcome_date <- apply(total_outcome[,grepl("date", colnames(total_outcome))], 1, min, na.rm=T)
  outcome_bool <- rep(0, length(outcome_date))
  outcome_bool[!is.na(outcome_date)] <- 1

  outcome_date_nocm <- apply(total_outcome[,grepl("date", colnames(total_outcome)) & !grepl("cardiomyopathy", colnames(total_outcome))], 1, min, na.rm=T)
  outcome_bool_nocm <- rep(0, length(outcome_date_nocm))
  outcome_bool_nocm[!is.na(outcome_date_nocm)] <- 1

  colnames(total_outcome) <- paste0("outcome_", gsub(".", "_", colnames(total_outcome), fixed=T))
  total_outcome$outcome_composite_date <- outcome_date
  total_outcome$outcome_composite_status <- outcome_bool
  total_outcome$outcome_composite_nocm_date <- outcome_date_nocm
  total_outcome$outcome_composite_nocm_status <- outcome_bool_nocm

  primary <- primary_list[[i]]
  colnames(primary) <- paste0("primary_", names(primary_list)[i], "_", colnames(primary))
  if(names(primary_list)[[i]] == "af"){
    primary$prev_af_status <- prev_af$status
    primary$prev_af_date <- prev_af$date
  }

  total_df <- cbind(survey, comor_df, total_outcome, primary)
  saveRDS(total_df, paste0(names(primary_list)[i], "_pheno_mat.RDS"))
}


