library(vroom)

################################################################
#                    GET PHENO
###############################################################

REAL_age_ind = 5

pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]



#######################################################
#                          GET COVAR
######################################################

csv_file <- as.data.frame(vroom("~/athena/ukbiobank/phenotypes/ukb26867.csv.gz", delim = ","))
brit_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/brit_fam", stringsAsFactors=F)
alive_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/qc_alive_fam", stringsAsFactors=F)

pcs_ind <- 1310:1319
age_ind <- 6
sex_ind <- 5
bmi_ind <- 1291
smoke_ind <- 1265 #NA to 0
alcohol_ind <- 217

covar_file <- csv_file[,c(1, 1, pcs_ind, age_ind, sex_ind, bmi_ind, smoke_ind, alcohol_ind)]
covar_file[,13] <- covar_file[,13] + csv_file[,40]/12
covar_file[,13] <- 2010 - covar_file[,13]
covar_file[is.na(covar_file[,16]), 16] <- 0 #convert smoke 0 to NA
covar_file[covar_file[,17] == 3 & !is.na(covar_file[,17]), 17] <- NA
covar_file[,17] <- abs(-1*covar_file[,17])
covar_file[is.na(covar_file[,17]), 17] <- mean(covar_file[,17], na.rm=T)

rm(csv_file)

covar_file <- covar_file[covar_file[,1] %in% brit_fam[,1] & covar_file[,1] %in% alive_fam[,1],]
simple_covar_file <- data.frame(covar_file[,1:16])

simple_covar_file[is.na(simple_covar_file[,ncol(simple_covar_file)]), ncol(simple_covar_file)] <- 0 #NA pack year smoking to 0
#simple_covar_file[is.na(simple_covar_file[,15]),15] <- mean(simple_covar_file[,15], na.rm = TRUE)
simple_covar_file <- simple_covar_file[complete.cases(simple_covar_file),]

covar_file <- covar_file[complete.cases(covar_file),]

for(i in 3:ncol(covar_file)){
  if(max(covar_file[,i]) == 1 & min(covar_file[,i]) == 0){
    print("skip")
  } else {
    covar_file[,i] <- (covar_file[,i] - min(covar_file[,i]))/(max(covar_file[,i]) - min(covar_file[,i]))
  }
}

for(i in 3:ncol(simple_covar_file)){
  if(max(simple_covar_file[,i]) == 1 & min(simple_covar_file[,i]) == 0){
    print("skip")
  } else {
    simple_covar_file[,i] <- (simple_covar_file[,i] - min(simple_covar_file[,i]))/(max(simple_covar_file[,i]) - min(simple_covar_file[,i]))
  }
}

############################################################
#                           FINISH
###########################################################

timed_files <- c("kr_covar.knee_injury.eid")
for(f in timed_files){
  spec_covar <- read.table(paste0("../get_timed_covars/", f), stringsAsFactors=F)
  var_name <- strsplit(f, ".", fixed = T)[[1]][2]
  covar_file[var_name] <- 0
  covar_file[var_name][covar_file$eid %in% spec_covar[,1],1] <- 1
}





#use_eids <- read.table("../comboortho/ptha_all_age.eids", stringsAsFactors=F)
all_use_eids <- list()
for(kk in 1:5){
  all_use_eids[[kk]] <- read.table(paste0("eid_gwas_", kk), stringsAsFactors=F)
}

pheno_df <- data.frame("FID" = covar_file$eid, "IID" = covar_file$eid, "P1" = 0)
pheno_df_simple <- data.frame("FID" = simple_covar_file$eid, "IID" = simple_covar_file$eid, "P1" = 0)


colnames(covar_file) <- c("FID", "IID", paste0("V", 1:(ncol(covar_file)-2)))
colnames(simple_covar_file) <- c("FID", "IID", paste0("V", 1:(ncol(simple_covar_file)-2)))



total_use_eids <- read.table("use_all_eid", stringsAsFactors=F)
pheno_df_simple$P1[pheno_df_simple$FID %in% total_use_eids[,1]] <- 1

write.table(simple_covar_file, paste0("use_covar_simple"), row.names = F, col.names = T, sep = "\t", quote = F)
write.table(pheno_df_simple, paste0("use_pheno_simple"), row.names = F, col.names = T, sep = "\t", quote = F)

pheno_df_simple$P1 <- 0

for(kk in 1:5){
  if(kk == REAL_age_ind) {
    pheno_df$P1[pheno_df$FID %in% all_use_eids[[kk]][,1]] <- 1
    pheno_df_simple$P1[pheno_df_simple$FID %in% all_use_eids[[kk]][,1]] <- 1
  } else {
    pheno_df <- pheno_df[!(pheno_df[,1] %in% all_use_eids[[kk]][,1]),]
    pheno_df_simple <- pheno_df_simple[!(pheno_df_simple[,1] %in% all_use_eids[[kk]][,1]),]

    covar_file <- covar_file[!(covar_file[,1] %in% all_use_eids[[kk]][,1]),]
    simple_covar_file <- simple_covar_file[!(simple_covar_file[,1] %in% all_use_eids[[kk]][,1]),]
  }
}





##################################################





#write.table(output, "use_pheno", row.names = F, col.names = T, sep = "\t", quote = F)
#write.table(covar_file, "use_covar_full", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(simple_covar_file, paste0("use_covar_simple_", REAL_age_ind), row.names = F, col.names = T, sep = "\t", quote = F)
#write.table(pheno_df, "use_pheno", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(pheno_df_simple, paste0("use_pheno_simple_", REAL_age_ind), row.names = F, col.names = T, sep = "\t", quote = F)
