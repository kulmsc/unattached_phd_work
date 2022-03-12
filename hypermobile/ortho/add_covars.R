library(vroom)

################################################################
#                    GET PHENO
###############################################################


pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]


the_files <- list.files(".", "eid_gwas")
output <- data.frame(pheno_eids, pheno_eids, matrix(0, nrow = length(pheno_eids), ncol = length(the_files)), stringsAsFactors = F)
for(i in 1:length(the_files)){
  use_eid <- read.table(paste0("./", the_files[i]), stringsAsFactors=F)
  output[output[,1] %in% use_eid[,1], i+2] <- 1
}


#######################################################
#                          GET COVAR
######################################################

csv_file <- as.data.frame(vroom("~/athena/ukbiobank/phenotypes/ukb26867.csv.gz", delim = ","))
brit_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/brit_fam", stringsAsFactors=F)
alive_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/qc_alive_fam", stringsAsFactors=F)

pcs_ind <- 1310:1319
age_ind <- 1300
sex_ind <- 5
bmi_ind <- 1291
smoke_ind <- 1265 #NA to 0

covar_file <- csv_file[,c(1, 1, pcs_ind, age_ind, sex_ind, bmi_ind, smoke_ind)]
covar_file[is.na(covar_file[,ncol(covar_file)]), ncol(covar_file)] <- 0
rm(csv_file)

covar_file <- covar_file[covar_file[,1] %in% brit_fam[,1] & covar_file[,1] %in% alive_fam[,1],]
covar_file <- covar_file[complete.cases(covar_file),]

for(i in 3:ncol(covar_file)){
  if(max(covar_file[,i]) == 1 & min(covar_file[,i]) == 0){
    print("skip")
  } else {
    covar_file[,i] <- (covar_file[,i] - min(covar_file[,i]))/(max(covar_file[,i]) - min(covar_file[,i]))
  }
}

############################################################
#                           FINISH
###########################################################

output <- output[output[,1] %in% covar_file[,1],]
covar_file <- covar_file[covar_file[,1] %in% output[,1],]

output <- output[order(output[,1])[rank(covar_file[,1])],]

colnames(output) <- c("FID", "IID", paste0("Y", 1:(ncol(output)-2)))
colnames(covar_file) <- c("FID", "IID", paste0("V", 1:(ncol(covar_file)-2)))

diabetes_eids <- read.table("diabetes_eids", stringsAsFactors=F)
hypo_eids <- read.table("hypo_eids", stringsAsFactors=F)
covar_file$diabetes <- 0
covar_file$hypo <- 0
covar_file$diabetes[covar_file$FID %in% diabetes_eids[,1]] <- 1
covar_file$hypo[covar_file$FID %in% hypo_eids[,1]] <- 1

#write.table(output, "use_pheno", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(covar_file, "use_diabetes_hypo_covar", row.names = F, col.names = T, sep = "\t", quote = F)

