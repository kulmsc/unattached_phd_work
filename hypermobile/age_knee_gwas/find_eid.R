pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]

file_names <- list.files("../comboortho/pheno_files/", "time")
code_names <- unlist(lapply(strsplit(list.files("../comboortho/pheno_files/", "time"), ".", fixed = T), function(x) x[2]))
list_data <- list()
for(i in 1:length(file_names)){
  list_data[[i]] <- read.table(paste0("../comboortho/pheno_files/", file_names[i]), stringsAsFactors = F)

  for(j in 1:5){
    if(any(grepl("/", unique(list_data[[i]][list_data[[i]][,j] != "__________", j])))){
      list_data[[i]][,j] <- as.Date(list_data[[i]][,j], "%d/%m/%Y")
    } else {

      list_data[[i]][,j] <- as.Date(list_data[[i]][,j], "%Y-%m-%d")
    }
  }

  if(grepl("icd", code_names[i])){
    list_data[[i]] <- list_data[[i]][,3]
  } else {
    list_data[[i]] <- list_data[[i]][,4]
  }
}
names(list_data) <- code_names

system("zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -f1,6,40 -d',' > age_data")
age_data <- read.table("age_data", stringsAsFactors=F, header=T, sep=",")
age_data$dob <- as.Date(paste0("15-", age_data[,3], "-", age_data[,2]), "%d-%m-%Y")
age_data <- age_data[age_data$eid %in% pheno_eids,]
age_data <- age_data[order(age_data$eid)[rank(pheno_eids)],]


#######################################################################################################################
exit()
system("zcat ~/athena/ukbiobank/phenotypes/ukb26867.csv.gz | cut -f1,1291 -d',' > bmi_data")
bmi_data <- read.table("bmi_data", stringsAsFactors=F, header=T, sep=",")
bmi_data <- bmi_data[bmi_data$eid %in% pheno_eids,]
bmi_data <- bmi_data[order(bmi_data$eid)[rank(pheno_eids)],]
w401_1 <- which(!is.na(list_data[["opcs_w401"]]) & bmi_data[,2] <= 25)
w401_2 <- which(!is.na(list_data[["opcs_w401"]]) & bmi_data[,2] > 25 & bmi_data[,2] < 30)
w401_3 <- which(!is.na(list_data[["opcs_w401"]]) & bmi_data[,2] >= 30)

w401_1 <- pheno_eids[w401_1]
w401_2 <- pheno_eids[w401_2]
w401_3 <- pheno_eids[w401_3]

write.table(w401_1, "eid_gwas_6", row.names = F, col.names = F)
write.table(w401_2, "eid_gwas_7", row.names = F, col.names = F)
write.table(w401_3, "eid_gwas_8", row.names = F, col.names = F)



#401, 411, 422
#w401_1 <- which((list_data[["opcs_w401"]] - age_data$dob)/365 < 55)
#w401_2 <- which((list_data[["opcs_w401"]] - age_data$dob)/365 >= 55 & (list_data[["opcs_w401"]] - age_data$dob)/365 < 60)
#w401_3 <- which((list_data[["opcs_w401"]] - age_data$dob)/365 >= 60 & (list_data[["opcs_w401"]] - age_data$dob)/365 < 65)
#w401_4 <- which((list_data[["opcs_w401"]] - age_data$dob)/365 >= 65 & (list_data[["opcs_w401"]] - age_data$dob)/365 < 70)
#w401_5 <- which((list_data[["opcs_w401"]] - age_data$dob)/365 > 70)

#w401_1 <- pheno_eids[w401_1]
#w401_2 <- pheno_eids[w401_2]
#w401_3 <- pheno_eids[w401_3]
#w401_4 <- pheno_eids[w401_4]
#w401_5 <- pheno_eids[w401_5]
#w401_all <- pheno_eids[which(!is.na(list_data[["opcs_w401"]]))]

#write.table(w401_1, "eid_gwas_1", row.names = F, col.names = F)
#write.table(w401_2, "eid_gwas_2", row.names = F, col.names = F)
#write.table(w401_3, "eid_gwas_3", row.names = F, col.names = F)
#write.table(w401_4, "eid_gwas_4", row.names = F, col.names = F)
#write.table(w401_5, "eid_gwas_5", row.names = F, col.names = F)

#write.table(w401_all, "use_all_eid", row.names = F, col.names = F)

