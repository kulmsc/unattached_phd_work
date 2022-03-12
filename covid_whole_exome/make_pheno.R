library(openxlsx)

#fam <- read.table("all_merged.fam", stringsAsFactors = F)
fam <- read.table("all_qc.psam", stringsAsFactors=F, comment.char="%", header = T)
manifest <- read.xlsx("../Manifes_covid_exome.xlsx", sheet=2, detectDates=T)

manifest[,2] <- as.Date(manifest[,2])
manifest[,4] <- as.Date(manifest[,4])

manifest <- manifest[manifest[,6] %in% fam[,1],]
manifest <- manifest[order(manifest[,6])[rank(fam[,1])],]

#fam[manifest[,7] == 1, 6] <- 2
#fam[manifest[,7] == 0, 6] <- 1

fam[manifest[,3] == "Male", 2] <- 1
fam[manifest[,3] == "Female", 2] <- 2

write.table(fam, "all_qc.new.psam", row.names = F, col.names = F, sep = " ", quote = F)

pheno <- data.frame(0, fam[,1],  manifest[,7]+1)
write.table(pheno, "pheno_for_gwas", row.names = F, col.names = F, quote = F, sep = '\t')

age <- as.numeric(manifest[,4] - manifest[,2])/365
age[age < 1] <- mean(age[age > 1])
covar <- cbind(fam[,1], fam[,1], age)
write.table(covar, "covar_for_gwas", row.names = F, col.names = F, quote = F, sep = '\t')

tfam <- data.frame(pheno[,1], pheno[,1], 0, 0, fam[,2], pheno[,3])
write.table(tafm, "tfam_for_snpeff", row.names = F, col.names = F, quote = F, sep = '\t')
