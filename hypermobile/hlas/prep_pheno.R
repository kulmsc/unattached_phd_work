args <- commandArgs(trailingOnly=TRUE)

pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]

if(args[1] > 2){
  use_arg <- as.numeric(args[1]) - 2
} else {
  use_arg <- as.numeric(args[1])
}

use_eid <- read.table(paste0("../eid_gwas_", use_arg), stringsAsFactors=F)
output <- data.frame(pheno_eids, pheno_eids, 1, stringsAsFactors = F)
output[output[,1] %in% use_eid[,1], 3] <- 2

if(args[1] == 3 | args[1] == 4){
  other_eid <- read.table("../eid_ankspond_3", stringsAsFactors=F)
  output <- output[output[,1] %in% other_eid[,1] | output[,2] %in% use_eid[,1],]
  write.table(output[,1], "use_eid", row.names = F, col.names = F, quote = F)
}

#train <- read.table("~/athena/doc_score/qc/cv_files/train_eid.0.5.txt", stringsAsFactors=F)
#output <- output[output[,1] %in% train[,1],]

write.table(output, "use_pheno", row.names = F, col.names = F, sep = "\t", quote = F)
