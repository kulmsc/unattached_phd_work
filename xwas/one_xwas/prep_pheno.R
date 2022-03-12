args <- commandArgs(trailingOnly=TRUE)

pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]

pheno <- read.table(paste0("../get_pheno/pheno_defs/diag.", args[1], ".txt.gz"), stringsAsFactors=F, header=F)
output <- data.frame(eid = pheno_eids, eid2 = pheno_eids, pheno = apply(pheno, 1, function(x) 1 %in% x[1:4]) * 1 + 1)

train <- read.table("~/athena/doc_score/qc/cv_files/train_eid.0.5.txt", stringsAsFactors=F)
output <- output[output[,1] %in% train[,1],]

write.table(output, "use_pheno", row.names = F, col.names = F, sep = "\t", quote = F)
