
output <- read.table("use_pheno")

fam <- read.table("toxwas.fam", stringsAsFactors=F, header=F)
output <- output[output[,1] %in% fam[,1],]
output <- output[order(output[,1])[rank(fam[,1])],]
fam[,6] <- output[,3]

write.table(fam, "toxwas.fam", row.names = F, col.names = F, sep = "\t", quote = F)

bim <- read.table("toxwas.bim", stringsAsFactors=F)
bim[,1] <- 23
write.table(bim, "toxwas.bim", row.names = F, col.names = F, quote = F, sep = "\t")


