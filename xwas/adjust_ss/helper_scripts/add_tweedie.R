args <- commandArgs(trailingOnly=TRUE)
lowauthor <- tolower(args[1])

bim <- read.table(paste0("geno_files/", lowauthor, ".", args[2], ".bim"), stringsAsFactors=F)
ss <- read.table(paste0("temp_files/ss.", lowauthor, ".", args[2]), stringsAsFactors=F, header=T)
ss <- ss[ss$RSID %in% bim[,2],]
frq <- read.table(paste0(args[3], "/freq.frq"), stringsAsFactors=F, header=T)

ss$frq <- rep(0, nrow(frq))
ss$frq[ss$A1 == frq$A1] <- frq$MAF[ss$A1 == frq$A1]
ss$frq[ss$A2 == frq$A1] <- frq$MAF[ss$A2 == frq$A1]

clump <- read.table(paste0(args[3], "/out.clumped"), stringsAsFactors=F, header = T)
ss <- ss[ss$RSID %in% clump$SNP,]

write.table(ss, paste0(args[3], "/ss"), row.names = F, col.names = T, sep = '\t', quote = F)
