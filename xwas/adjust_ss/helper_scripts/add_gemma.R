args <- commandArgs(trailingOnly=TRUE)
lowauthor <- tolower(args[1])

bim <- read.table(paste0("geno_files/", lowauthor, ".", args[2], ".bim"), stringsAsFactors=F)
ss <- read.table(paste0("temp_files/ss.", lowauthor, ".", args[2]), stringsAsFactors=F, header=T)
ss <- ss[ss$RSID %in% bim[,2],]
frq <- read.table(paste0("temp_files/", lowauthor, ".", args[2], ".frq"), stringsAsFactors=F, header=T)

ss$frq <- rep(0, nrow(frq))
ss$frq[ss$A1 == frq$A1] <- frq$MAF[ss$A1 == frq$A1]
ss$frq[ss$A2 == frq$A1] <- frq$MAF[ss$A2 == frq$A1]

ss$nmis <- 0

ss <- ss[,c(1,3,2,11,9,4,5,10,7,6,8)]
colnames(ss) <- c("chr", "rs", "ps", "n_mis", "n_obs", "allele1", "allele0", "af", "beta", "se", "p_wald")
ss <- ss[ss$beta != 0 & ss$se != 0,]
write.table(ss, paste0("temp_files/ss.", lowauthor, ".", args[2],".gemma"), row.names = F, col.names = T, sep = '\t', quote = F)
