#arguments are: dir#, ldpred file name, author, chr #, score #
args <- commandArgs(trailingOnly=TRUE)
low_author <- tolower(args[3])

ldpred_ss <- read.table(args[2], stringsAsFactors=F, header = T)
raw_ss <- read.table(paste0("temp_files/ss.", low_author, ".", args[4]), stringsAsFactors=F, header=T)

raw_ss <- raw_ss[raw_ss$RSID %in% ldpred_ss$sid,]
ldpred_ss <- ldpred_ss[order(ldpred_ss$sid)[rank(raw_ss$RSID)],]

raw_ss$BETA <- ldpred_ss$ldpred_beta
write.table(raw_ss, paste0("~/athena/doc_score/mod_sets/", args[3], "/", low_author, ".", args[4],  ".ldpred.", args[5] ,".ss"), row.names = F, col.names = T, sep = '\t', quote = F)
