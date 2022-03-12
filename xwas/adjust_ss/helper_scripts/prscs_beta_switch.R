#arguments are: d, i, author, chr 
args <- commandArgs(trailingOnly=TRUE)
d <- args[1]
i <- args[2]
author <- args[3]
chr <- args[4]
low_author <- tolower(author)

print("starting r script")

file_name <- list.files(d, paste0("output.", i, "_pst"))
beta <- read.table(paste0(d, "/", file_name), stringsAsFactors = F)
ss <- read.table(paste0("temp_files/ss.", low_author, ".", chr), stringsAsFactors=F, header=T)

ss <- ss[ss$RSID %in% beta[,2],]
beta <- beta[order(beta[,2])[rank(ss$RSID)],]
ss$BETA <- beta[,6]

print("writing r script")
print(head(ss))

write.table(ss, paste0("~/athena/xwas/mod_sets/", author, "/", low_author, ".", chr,  ".prscs.", i ,".ss"), row.names = F, col.names = T, sep = '\t', quote = F)

