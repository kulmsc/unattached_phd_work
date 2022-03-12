chr <- commandArgs(trailingOnly=TRUE)

into <- read.table(paste0("into.", chr, ".vcf"), stringsAsFactors=F)
bed <- read.table(paste0("out.", chr, ".bed"), stringsAsFactors=F)

bed <- bed[bed[,1] == paste0("chr", chr),]

into <- into[into[,3] %in% bed[,4],]
into <- into[order(into[,3])[rank(bed[,4])],]
into[,2] <- bed[,2]

write.table(into, "temp", row.names = F, col.names = F, quote = F, sep = '\t')
