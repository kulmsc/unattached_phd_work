library(data.table)

args <- commandArgs(trailingOnly=TRUE) #$chrom $y_index $dir_name

res <- as.data.frame(fread(paste0("../", args[3], "/full_res_Y", args[2], ".txt.gz"), stringsAsFactors=F, header=T))

res$P <- 10^-res$LOG10P

res <- res[res$P < 5e-8 & !is.na(res$P),]

keep_snps <- c()

while(nrow(res) > 0){

  keep_snps <- c(keep_snps, res$ID[which.min(res$P)])

  bad_inds <- which(res$GENPOS < res$GENPOS[which.min(res$P)]+250000 & res$GENPOS > res$GENPOS[which.min(res$P)]-250000 & res$CHROM == res$CHROM[which.min(res$P)])

  res <- res[-bad_inds,]
}

if(length(keep_snps) > 0){
  write.table(keep_snps, paste0("peak_snps.", args[3], ".", args[2]), row.names = F, col.names = F, quote = F)
}
