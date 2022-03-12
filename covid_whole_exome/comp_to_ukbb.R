
bim_file <- read.table("plink_files/all_merged.bim", stringsAsFactors=F)

for(chr in as.character(1:22)){
  sub_bim <- bim_file[bim_file[,1] == chr,]
  sub_bim <- sub_bim[order(sub_bim[,3]),]
  sub_bim$rsid <- "none"
  bad_pos <- sub_bim[duplicated(sub_bim[,4]), 4]
  sub_bim <- sub_bim[!(sub_bim[,4] %in% bad_pos),]

  ukbb <- readRDS(paste0("~/athena/ukbiobank/qc/imputed/chr", chr, ".RDS"))
  bad_pos <- ukbb[duplicated(ukbb[,3]), 3]
  ukbb <- ukbb[!(ukbb[,3] %in% bad_pos),]

  sub_bim$rsid[sub_bim[,4] %in% ukbb[,3]] <- ukbb[ukbb[,3] %in% sub_bim[,4], 2]
  write.table(sub_bim, paste0("plink_files/matched_", chr, ".bim"), row.names = F, col.names = F, sep = '\t', quote = F)
}
