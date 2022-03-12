
all_files <- list.files("eur_final_results/")

poss_impact_name <- rep(NA, length(all_files))
poss_impact_val <- rep(NA, length(all_files))

for(f in all_files){

  impact <- read.table(paste0("eur_final_results/", f), stringsAsFactors = F, header = T)
  ci_hi <- (impact[,3] + impact[,4] * 1.96)[nrow(impact)]
  ci_lo <- (impact[,3] - impact[,4] * 1.96)[nrow(impact)]
  if(ci_hi > 0 & ci_lo > 0){
    poss_impact_name[which(all_files == f)] <- f
    poss_impact_val[which(all_files == f)] <- impact[nrow(impact),3]
  }

}

best_impact <- poss_impact_name[which.max(poss_impact_val)]
best_impact <- as.numeric(strsplit(best_impact, ".", fixed = T)[[1]][2]) + 4

all_impact <- list()
for(i in 1:22){
  system(paste0("zcat ~/athena/exome_score/funct_anno/impact/IMPACT707_EUR_chr", i, ".annot.gz | cut -f3,", best_impact, " > temp"))
  all_impact[[i]] <- read.table("temp", stringsAsFactors = F, header = T)
}



all_impact <- do.call("rbind", all_impact)

write.table(all_impact[all_impact[,2] > quantile(all_impact[,2], 0.9),1], "eur.top.10.txt", row.names = F, col.names = F, quote = F)
write.table(all_impact[all_impact[,2] > quantile(all_impact[,2], 0.95),1], "eur.top.5.txt", row.names = F, col.names = F, quote = F)
write.table(all_impact[all_impact[,2] > quantile(all_impact[,2], 0.99),1], "eur.top.1.txt", row.names = F, col.names = F, quote = F)
