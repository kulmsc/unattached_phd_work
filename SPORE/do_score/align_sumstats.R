library(data.table)

args = commandArgs(trailingOnly=TRUE)

chr <- as.character(args[1])
author <- as.character(args[2])

mod_sets <- list()
subset_rsids <- list()

ss_names <- list.files(paste0("../mod_sets/", author, "/"), paste0(tolower(author), ".", chr))
ss_names <- grep(paste0(tolower(author), ".", chr, "."), ss_names, fixed = T, value = T)

for(i in 1:length(ss_names)){
  mod_sets[[i]] <- as.data.frame(fread(paste0("../mod_sets/", author, "/", ss_names[i])))
  if(!grepl("rs", mod_sets[[i]][1,3])){
    mod_sets[[i]] <- mod_sets[[i]][,c(1,2,5,3,4,6:9)]
  }
  subset_rsids[[i]] <- mod_sets[[i]][,3]
}



master_ms <- data.frame("rsid" = unlist(lapply(mod_sets, function(x) x[,3])), "A1" = unlist(lapply(mod_sets, function(x) x[,4])), stringsAsFactors=F)
master_ms <- master_ms[!duplicated(master_ms[,1]),]



betas <- list()
name_beta <- rep("", length(mod_sets) * (1+length(subset_rsids)))
for(i in 1:length(mod_sets)){
    print(paste(i, "of", length(mod_sets)))

    betas[[i]] <- rep(0, nrow(master_ms))
    rele_rsids <- master_ms[master_ms[,1] %in% mod_sets[[i]][,3],1]
    curr_mod_set <- mod_sets[[i]][order(mod_sets[[i]][,3])[rank(rele_rsids)],]
    betas[[i]][master_ms[,1] %in% curr_mod_set[,3]] <- curr_mod_set[,7]    
    betas[[i]][is.na(betas[[i]])] <- 0
    betas[[i]][is.infinite(betas[[i]])] <- 0
}


betas <- data.frame(do.call("cbind", betas))
colnames(betas) <- ss_names

master_ms <- data.frame(master_ms, betas)
write.table(master_ms, "big_mod_set", row.names = F, col.names = T, quote = F, sep = "\t")
