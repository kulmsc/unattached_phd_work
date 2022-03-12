
modsets <- read.table("temp_files/all_modsets", stringsAsFactors=F)
mergeset <- NULL
for(i in 1:nrow(modsets)){
  readset <- read.table(paste0("../mod_sets/", strsplit(modsets[i,1], ".", fixed=T)[[1]][1], "/", modsets[i,1]), stringsAsFactors = F, header = T)
  readset <- readset[,c(3,4,5,7)]
  colnames(readset)[4] <- paste0("BETA_", modsets[i,1])

  if(is.null(mergeset)){
    mergeset <- readset
  } else {
    mergeset <- merge(mergeset, readset, by = "RSID", all = T)
    mergeset$A1.x[is.na(mergeset$A1.x)] <- mergeset$A1[is.na(mergeset$A1.x)]
    mergeset$A2.x[is.na(mergeset$A2.x)] <- mergeset$A2[is.na(mergeset$A2.x)]
    mergeset <- mergeset[,c(1:3, grep("BETA", colnames(mergeset)))]
  }
}

colnames(mergeset)[2:3] <- c("A1", "A2")
for(i in grep("BETA", colnames(mergeset))){
  mergeset[is.na(mergeset[,i]),i] <- 0
}

write.table(mergeset, "temp_files/comboset.txt", row.names = F, col.names = T, quote = F, sep = "\t")
