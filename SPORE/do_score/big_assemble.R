library(data.table)
author <- "eastasian"

all_files <- list.files("big_small_scores/", paste0(author, ".sscore"))
read_files <- list()

for(i in 1:length(all_files)){
  read_files[[i]] <- as.data.frame(fread(paste0("big_small_scores/", all_files[i])))
  read_files[[i]] <- read_files[[i]][,grepl("_SUM", colnames(read_files[[i]]))]
}

u_colnames <- unique(unlist(lapply(read_files, colnames)))
u_colnames <- u_colnames[u_colnames != "NAMED_ALLELE_DOSAGE_SUM"]
new_cols <- unique(unlist(lapply(strsplit(u_colnames, ".", fixed = T), function(x) paste(x[3:4], collapse="."))))
score_holder <- data.frame(matrix(0, nrow = nrow(read_files[[1]]), ncol = length(new_cols)))


for(i in 1:length(read_files)){
  c_colnames <- unlist(lapply(strsplit(colnames(read_files[[i]]), ".", fixed = T), function(x) paste(x[3:4], collapse=".")))
  for(j in 1:length(c_colnames)){
    if(c_colnames[j] %in% new_cols){
      score_holder[,which(new_cols == c_colnames[j])] <- score_holder[,which(new_cols == c_colnames[j])] + read_files[[i]][,j]
    }
  }
}

colnames(score_holder) <- unlist(lapply(strsplit(paste(new_cols, author, sep="."), ".", fixed = T), function(x) paste(rev(x), collapse=".")))
saveRDS(score_holder, paste0("final_scores/all_score.", tolower(author), ".RDS"))
