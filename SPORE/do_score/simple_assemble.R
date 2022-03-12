library(data.table)
author <- "total"
#should add author and remove brit_eid

all_files <- list.files("small_score_files/", pattern = "profile")
all_files <- grep(tolower(author), all_files, value = T)
split_files <- strsplit(all_files, ".", fixed = T)
all_author <- unlist(lapply(split_files, function(x) x[2]))
all_chr <- unlist(lapply(split_files, function(x) x[3]))
all_ind <- unlist(lapply(split_files, function(x) x[4]))
all_method <- unlist(lapply(split_files, function(x) x[5]))

#brit_eid <- read.table("temp_files/brit_eid", stringsAsFactors=F)
npeople <- 487409
new_name <- unique(paste0(all_author, ".", all_ind, ".", all_method))
all_score <- matrix(0, nrow = npeople, ncol = length(new_name))
almost_new_name <- paste0(all_author, ".", all_ind, ".", all_method)
colnames(all_score) <- new_name



for(i in 1:length(all_files)){
      system(paste0("zstd -d small_score_files/", all_files[i]))
      sub_score <- as.data.frame(fread(paste0("small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4))))
      system(paste0("rm small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4)))
 
      all_score[,new_name == almost_new_name[i]] <- all_score[,new_name == almost_new_name[i]] + sub_score$SCORESUM 
}

all_score <- all_score[,colSums(all_score) != 0]

next_file <- paste0(tolower(author), "_CHECK")

write.table(all_score, paste0("final_scores/all_score.", next_file), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(all_score, paste0("final_scores/all_score.", next_file, ".RDS"))
write.table(colnames(all_score), paste0("final_scores/done_names.", next_file), row.names = F, col.names = F, quote = F)
system(paste0("gzip final_scores/all_score.", next_file))
