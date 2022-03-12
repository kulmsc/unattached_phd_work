
#meta_stats <- read.csv("../raw_ss/meta_stats", stringsAsFactors=F)

#for(author in meta_stats[,1]){

#  final_scores <- readRDS(paste0("final_scores/all_score.", tolower(author), ".RDS"))
#  final_scores <- colnames(final_scores)

#  poss_files <- list.files(paste0("../mod_sets/", author, "/"))
#  poss_files <- unique(unlist(lapply(strsplit(poss_files, ".", fixed = T), function(x) paste0(x[1], ".", x[4], ".", x[3]))))

#  not_done <- poss_files[!poss_files %in% final_scores]
#  write.table(not_done, paste0("check_it/not_done.", tolower(author), ".txt"), row.names = F, col.names = F, quote = F) 

#}


authors <- c("IMSGC", "Mahajan", "Malik", "Nikpay", "Okada")

for(author in authors){

done_files <- grep(tolower(author), list.files("small_score_files/", "zst"), value = T)
to_do_files <- list.files(paste0("../mod_sets/", author))

total_to_do <- 0
total_done <- 0

check_it <- read.table(paste0("check_it/not_done.", tolower(author), ".txt"), stringsAsFactors = F)
for(i in 1:nrow(check_it)){
  to_do_comp <- paste0(".", strsplit(check_it[i,], ".", fixed=T)[[1]][3], ".", strsplit(check_it[i,], ".", fixed=T)[[1]][2], ".ss")
  done_comp <- paste0(".", strsplit(check_it[i,], ".", fixed=T)[[1]][2], ".", strsplit(check_it[i,], ".", fixed=T)[[1]][3], ".profile")

  total_to_do <- total_to_do + length(grep(to_do_comp, to_do_files, fixed = T))
  total_done <-  total_done + length(grep(done_comp, done_files, fixed = T))
}

print(author)
print(total_done/total_to_do)
}
