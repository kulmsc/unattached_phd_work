library(data.table)
options(warn=2)
author <- "mdd"
#should add author and remove brit_eid

all_files <- list.files(paste0(author, "_small_score_files/"), pattern = "profile")
all_files <- grep(tolower(author), all_files, value = T)
split_files <- strsplit(all_files, ".", fixed = T)
all_author <- unlist(lapply(split_files, function(x) x[2]))
all_chr <- unlist(lapply(split_files, function(x) x[3]))
all_ind <- unlist(lapply(split_files, function(x) x[4]))
all_method <- unlist(lapply(split_files, function(x) x[5]))

brit_eid <- read.table("temp_files/fam23", stringsAsFactors=F, header=T)
new_name <- unique(paste0(all_author, ".", all_ind, ".", all_method))
all_score_22 <- matrix(0, nrow = nrow(brit_eid), ncol = length(new_name))
all_score_23_male <- matrix(0, nrow = nrow(brit_eid), ncol = length(new_name))
all_score_23_female <- matrix(0, nrow = nrow(brit_eid), ncol = length(new_name))
all_score_23_log1 <- matrix(0, nrow = nrow(brit_eid), ncol = length(new_name))
all_score_23_log2 <- matrix(0, nrow = nrow(brit_eid), ncol = length(new_name))
all_score_23_fish <- matrix(0, nrow = nrow(brit_eid), ncol = length(new_name))

almost_new_name <- paste0(all_author, ".", all_ind, ".", all_method)
colnames(all_score_22) <- new_name
colnames(all_score_23_male) <- new_name
colnames(all_score_23_female) <- new_name
colnames(all_score_23_log1) <- new_name
colnames(all_score_23_log2) <- new_name
colnames(all_score_23_fish) <- new_name

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

#################

for(i in 1:length(all_files)){
      system(paste0("zstd -d ", author, "_small_score_files/", all_files[i]))
      sub_score <- as.data.frame(fread(paste0(author, "_small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4))))
      system(paste0("rm ", author, "_small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4)))
 
      sub_score <- sub_score[sub_score[,1] %in% brit_eid[,1],]
      sub_score <- sub_score[order(sub_score[,1])[rank(brit_eid[,1])],]

      if(!grepl("23", all_files[i])){
        all_score_22[,new_name == almost_new_name[i]] <- all_score_22[,new_name == almost_new_name[i]] + sub_score$SCORESUM 
        all_score_23_male[,new_name == almost_new_name[i]] <- all_score_23_male[,new_name == almost_new_name[i]] + sub_score$SCORESUM
        all_score_23_female[,new_name == almost_new_name[i]] <- all_score_23_female[,new_name == almost_new_name[i]] + sub_score$SCORESUM
        all_score_23_log1[,new_name == almost_new_name[i]] <- all_score_23_log1[,new_name == almost_new_name[i]] + sub_score$SCORESUM
        all_score_23_log2[,new_name == almost_new_name[i]] <- all_score_23_log2[,new_name == almost_new_name[i]] + sub_score$SCORESUM
        all_score_23_fish[,new_name == almost_new_name[i]] <- all_score_23_fish[,new_name == almost_new_name[i]] + sub_score$SCORESUM
      }
      if(grepl("23", all_files[i])){
        if(grepl("23_male", all_files[i])){
          all_score_23_male[,new_name == almost_new_name[i]] <- all_score_23_male[,new_name == almost_new_name[i]] + sub_score$SCORESUM
        } else if(grepl("23_female", all_files[i])){
          all_score_23_female[,new_name == almost_new_name[i]] <- all_score_23_female[,new_name == almost_new_name[i]] + sub_score$SCORESUM
        } else if(grepl("23_log1", all_files[i])){
          all_score_23_log1[,new_name == almost_new_name[i]] <- all_score_23_log1[,new_name == almost_new_name[i]] + sub_score$SCORESUM
        } else if(grepl("23_log2", all_files[i])){
          all_score_23_log2[,new_name == almost_new_name[i]] <- all_score_23_log2[,new_name == almost_new_name[i]] + sub_score$SCORESUM
        } else if(grepl("23_fish", all_files[i])){
          all_score_23_fish[,new_name == almost_new_name[i]] <- all_score_23_fish[,new_name == almost_new_name[i]] + sub_score$SCORESUM
        }
      }
}


next_file <- tolower(author)

write.table(all_score_22, paste0("final_scores/all_score_22.", next_file), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(all_score_22, paste0("final_scores/all_score_22.", next_file, ".RDS"))
write.table(colnames(all_score_22), paste0("final_scores/done_names_22.", next_file), row.names = F, col.names = F, quote = F)
system(paste0("gzip final_scores/all_score_22.", next_file))

write.table(all_score_23_male, paste0("final_scores/all_score_23_male.", next_file), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(all_score_23_male, paste0("final_scores/all_score_23_male.", next_file, ".RDS"))
write.table(colnames(all_score_23_male), paste0("final_scores/done_names_23_male.", next_file), row.names = F, col.names = F, quote = F)
system(paste0("gzip final_scores/all_score_23_male.", next_file))

write.table(all_score_23_female, paste0("final_scores/all_score_23_female.", next_file), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(all_score_23_female, paste0("final_scores/all_score_23_female.", next_file, ".RDS"))
write.table(colnames(all_score_23_female), paste0("final_scores/done_names_23_female.", next_file), row.names = F, col.names = F, quote = F)
system(paste0("gzip final_scores/all_score_23_female.", next_file))

write.table(all_score_23_log1, paste0("final_scores/all_score_23_log1.", next_file), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(all_score_23_log1, paste0("final_scores/all_score_23_log1.", next_file, ".RDS"))
write.table(colnames(all_score_23_log1), paste0("final_scores/done_names_23_log1.", next_file), row.names = F, col.names = F, quote = F)
system(paste0("gzip final_scores/all_score_23_log1.", next_file))

write.table(all_score_23_log2, paste0("final_scores/all_score_23_log2.", next_file), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(all_score_23_log2, paste0("final_scores/all_score_23_log2.", next_file, ".RDS"))
write.table(colnames(all_score_23_log2), paste0("final_scores/done_names_23_log2.", next_file), row.names = F, col.names = F, quote = F)
system(paste0("gzip final_scores/all_score_23_log2.", next_file))

write.table(all_score_23_fish, paste0("final_scores/all_score_23_fish.", next_file), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(all_score_23_fish, paste0("final_scores/all_score_23_fish.", next_file, ".RDS"))
write.table(colnames(all_score_23_fish), paste0("final_scores/done_names_23_fish.", next_file), row.names = F, col.names = F, quote = F)
system(paste0("gzip final_scores/all_score_23_fish.", next_file))
