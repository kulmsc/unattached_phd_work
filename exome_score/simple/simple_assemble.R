library(data.table)
author <- "christophersen"
#should add author and remove brit_eid

all_files <- list.files("small_scores/", pattern = "profile")
all_files <- grep(tolower(author), all_files, value = T)
split_files <- strsplit(all_files, ".", fixed = T)
all_author <- unlist(lapply(split_files, function(x) x[1]))
all_method <- unlist(lapply(split_files, function(x) x[2]))
all_chr <- unlist(lapply(split_files, function(x) x[3]))
other_method <- unlist(lapply(split_files, function(x) x[4]))
all_method <- paste0(all_method, "_", other_method)



brit_eid <- read.table("~/athena/doc_score/do_score/temp_files/brit_eid", stringsAsFactors=F)
new_name <- unique(paste0(all_author, ".", all_method))
all_score <- matrix(0, nrow = nrow(brit_eid), ncol = length(new_name))
almost_new_name <- paste0(all_author, ".", all_method)
colnames(all_score) <- new_name

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

all_upauthor <- firstup(all_author)
all_upauthor[all_upauthor == "Imsgc"] <- "IMSGC"

###################################

for(i in 1:length(all_files)){
      system(paste0("zstd -d small_scores/", all_files[i]))
      sub_score <- as.data.frame(fread(paste0("small_scores/", substr(all_files[i], 1, nchar(all_files[i])-4))))
      system(paste0("rm small_scores/", substr(all_files[i], 1, nchar(all_files[i])-4)))
 
      rel_eids <- brit_eid[brit_eid[,1] %in% sub_score[,1],1]
      sub_score <- sub_score[sub_score[,1] %in% rel_eids,]
      sub_score <- sub_score[order(sub_score[,1])[rank(rel_eids)],]

      all_score[brit_eid[,1] %in% sub_score[,1], new_name == almost_new_name[i]] <- all_score[brit_eid[,1] %in% sub_score[,1], new_name == almost_new_name[i]] + sub_score$SCORESUM 
}

all_score <- all_score[,colSums(all_score) != 0]

next_file <- tolower(author)


write.table(all_score, paste0("final_scores/all_score.", next_file), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(all_score, paste0("final_scores/all_score.", next_file, ".RDS"))
write.table(colnames(all_score), paste0("final_scores/done_names.", next_file), row.names = F, col.names = F, quote = F)
system(paste0("gzip final_scores/all_score.", next_file))
