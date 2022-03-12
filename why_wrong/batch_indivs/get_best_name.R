all_author <- unlist(lapply(strsplit(list.files("res", "base"), ".", fixed = T), function(x) x[2]))

a_type <- c("all_around", "mean", "min", "count")
b_type <- c("all_remove", "check_pheno", "explicit")

good_files <- rep("", length(all_author) * length(a_type))
i <- 1

for(author in all_author){
  for(a in a_type){

    all_auc <- rep(0, length(b_type))
    for(b in b_type){
      if(file.exists(paste0("res/togo.", a, ".", author, ".", b, ".RDS"))){
        print(paste0("res/togo.", a, ".", author, ".", b, ".RDS"))
        print(paste(x[[4]], x[[5]]))
       
      x <- readRDS(paste0("res/togo.", a, ".", author, ".", b, ".RDS")) 
      all_auc[which(b_type == b)] <- max(x[[3]])
      }
    }

    print(all_auc)
    good_files[i] <- paste0("res/togo.", a, ".", author, ".", b_type[which.max(all_auc)], ".RDS")
    i <- i + 1
  }

}
