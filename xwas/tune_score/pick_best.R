
#Read in the starting information
all_author <- unlist(lapply(strsplit(list.files("tune_results/", "res"), "_"), function(x) x[1]))
for(author in all_author){
  all_res <- readRDS(paste0("tune_results/", author, "_res.RDS"))

  #Average all of the results
  mean_conc <- Reduce("+", all_res[["conc"]])/length(all_res[["conc"]])

  mean_survfit <- matrix(0, nrow = length(all_res[["survfit"]][[1]]), ncol = 6)
  for(i in 1:length(all_res[["survfit"]])){
    for(j in 1:length(all_res[["score_names"]])){
      mean_survfit[j,] <- mean_survfit[j,] + as.numeric(tail(all_res[["survfit"]][[i]][[j]],1)[-1])
    }
  }
  mean_survfit <- mean_survfit/length(all_res[["survfit"]])

  mean_auc <- Reduce("+", all_res[["auc"]])/length(all_res[["auc"]])

  mean_or <- Reduce("+", all_res[["or"]])/length(all_res[["or"]])


  #Get the corresponding best score name
  conc_best_name <- all_res[["score_names"]][which.max(mean_conc[,2])]

  survfit_best_name <- all_res[["score_names"]][which.max(mean_survfit[,3]/mean_survfit[,1])]

  auc_best_name <- all_res[["score_names"]][which.max(mean_auc[,2])]

  or_best_name <- all_res[["score_names"]][which.max(mean_or[,2])]


  #Write the result
  to_write <- rbind(conc_best_name, survfit_best_name, auc_best_name, or_best_name)
  write.table(to_write, paste0("tune_results/", tolower(author), ".best.ss"), row.names = T, col.names = F, quote = F, sep = '\t')

}
