
basic_data <- read.table("../get_pheno/basic_pheno_data", stringsAsFactors=F, header=T, sep=",")
pheno <- read.table("use_pheno", stringsAsFactors=F)

basic_data[is.na(basic_data[,4]),4] <- mean(basic_data[,4], na.rm=T)
basic_data[is.na(basic_data[,5]),5] <- mean(basic_data[,5], na.rm=T)
basic_data[,4] <- (basic_data[,4] - mean(basic_data[,4]))/sd((basic_data[,4] - mean(basic_data[,4])))
basic_data[,5] <- (basic_data[,5] - mean(basic_data[,5]))/sd((basic_data[,5] - mean(basic_data[,5])))

case_data <- basic_data[basic_data[,1] %in% pheno[pheno[,3] == 2,1],]
control_data <- basic_data[basic_data[,1] %in% pheno[pheno[,3] == 1,1],]

done_eid <- rep(1, nrow(control_data)*3)
start_ind <- 1
for(i in 1:nrow(case_data)){
  sub_cc <- control_data[control_data[,2] == case_data[i,2] & control_data[,3] == case_data[i,3] & !(control_data[,1] %in% done_eid[1:start_ind]),]
  if(nrow(sub_cc) > 3){
    diff_val <- abs(sub_cc[,4] - case_data[i,4]) + abs(sub_cc[,5] - case_data[i,5])
    done_eid[start_ind:(start_ind+2)] <- sub_cc[order(diff_val)[1:3],1]
    start_ind <- start_ind + 3
  } else if(nrow(sub_cc) > 0){
    done_eid[start_ind:(start_ind+nrow(sub_cc)-1)] <- sub_cc[,1]
    start_ind <- start_ind + nrow(sub_cc)
  }
}

done_eid <- c(done_eid[done_eid != 1], case_data[,1])
write.table(done_eid, "use_eid", row.names = F, col.names = F, quote = F)
write.table(pheno[pheno[,1] %in% done_eid,], "use_pheno", row.names = F, col.names = F, quote = F)
