library(vroom)
options(warn=2)

codes <- read.csv("features_improved_3.csv", stringsAsFactors=F)


get_pheno <- function(pull_ind, do_job){
 
  if(pull_ind == 1){
    pheno_1 <- as.data.frame(vroom("~/athena/ukbiobank/phenotypes/ukb26867.csv.gz", delim = ","))
  } else if(pull_ind == 2){
    pheno_1 <- as.data.frame(vroom("~/athena/ukbiobank/phenotypes/ukb33822.csv.gz", delim = ","))
  } else if(pull_ind == 3){
    pheno_1 <- as.data.frame(vroom("~/athena/ukbiobank/setup_thirdPhenos/ukb41972.csv.gz", delim = ","))
  } else if(pull_ind == 4){
    pheno_1 <- as.data.frame(vroom("~/athena/ukbiobank/setup_fourthPhenos/ukb42385.csv.gz", delim = ","))
  }

  eid_1 <- pheno_1$eid
  pheno <- pheno_1
  conv_codes <- codes[codes[,1] %in% colnames(pheno),]
  use_pheno_1 <- matrix(0, nrow = nrow(pheno), ncol = nrow(conv_codes))


  for(i in 1:ncol(use_pheno_1)){
    print(i)
    print(conv_codes$UDI[i])

    #Just map the direct column over
    if(is.na(conv_codes$Coding[i])){
      use_pheno_1[,i] <- pheno[,which(colnames(pheno) == conv_codes$UDI[i])]
    } else {

      #if there is a comma there are multiple codes to be thought about, need to serperate them
      if(grepl(",", conv_codes$Coding[i])){
        code_vec <- as.numeric(strsplit(conv_codes$Coding[i], ",")[[1]])
      } else {
        code_vec <- as.numeric(conv_codes$Coding[i])
      }

      #if there is a "set to NA" then we use the code_vec to set to NA, otherwise we set to 1
      if(conv_codes$Coding.Description[i] == "set to NA"){
        use_pheno_1[,i] <- pheno[,which(colnames(pheno) == conv_codes$UDI[i])]
        use_pheno_1[use_pheno_1[,i] == -10,i] <- 0
        use_pheno_1[use_pheno_1[,i] %in% code_vec,i] <- NA 
        
      } else {
    
        print(conv_codes$Ending.UDI[i])
        if(conv_codes$Ending.UDI[i] == ""){
          use_pheno_1[pheno[,which(colnames(pheno) == conv_codes$UDI[i])] %in% code_vec,i] <- 1
        } else {
          start_ind <- which(colnames(pheno) == conv_codes$UDI[i])
          end_ind <- which(colnames(pheno) == conv_codes$Ending.UDI[i])
          for(j in start_ind:end_ind){
            use_pheno_1[pheno[,j] %in% code_vec,i] <- 1
          }
        }
      }
    }

    if(conv_codes$remove.minmax[i] == "both"){
      print("minmax")
      print(sum(is.na(use_pheno_1[,i])))
      max_val <- max(use_pheno_1[!is.na(use_pheno_1[,i]),i])
      min_val <- min(use_pheno_1[!is.na(use_pheno_1[,i]),i])      

      use_pheno_1[use_pheno_1[,i] == max_val & !is.na(use_pheno_1[,i]), i] <- NA
      use_pheno_1[use_pheno_1[,i] == min_val & !is.na(use_pheno_1[,i]), i] <- NA
      print(sum(is.na(use_pheno_1[,i])))
    }

  }


  temp <- data.frame(matrix(0, nrow = ncol(use_pheno_1), ncol = 6))
  for(i in 1:nrow(temp)){
    temp[i,1:6] <- as.numeric(summary(use_pheno_1[,i]))[1:6]
    temp[i,7] <- sum(is.na(use_pheno_1[,i]))
    temp[i,8] <- conv_codes$All.Description[i]
  }
  write.table(temp, "check_it", sep = '\t', row.names =F, col.names=F, quote = F)

  if(do_job){
    #should here add one-hot encoding for jobs 132-0.0
    job_vec <- pheno_1[,which(colnames(pheno_1) == "132-0.0")]
    tab_job_vec <- table(job_vec)
    tab_job_vec <- tab_job_vec[tab_job_vec > 100]
    one_hot_vec <- matrix(0, nrow = nrow(use_pheno_1), ncol = length(tab_job_vec))
    for(kk in 1:length(tab_job_vec)){
      one_hot_vec[job_vec == names(tab_job_vec)[kk], kk] <- 1 
    }
    use_pheno_1 <- cbind(use_pheno_1, one_hot_vec)
    new_conv_codes <- paste0("job_", names(tab_job_vec))
    add <- data.frame(matrix("", ncol = ncol(conv_codes), nrow = length(new_conv_codes)))
    colnames(add) <- colnames(conv_codes)
    add$UDI <- new_conv_codes
    conv_codes <- rbind(conv_codes, add)
  }

  return(list(use_pheno_1, conv_codes, eid_1))

}


use_pheno_1 <- get_pheno(1, TRUE)
use_pheno_2 <- get_pheno(2, FALSE)
use_pheno_3 <- get_pheno(3, FALSE)
use_pheno_4 <- get_pheno(4, FALSE)


all_codes <- rbind(use_pheno_1[[2]], use_pheno_2[[2]], use_pheno_3[[2]], use_pheno_4[[2]])

eid_1 <- use_pheno_1[[3]]
eid_2 <- use_pheno_2[[3]]
eid_3 <- use_pheno_3[[3]]
eid_4 <- use_pheno_4[[3]]

use_pheno_1 <- use_pheno_1[[1]]
use_pheno_2 <- use_pheno_2[[1]]
use_pheno_3 <- use_pheno_3[[1]]
use_pheno_4 <- use_pheno_4[[1]]

use_pheno_1 <- use_pheno_1[eid_1 %in% eid_2 & eid_1 %in% eid_3 & eid_1 %in% eid_4,]
use_pheno_2 <- use_pheno_2[eid_2 %in% eid_1 & eid_2 %in% eid_3 & eid_2 %in% eid_4,]
use_pheno_3 <- use_pheno_3[eid_3 %in% eid_2 & eid_3 %in% eid_1 & eid_3 %in% eid_4,]
use_pheno_4 <- use_pheno_4[eid_4 %in% eid_1 & eid_4 %in% eid_2 & eid_4 %in% eid_3,]

eid_1 <- eid_1[eid_1 %in% eid_2 & eid_1 %in% eid_3 & eid_1 %in% eid_4]
eid_2 <- eid_2[eid_2 %in% eid_1 & eid_2 %in% eid_3 & eid_2 %in% eid_4]
eid_3 <- eid_3[eid_3 %in% eid_2 & eid_3 %in% eid_1 & eid_3 %in% eid_4]
eid_4 <- eid_4[eid_4 %in% eid_1 & eid_4 %in% eid_2 & eid_4 %in% eid_3]

use_pheno_2 <- use_pheno_2[order(eid_2)[rank(eid_1)],]
use_pheno_3 <- use_pheno_3[order(eid_3)[rank(eid_1)],]
use_pheno_4 <- use_pheno_4[order(eid_4)[rank(eid_1)],]

final <- cbind(use_pheno_1, use_pheno_2, use_pheno_3, use_pheno_4)


codes$Coding.Description[is.na(codes$Coding.Description)] <- "NA"
all_codes$Coding.Description[is.na(all_codes$Coding.Description)] <- "NA"
send_out <- data.frame(matrix(0, nrow = length(eid_1), ncol = nrow(codes)))
colnames(send_out) <- codes[,1]
for(i in 1:nrow(codes)){
  print(i)
  send_out[,i] <- final[,which(codes$UDI[i] == all_codes$UDI & codes$Coding.Description[i] == all_codes$Coding.Description)[1]]
}
send_out$eid <- eid_1

write.table(send_out, "imptoved_features.txt", row.names = F, col.names = T, sep = '\t', quote = F)
