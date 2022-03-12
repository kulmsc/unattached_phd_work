library(vroom)
library(data.table)

#author <- "gormley"
author <- commandArgs(trailingOnly=TRUE)[1]

############################# READ IN AND SET UP #########################################

train_df <- readRDS(paste0("../get_data/df_output/df_train.", author, ".RDS"))
test_df <- readRDS(paste0("../get_data/df_output/df_test.", author, ".RDS"))
surv_train_df <- readRDS(paste0("../get_data/df_output/survdf_train.", author, ".RDS"))
surv_test_df <- readRDS(paste0("../get_data/df_output/survdf_test.", author, ".RDS"))

big_data <- readRDS("../get_data/big_data.RDS")
colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))] <- paste0("X", colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))])

big_data <- big_data[big_data$eid %in% train_df$eid,]
train_df <- train_df[train_df$eid %in% big_data$eid,]
big_data <- big_data[order(big_data$eid)[rank(train_df$eid)],]


if(mean(train_df$score[train_df$pheno == 1]) < mean(train_df$score[train_df$pheno == 0])){
  train_df$score <- train_df$score * -1
}


#need to remove phenotypes that occured before enrollment
#note that the big_data only has ICDs, etc. that occured before enrollment
#can compare phenotypes such that where a big data col matches with the phenotype call we remove the person from the analysis
#then also remove all relevant cols just to make sure



pheno_decoder <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_defs", stringsAsFactors=F, header=T)
useful_decoder <- as.character(pheno_decoder[tolower(pheno_decoder[,1]) == author,4:9][1,])
real_names <- c("CANCER", "NONCANCER", "ICD9", "ICD10", "OPCS", "MEDS")
bad_names <- c()
bad_people <- list()
j <- 1

#remove individuals who were already diagnosed
for(i in 1:length(useful_decoder)){
  if(!is.na(useful_decoder[i])){
    if(grepl("|", useful_decoder[i])){
      for(subcode in strsplit(useful_decoder[i], "|", fixed = T)[[1]]){
        bad_names <- c(bad_names, paste0(real_names[i], "_", subcode))   
        if(real_names[i] %in% c("ICD9", "ICD10", "OPCS")){
          bad_people[[j]] <- big_data$eid[big_data[,which(colnames(big_data) == paste0(real_names[i], "_", subcode))] == 1]
          j <- j + 1
        }
      }

    } else {
      bad_names <- c(bad_names, paste0(real_names[i], "_", useful_decoder[i]))
      bad_people[[j]] <- big_data$eid[big_data[,which(colnames(big_data) == paste0(real_names[i], "_", useful_decoder[i]))] == 1]
      j <- j + 1
    }
    
  }
}

bad_people <- unique(unlist(bad_people))
big_data <- big_data[!(big_data$eid %in% bad_people),]
train_df <- train_df[!(train_df$eid %in% bad_people),]

big_data <- big_data[,-which(colnames(big_data) %in% bad_names)]
minor_counts <- readRDS("minor_counts.RDS")
big_data <- big_data[,!(colnames(big_data) %in% names(minor_counts)[minor_counts < 250])]

################################### COMPARE RESIDUALS AND INTERACTIONS ###########################################

base_mod <- glm(pheno ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, train_df, family = "binomial")
resids <- base_mod$residuals


#exit() #already did this
tryit <- function(x){ 
  if(length(unique(x)) == 2){
    return(min(table(x)))
  } 
  return(10000)
}

na_counts <- apply(big_data, 2, function(x) sum(is.na(x)))
big_data <- big_data[,na_counts < nrow(big_data)*0.6]

rho_val <- rep(NA, ncol(big_data))
p_val <- rep(NA, ncol(big_data))
inter_p <- rep(NA, ncol(big_data))

for(i in 1:ncol(big_data)){
  print(paste(i, " - resid"))
  if(sum(!is.na(big_data[,i])) > 100){

    test_ans <- cor.test(base_mod$residuals[!is.na(big_data[,i])], big_data[!is.na(big_data[,i]), i])
    rho_val[i] <- test_ans$est
    p_val[i] <- test_ans$p.value

    train_df$try_out <- big_data[,i]
    inter_mod <- glm(pheno ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score + score:try_out, train_df, family = "binomial")
    if(rownames(summary(inter_mod)$coef)[nrow(summary(inter_mod)$coef)] == "score:try_out"){
      inter_p[i] <- summary(inter_mod)$coef[nrow(summary(inter_mod)$coef), 4]
    }

  }
}

################################### GET CORR BETWEEN PHENO ###########################################

mod_p <- rep(NA, ncol(big_data))
mod_coef <- rep(NA, ncol(big_data))
test_p <- rep(NA, ncol(big_data))

for(i in 1:ncol(big_data)){
  print(paste(i, " - corr"))
  train_df$try_out <- big_data[,i]

  if(sum(is.na(train_df$try_out)) < 50){
    if(length(unique(train_df$try_out[!is.na(train_df$try_out)])) > 1){
      mod <- glm(pheno ~ age + try_out, family = "binomial", data = train_df)
      mod_coef[i] <- summary(mod)$coef[2,1]
      mod_p[i] <- summary(mod)$coef[length(mod$coef),4]
    }

    if(length(unique(train_df$try_out[!is.na(train_df$try_out)])) > 6){ #Continous
      test_p[i] <- wilcox.test(train_df$try_out[train_df$pheno == 0], train_df$try_out[train_df$pheno == 1])$p.value

    } else if(length(unique(train_df$try_out[!is.na(train_df$try_out)])) == 2){ #Binary
      test_p[i] <- fisher.test(table(train_df$try_out, train_df$pheno))$p.value

    } else if(length(unique(train_df$try_out[!is.na(train_df$try_out)])) == 1){
      test_p[i] <- NA

    } else { #Categorical
      test_p[i] <- chisq.test(table(train_df$try_out, train_df$pheno), simulate.p = T)$p.value

    }
  }
}




####################

save_coef <- data.frame("rho" = rho_val, "cor_p" = p_val, "inter" = inter_p,
                        "mod_p" = mod_p, "mod_coef" = mod_p, "test_p" = test_p)


saveRDS(save_coef, paste0("resid_res/", author, ".RDS"))


