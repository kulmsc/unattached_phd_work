library(vroom)
library(data.table)

author <- "gormley"
#author <- commandArgs(trailingOnly=TRUE)[1]

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


hi_prs_no_disease_eid <- train_df$eid[train_df$pheno == 0][train_df$score[train_df$pheno == 0] >= sort(train_df$score[train_df$pheno == 0], decreasing=T)[100]]
hi_prs_yes_disease_eid <- train_df$eid[train_df$pheno == 1][train_df$score[train_df$pheno == 1] >= sort(train_df$score[train_df$pheno == 1], decreasing=T)[100]]
lo_prs_no_disease_eid <- train_df$eid[train_df$pheno == 0][train_df$score[train_df$pheno == 0] <= sort(train_df$score[train_df$pheno == 0], decreasing=F)[100]]
lo_prs_yes_disease_eid <- train_df$eid[train_df$pheno == 1][train_df$score[train_df$pheno == 1] <= sort(train_df$score[train_df$pheno == 1], decreasing=F)[100]]



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


################################### RUN SIMPLE TESTS ###########################################

############ Univariate

#run_regressions <- function(ppl_1, ppl_2, comp){
ppl_1 <- hi_prs_no_disease_eid
ppl_2 <- hi_prs_yes_disease_eid
comp <- "hi_lo_all"

run_df <- train_df[train_df$eid %in% c(ppl_1, ppl_2),]
run_df$y <- 0
run_df$y[run_df$eid %in% ppl_1] <- 1 ################################# NO DISEASE = 1
small_data <- big_data[big_data$eid %in% c(ppl_1, ppl_2),]

mod_p <- rep(NA, ncol(small_data))
mod_coef <- rep(NA, ncol(small_data))
test_p <- rep(NA, ncol(small_data))

for(i in 1:ncol(small_data)){
  run_df$x <- small_data[,i]

  if(sum(is.na(run_df$x)) < 50){
    if(length(unique(run_df$x[!is.na(run_df$x)])) > 1){
      mod <- glm(y ~ age + x, family = "binomial", data = run_df)
      mod_coef[i] <- summary(mod)$coef[2,1]
      mod_p[i] <- summary(mod)$coef[length(mod$coef),4]
    }

    if(length(unique(run_df$x[!is.na(run_df$x)])) > 6){ #Continous
      test_p[i] <- wilcox.test(run_df$x[run_df$y == 0], run_df$x[run_df$y == 1])$p.value

    } else if(length(unique(run_df$x[!is.na(run_df$x)])) == 2){ #Binary
      test_p[i] <- fisher.test(table(run_df$x, run_df$y))$p.value

    } else if(length(unique(run_df$x[!is.na(run_df$x)])) == 1){
      test_p[i] <- NA

    } else { #Categorical
      test_p[i] <- chisq.test(table(run_df$x, run_df$y), simulate.p = T)$p.value

    }
  }
}

############### Multivariate

multi_by_p <- TRUE
if(multi_by_p){
  test_p[is.na(test_p)] <- 1
  mod_p[is.na(mod_p)] <- 1

  use_p <- as.numeric(mod_p)
  use_p[test_p > mod_p] <- mod_p[test_p > mod_p]
  #use_p <- test_p * mod_p
  
  multi_data <- small_data[,use_p < sort(use_p)[11]]
  #test_p <- test_p[use_p < sort(use_p)[11]]
  #mod_p <- mod_p[use_p < sort(use_p)[11]]
  use_p <- use_p[use_p < sort(use_p)[11]]
  print("DIM-1")
  print(length(use_p))


  for(i in 1:ncol(multi_data)){
    multi_data[,i][is.na(multi_data[,i])] <- median(multi_data[,i], na.rm=T)
  }
  multi_df <- cbind(run_df, multi_data)
  multi_df <- multi_df[,colnames(multi_df) %in% c("y", "age", colnames(multi_data))]

} else {

  multi_data <- small_data[,apply(small_data, 2, function(x) length(unique(x))) > 1]
  multi_data <- multi_data[,apply(multi_data, 2, function(x) sum(is.na(x)) < 10)]
  for(i in 1:ncol(multi_data)){
    multi_data[,i][is.na(multi_data[,i])] <- median(multi_data[,i], na.rm=T)
  }

  mod <- glmnet(y = run_df$y, x = as.matrix(multi_data), lambda = 0.05)
  multi_df <- cbind(run_df, multi_data[,as.numeric(mod$beta) != 0])
  multi_df <- multi_df[,colnames(multi_df) %in% c("y", "age", colnames(multi_data[,as.numeric(mod$beta) != 0]))]
}

#remove high corr vals (within the multi_df)
corr_mat <- cor(multi_df[,-2]) #remove y
diag(corr_mat) <- 0
cols_to_remove <- c()
for(i in 1:nrow(corr_mat)){
  if(!colnames(corr_mat)[i] %in% cols_to_remove){
    cols_to_remove <- c(cols_to_remove, colnames(corr_mat)[abs(corr_mat[i,]) > 0.9])
  }
}

#normalize
if(length(multi_df[,-which(colnames(multi_df) %in% cols_to_remove)]) > 0){
  multi_df <- multi_df[,-which(colnames(multi_df) %in% cols_to_remove)]
}
for(i in 1:ncol(multi_df)){
  if(!all(unique(multi_df[,i]) %in% c(0,1))){
    multi_df[,i] <- (multi_df[,i] - min(multi_df[,i]))/(max(multi_df[,i]) - min(multi_df[,i]))
  }
}

exit()
mod <- glm(y ~ ., family = "binomial", data = multi_df)
save_coef <- as.data.frame(summary(mod)$coef)
save_coef$uni_mod_p <- NA
save_coef$uni_mod_coef <- NA
save_coef$uni_test_p <- NA

save_coef$uni_mod_p[rownames(save_coef) %in% colnames(small_data)] <- mod_p[colnames(small_data) %in% rownames(save_coef)]
save_coef$uni_mod_coef[rownames(save_coef) %in% colnames(small_data)] <- mod_coef[colnames(small_data) %in% rownames(save_coef)]
save_coef$uni_test_p[rownames(save_coef) %in% colnames(small_data)] <- test_p[colnames(small_data) %in% rownames(save_coef)]

############### Need to change the names of the rownames #################
data_type <- c("ICD10", "ICD9", "NONCANCER", "CANCER", "MEDS", "OPCS", "JOB")
coding_file <- c("coding19.tsv", "coding87.tsv", "coding6.tsv", "coding3.tsv", "coding4.tsv", "coding240.tsv", "coding2.tsv")
#survey_coding <- as.data.frame(vroom("../get_data/features_improved_3.csv", delim=","))
survey_coding <- as.data.frame(fread("../get_data/features_improved_3.csv", sep=","))


real_coefs <- grep("_", rownames(save_coef), value = T)
type_coef <- unlist(lapply(strsplit(real_coefs, "_"), function(x) x[[1]]))
translate_coef <- rep(" ", length(type_coef))
for(curr_type in unique(type_coef)){
  print(curr_type)
  if(curr_type == "SURVEY"){
    for(curr_coef in real_coefs[type_coef == curr_type]){
      new_coef <- strsplit(strsplit(curr_coef, "_X")[[1]][2], ".", fixed = T)[[1]][1]
      poss_names <- survey_coding[survey_coding[,1] == paste0(new_coef, "-0.0"), ]
      if(length(strsplit(curr_coef, ".", fixed = T)[[1]]) == 4){
        one_hot_ind <- as.numeric(strsplit(curr_coef, ".", fixed = T)[[1]][4]) + 1
        final_name <- paste(poss_names[one_hot_ind,3:4], collapse = " ")
      } else {
        final_name <- poss_names[1,3]
      }
      translate_coef[which(real_coefs == curr_coef)] <- final_name
    }

  } else if(curr_type == "CENSUS"){
    for(curr_coef in real_coefs[type_coef == curr_type]){
      translate_coef[which(real_coefs == curr_coef)] <- strsplit(curr_coef, "_")[[1]][2]
    }

  } else {
    print("curr_type")
    print(curr_type)
    coding_dict <- as.data.frame(fread(paste0("~/athena/prep_uk_data/", coding_file[data_type == curr_type]), sep="\t"))
    for(curr_coef in real_coefs[type_coef == curr_type]){
      translate_coef[which(real_coefs == curr_coef)] <- coding_dict[coding_dict[,1] == strsplit(curr_coef, "_")[[1]][2],2][1]
    }

  }
}

print("length-2")
print(length(translate_coef))
save_coef$code_name <- rownames(save_coef)
save_coef$english_name <- NA
save_coef$english_name[grepl("_", rownames(save_coef))] <- translate_coef

saveRDS(save_coef, paste0("reg_res/", author, ".", comp,  ".coefs.RDS"))

#}

print("one")
run_regressions(hi_prs_no_disease_eid, hi_prs_yes_disease_eid, "disease_in_hi_prs")
print("two")
run_regressions(lo_prs_no_disease_eid, lo_prs_yes_disease_eid, "disease_in_lo_prs")
print("three")
run_regressions(hi_prs_no_disease_eid, lo_prs_yes_disease_eid, "hi_lo_all")
#print("four")
#run_regressions(lo_prs_no_disease_eid, hi_prs_yes_disease_eid, "hi_lo_all")

#What do I actually want out of this process?
#Check to see if accuracy improves? - don't really care and this would mean overfitting
#Just want features that are significant, form a forest plot from the results
#data viz does not need to be fancy
