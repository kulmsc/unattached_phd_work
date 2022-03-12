library(vroom)
library(data.table)

author <- "onengut"
#author <- commandArgs(trailingOnly=TRUE)[1]

############################# READ IN AND SET UP #########################################

train_df <- readRDS(paste0("../get_data/df_output/df_train.", author, ".RDS"))
test_df <- readRDS(paste0("../get_data/df_output/df_test.", author, ".RDS"))
surv_train_df <- readRDS(paste0("../get_data/df_output/survdf_train.", author, ".RDS"))
surv_test_df <- readRDS(paste0("../get_data/df_output/survdf_test.", author, ".RDS"))

big_data <- readRDS("../get_data/big_data.RDS")
colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))] <- paste0("X", colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))])

night <- read.csv("~/athena/ukbiobank/setup_sixthPhenos/ukb46898.csv.gz", stringsAsFactors=F, header=T)
eid_night <- night$eid
night <- night[,14:348]
night <- night[,grepl(".0.0", colnames(night))]
eid_night <- eid_night[!is.na(night[,1])]
night <- night[!is.na(night[,1]),]


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

run_regressions <- function(ppl_1, ppl_2, comp){
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

  names(test_p) <- colnames(small_data)
  names(mod_coef) <- colnames(small_data)
  names(mod_p) <- colnames(small_data)
  save_coef <- list("test_p" = test_p, "mod_coef" = mod_coef, "mod_p" = mod_p)

  saveRDS(save_coef, paste0("reg_res/", author, ".", comp,  ".coefs.RDS"))

}

check_night <- function(ppl1, ppl2, comp){
  use_night <- rbind(night[eid_night %in% ppl1,], night[eid_night %in% ppl2,])
  pheno_night <- c(rep(0, sum(eid_night %in% ppl1)), rep(1, sum(eid_night %in% ppl2)))

  p_glm_night <- rep(NA, ncol(night))
  p_wilcox_night <- rep(NA, ncol(night))
  for(i in 1:ncol(use_night)){
    mod <- glm(pheno_night ~ use_night[,i], family = "binomial")
    p_glm_night[i] <- summary(mod)$coef[2,4]
    p_wilcox_night[i] <- wilcox.test(use_night[pheno_night == 0,i], use_night[pheno_night == 1,i], simulate.p = TRUE)$p.value
  }

  stats <- data.frame("cols" = colnames(night), "wilcox" = p_wilcox_night, "glm" = p_glm_night)
  saveRDS(stats, paste0("reg_res/", author, ".", comp,  ".night.RDS"))
}


exit()
print("one")
#run_regressions(hi_prs_no_disease_eid, hi_prs_yes_disease_eid, "disease_in_hi_prs")
check_night(hi_prs_no_disease_eid, hi_prs_yes_disease_eid, "disease_in_hi_prs")
print("two")
#run_regressions(lo_prs_no_disease_eid, lo_prs_yes_disease_eid, "disease_in_lo_prs")
check_night(lo_prs_no_disease_eid, lo_prs_yes_disease_eid, "disease_in_lo_prs")
print("three")
#run_regressions(hi_prs_no_disease_eid, lo_prs_yes_disease_eid, "hi_lo_all")
check_night(hi_prs_no_disease_eid, lo_prs_yes_disease_eid, "hi_lo_all")
#print("four")
#run_regressions(lo_prs_no_disease_eid, hi_prs_yes_disease_eid, "hi_lo_all")

#What do I actually want out of this process?
#Check to see if accuracy improves? - don't really care and this would mean overfitting
#Just want features that are significant, form a forest plot from the results
#data viz does not need to be fancy
