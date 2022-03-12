library(glmnet)
library(survival)
library(Hmisc)

#author <- commandArgs(trailingOnly=TRUE)[1]
author <- "christophersen"

#train_data <- readRDS(paste0("../batch_indivs/data/use_big_data.train.", author, ".RDS"))
#test_data <- readRDS(paste0("../batch_indivs/data/use_big_data.test.", author, ".RDS"))

#train_timeline <- readRDS("../comp_outcomes/timeline_res/all_eid.time.RDS")
#test_timeline <- readRDS("../comp_outcomes/timeline_res/test_eid.time.RDS")

surv_data <- readRDS(paste0("../batch_indivs/surv_data/survdf.", author, ".RDS"))
surv_data$age <- (as.Date(surv_data$final_rec_date) - as.Date(surv_data$dob))/365

#final_rec_date is just the end_date, either date of death, date of disease diagnosis, date lost to follow up, or date of censorship, whichever comes first
#should do this right, read in big data, find bad people and cols, remove minor counts, correlates, split by external eid lists, then impute




###########################################################
#    Get the known incidence rate
############################################################
#the rate, val is the column head, is number of individuals per 100k

translator <- read.table("~/athena/pgs_working_group/absrisk/author_translator.txt", stringsAsFactors=F, sep="\t")
gen_inci <- read.table("disease_prev/gbd_inci.csv", stringsAsFactors=F, header=T, sep=",")
gen_inci <- gen_inci[gen_inci$cause == translator[translator[,2] == author,1] & gen_inci$sex == "Both",]

#exit()
gen_prev <- read.table("disease_prev/gbd_prev.csv", stringsAsFactors=F, header=T, sep=",")
gen_prev <- gen_prev[gen_prev$cause == translator[translator[,2] == author,1] & gen_prev$sex == "Both",]
#let's try to get these to line up
#gen_inci <- gen_inci[gen_inci$year == 2019,c(4,8)]
#gen_prev <- gen_prev[gen_prev$year == 2019,c(4,8)]
#gen_prev$inci <- gen_inci$val
#gen_prev$diff_prev <- c(0, (c(gen_prev$val[-1], gen_prev$val[nrow(gen_prev)]) - gen_prev$val)[-nrow(gen_prev)])


#WAIT! Is the incidence all wrong if it is stated to occur over 5 year ages
#Incidence = New cases per 100,000 population
#I might think these numbers need to be divided by 5 to spread it out over the years


u_year <- sort(unique(gen_inci$year))

gbd_inci <- list()
gbd_prev <- list()
for(i in 1:length(u_year)){
  #inci
  sub_inci <- gen_inci[gen_inci$year == u_year[i],]
  sub_inci$age <- seq(2, 92, by = 5)
  sub_inci <- sub_inci[sub_inci$age > 37,]

  new_inci <- data.frame("age" = min(sub_inci$age):max(sub_inci$age))
  
  funcy <- splinefun(sub_inci$age, sub_inci$val, method = "fmm")
  new_inci$val <- funcy(new_inci$age)

  funcy <- splinefun(sub_inci$age, sub_inci$upper, method = "fmm")
  new_inci$upper <- funcy(new_inci$age)

  funcy <- splinefun(sub_inci$age, sub_inci$lower, method = "fmm")
  new_inci$lower <- funcy(new_inci$age)

  new_inci[,2:4] <- new_inci[,2:4]/100000
  #new_inci[,2:4] <- new_inci[,2:4]/5 #need to adjust for a single age year
  gbd_inci[[i]] <- new_inci

  #prev
  sub_prev <- gen_prev[gen_prev$year == u_year[i],]
  sub_prev$age <- seq(2, 92, by = 5)
  sub_prev <- sub_prev[sub_prev$age > 37,]

  new_prev <- data.frame("age" = min(sub_prev$age):max(sub_prev$age))

  funcy <- splinefun(sub_prev$age, sub_prev$val, method = "fmm")
  new_prev$val <- funcy(new_prev$age)

  funcy <- splinefun(sub_prev$age, sub_prev$upper, method = "fmm")
  new_prev$upper <- funcy(new_prev$age)

  funcy <- splinefun(sub_prev$age, sub_prev$lower, method = "fmm")
  new_prev$lower <- funcy(new_prev$age)

  #pop <- read.table("disease_prev/IHME_GBD_2019_POP_2019_Y2020M10D15.CSV.gz", stringsAsFactors=F, header=T, sep=",")
  #pop <- pop[pop$location_name == "United Kingdom" & pop$age_group_name %in% gen_prev$age & pop$sex_name == "both",]

  new_prev[,2:4] <- new_prev[,2:4]/100000
  #new_prev[,2:4] <- new_prev[,2:4]/5 #need to adjust for a single age year
  gbd_prev[[i]] <- new_prev
}

names(gbd_inci) <- u_year
names(gbd_prev) <- u_year

#need to adjust prevalence for the total population

##################################################
#   Get the average empirical incidence rate
####################################################

#case_surv_data <- surv_data[surv_data$pheno == 1,]
#case_surv_data$end_date <- as.Date(as.character(case_surv_data$end_date))
#case_surv_data$start_date <- as.Date(as.character(case_surv_data$start_date))
#case_surv_data$age_at_diag <- round(as.numeric((case_surv_data$end_date -  case_surv_data$dob)/365))

#diag_years <- unlist(lapply(strsplit(as.character(case_surv_data$end_date), "-", fixed = T), function(x) x[1]))
#diag_fracs <- table(diag_years)/length(diag_years)

#comp_gbd_prev <- data.frame(matrix(0, nrow = nrow(gbd_prev[[1]]), ncol = 4))
#colnames(comp_gbd_prev) <-  colnames(gbd_prev[[1]])

#for(i in 1:length(diag_fracs)){
#  comp_gbd_prev <- comp_gbd_prev + gbd_prev[[which(names(gbd_prev) == names(diag_fracs)[1])]] * diag_fracs[i]
#}




#####################################################3
# Get chance each person has disease
####################################################

#Assume each person could have ICD starting at year 2000
#Remove individuals who are prevalent for disease


##################################################
#   Get the UKBB incidence rate
####################################################

surv_data$end_date <- as.Date(as.character(surv_data$end_date))
surv_data$start_date <- as.Date(as.character(surv_data$start_date))


all_years <- 2010:2019
ukb_inci <- list()
ukb_prev <- list()

#in year k
#take the current age
#

for(k in 1:length(all_years)){

  comp_ukb_inci <- data.frame(matrix(0, nrow = nrow(gbd_inci[[i]]), ncol = 5))
  colnames(comp_ukb_inci) <-  c(colnames(gbd_inci[[1]]), "n")
  comp_ukb_inci[,1] <- gbd_inci[[1]][,1]

  comp_ukb_prev <- data.frame(matrix(0, nrow = nrow(gbd_prev[[i]]), ncol = 5))
  colnames(comp_ukb_prev) <-  c(colnames(gbd_prev[[1]]), "n")
  comp_ukb_prev[,1] <- gbd_prev[[1]][,1]


  #current age is measured from July 1 - with floor function
  surv_data$current_age <- floor(as.numeric(as.Date(paste0(all_years[k], "-07-01")) - surv_data$dob)/365)
  #surv_data$current_age <- all_years[k] - as.numeric(substr(surv_data$dob, 1, 4))
  for(j in 1:nrow(gbd_inci[[1]])){
    age_surv_data <- surv_data[surv_data$current_age == gbd_inci[[1]][j,1],]
    #have to lived more than half the year
    bad_eid <- age_surv_data$eid[age_surv_data$event_type == "death"][age_surv_data$end_date[age_surv_data$event_type == "death"] > as.Date(paste0(all_years[k], "-07-01"))]
    age_surv_data <- age_surv_data[!(age_surv_data$eid %in% bad_eid),]

    if(nrow(age_surv_data) > 0){
      #this should be the number who are currently this exact age and have ever been diagnosed
      #comp_ukb_prev[j,2] <- sum(age_surv_data$end_date[age_surv_data$pheno == 1] < as.Date(paste0(all_years[k], "-06-15")))/nrow(age_surv_data)
      comp_ukb_inci[j,2] <- sum(substr(age_surv_data$end_date[age_surv_data$pheno == 1], 1, 4) == all_years[k])/nrow(age_surv_data)
      comp_ukb_inci[j,5] <- nrow(age_surv_data)

      comp_ukb_prev[j,2] <- sum(substr(age_surv_data$end_date[age_surv_data$pheno == 1], 1, 4) <= all_years[k])/nrow(age_surv_data)
      comp_ukb_prev[j,5] <- nrow(age_surv_data)

      if(comp_ukb_inci[j,2] > 0){
        pro_se <- sqrt((comp_ukb_inci[j,2] * (1 - comp_ukb_inci[j,2]))/nrow(age_surv_data))
        comp_ukb_inci[j,3] <- (comp_ukb_inci[j,2] + pro_se)
        comp_ukb_inci[j,4] <- (comp_ukb_inci[j,2] - pro_se) 
        comp_ukb_inci[j,2] <- comp_ukb_inci[j,2]
      }

      if(comp_ukb_prev[j,2] > 0){
        pro_se <- sqrt((comp_ukb_prev[j,2] * (1 - comp_ukb_prev[j,2]))/nrow(age_surv_data))
        comp_ukb_prev[j,3] <- (comp_ukb_prev[j,2] + pro_se)
        comp_ukb_prev[j,4] <- (comp_ukb_prev[j,2] - pro_se)
        comp_ukb_prev[j,2] <- comp_ukb_prev[j,2]
      }

    }
  }

  ukb_inci[[k]] <- comp_ukb_inci
  ukb_prev[[k]] <- comp_ukb_prev
}

names(ukb_inci) <- all_years
names(ukb_prev) <- all_years

#exit()

#now I need an age range for each year to determine which varieties of groups I am really looking at
#maybe for ukbb also make a column of sample size
#dimensions of all_missing are age by year
all_missing_inci <- data.frame(matrix(0, nrow = nrow(ukb_inci[[1]]), ncol = length(ukb_inci)+1 ))
all_missing_inci[,1] <- ukb_inci[[1]]$age
colnames(all_missing_inci) <- c("age", paste0("year_", names(ukb_inci)))

all_missing_prev <- data.frame(matrix(0, nrow = nrow(ukb_prev[[1]]), ncol = length(ukb_prev)+1 ))
all_missing_prev[,1] <- ukb_prev[[1]]$age
colnames(all_missing_prev) <- c("age", paste0("year_", names(ukb_prev)))



for(i in 1:length(ukb_inci)){
  curr_uk <- ukb_inci[[i]][ukb_inci[[i]]$n > 1000,]
  curr_gbd <- gbd_inci[[which(names(gbd_inci) == names(ukb_inci)[i])]][gbd_inci[[which(names(gbd_inci) == names(ukb_inci)[i])]]$age %in% curr_uk$age,]

  all_missing_inci[all_missing_inci$age %in% curr_uk$age,i+1] <- (curr_gbd$val - curr_uk$val)

  curr_uk <- ukb_prev[[i]][ukb_prev[[i]]$n > 1000,]
  curr_gbd <- gbd_prev[[which(names(gbd_prev) == names(ukb_prev)[i])]][gbd_prev[[which(names(gbd_prev) == names(ukb_prev)[i])]]$age %in% curr_uk$age,]

  all_missing_prev[all_missing_prev$age %in% curr_uk$age,i+1] <- (curr_gbd$val - curr_uk$val)
}


#the comparison between incidence and prevalence actually looks good, although not perfect
#the problem must be how I calculate the UKBB incidence
#it is more instantaneous and therefore more sensitive to randomness (bad cutoffs?)
#so just stick with prevalence here and incidence in the absrisk calculation



#comp_ukb_prev <- data.frame(matrix(0, nrow = nrow(gbd_prev[[i]]), ncol = 4))
#colnames(comp_ukb_prev) <-  colnames(gbd_prev[[1]])
#comp_ukb_prev[,1] <- gbd_prev[[1]][,1]

#for(j in 1:nrow(comp_ukb_prev)){
#  samp_size <- sum(surv_data$age >= comp_ukb_prev[j,1])
#  sub_surv_data <- case_surv_data[case_surv_data$age >= comp_ukb_prev[j,1],]
#  if(nrow(sub_surv_data) > 0){
#    #this should be the number who are currently this exact age and have ever been diagnosed
#    comp_ukb_prev[j,2] <- (sum(sub_surv_data$age_at_diag <= comp_ukb_prev[j,1]) / samp_size ) 
#
#    if(comp_ukb_prev[j,2] > 0){
#      pro_se <- sqrt((comp_ukb_prev[j,2] * (1 - comp_ukb_prev[j,2]))/samp_size)
#      comp_ukb_prev[j,3] <- (comp_ukb_prev[j,2] + pro_se) * 100000
#      comp_ukb_prev[j,4] <- (comp_ukb_prev[j,2] - pro_se) * 100000
#      comp_ukb_prev[j,2] <- comp_ukb_prev[j,2] * 100000
#    }
#  }
#}


#comp_gbd_prev <- comp_gbd_prev[22:54,]
#comp_ukb_prev <- comp_ukb_prev[22:54,]




###############################################
#   Make models
#################################################

#want to make a predictive model for each year after assessment (assume survey questions stay the same for 10 years)
#train on individuals who have an ICD diagnosis past the age denoted
#test on individuals who do not have any ICD diagnosis paste the age denoted
#note that all available cases will be in training set

#say I'm looking to make a model for age 60
#split into train and test, case and control
#I am changing the case and control based on the age, so shouldnt I also change the feature set
#ICD codes that have not occured until a time should be disregarded
#can assume that other non-ICD codes are constant throughout time (or at least from 40 onwards)

#currently set up to make model at one timepoint
#will still have to use multiple prev values


#iterating over each time
# step 1 - calculate how many people we know are diagnosed
# step 2 - calculate how many more people are needed (at a few thresholds)
# step 3 - identify people who may fall into this diagnosis group
# step 4 - make an enet model, fit on people we do know, assessing those we do not
# step 5 - with predictions, identify groups of people who are either likely with disease (at a few threshold) or we do not know


#across all times find consensus groups of people with disease or need to be removed?
#incorparate those people into a survival model rather than just a binary model

#ensure that only individuals we know are in the training group
#specifically, 80% of UKBB is subset up top
#then, X% who have recent ICD diagnosis become train, remaining are test
#on the X% do CV to determine best alpha, predict with best alpha on test
#Repeat last 2 steps but the train group is random 


####
#Problem is that because of volunteer bias UKBB is healthier than average, so they will have less prevalence than the larger populace
#how will I know what the true prevalence for my population should be?
#can try to scale the prevalence based on the self-reported health
#can forget the GBD prevalence and simply compare prevalence of those with and without recent EHR
#but individuals with limited EHR will likely be healthier and have less prevalence than the group that has had a recent EHR event
#the cox model predictions should take into account whether someone is at low or high risk on a scale relative to the population


big_data <- readRDS(paste0("use_big_data.", author, ".RDS"))
eid_big_data <- readRDS("eid_big_data.RDS")


#first a random group
train_surv <- surv_data[1:100000,]
train_data <- big_data[eid_big_data %in% train_surv$eid,]
train_eid <- eid_big_data[eid_big_data %in% train_surv$eid]
train_surv <- train_surv[train_surv$eid %in% train_eid,]

train_data <- train_data[order(train_eid)[rank(train_surv$eid)],]

test_surv <- surv_data[!(surv_data$eid %in% train_surv$eid),]
test_data <- big_data[eid_big_data %in% test_surv$eid,]
test_eid <- eid_big_data[eid_big_data %in% test_surv$eid]
test_surv <- test_surv[test_surv$eid %in% test_eid,]

test_data <- test_data[order(test_eid)[rank(test_surv$eid)],]



y_data <- train_surv[,1:2]
colnames(y_data) <- c("time", "status")

write.table(train_data, "x_data", row.names = F, quote = F, col.names = T, sep = ",")
write.table(y_data, "y_data", row.names = F, quote = F, col.names = T, sep = ",")

#expect to return a set of coefficients that can be used to predict the test data
system("python enet_cox.py")
system("rm x_data y_data")

best_coef <- read.csv("best_coef.csv", stringsAsFactors=F, header=T)

train_preds <- rowSums(t(t(as.matrix(train_data)) * best_coef[,2]))
test_preds <- rowSums(t(t(as.matrix(test_data)) * best_coef[,2]))


exit()





#Then a group of people we know
all_timeline <- readRDS("../comp_outcomes/timeline_res/bentham.new_all_eid.time.RDS")
all_timeline$last_diag <- as.Date(all_timeline$last_diag, "1970-01-01")
all_timeline$last_icd_diag <- as.Date(all_timeline$last_icd_diag, "1970-01-01")
all_timeline$last_oper_diag <- as.Date(all_timeline$last_oper_diag, "1970-01-01")
all_timeline$first_diag <- as.Date(all_timeline$first_diag, "1970-01-01")
all_timeline$first_icd_diag <- as.Date(all_timeline$first_icd_diag, "1970-01-01")
all_timeline$first_oper_diag <- as.Date(all_timeline$first_oper_diag, "1970-01-01")
















#------------------------------------------------------
# Enet model over all features
#-----------------------------------------------------------



cvmod <- cv.glmnet(y = train_df$pheno, x = as.matrix(big_data), alpha = 1, family = "binomial")
#mod <- glmnet(y = train_df$pheno, x = as.matrix(big_data), alpha = 1, lambda = cvmod$lambda.min, family = "binomial")
mod <- glmnet(y = train_df$pheno, x = as.matrix(big_data), alpha = 1, lambda = 0.0001, family = "binomial")

enet_preds <- predict(mod, as.matrix(big_data))
test_preds <- predict(mod, as.matrix(test_data))




wrong_list <- list(train_df$eid[enet_preds > quantile(enet_preds, 0.99)],
                   train_df$eid[enet_preds > quantile(enet_preds, 0.995)],
                   train_df$eid[enet_preds > quantile(enet_preds, 0.9995)],
                   train_df$eid[enet_preds > quantile(enet_preds, 0.9999)])

test_wrong_list <- list(valid_df$eid[test_preds > quantile(enet_preds, 0.99)],
                        valid_df$eid[test_preds > quantile(enet_preds, 0.995)],
                        valid_df$eid[test_preds > quantile(enet_preds, 0.9995)],
                        valid_df$eid[test_preds > quantile(enet_preds, 0.9999)])

control_wrong_list <- list(train_df$eid[enet_preds < quantile(enet_preds, 0.01)],
                           train_df$eid[enet_preds < quantile(enet_preds, 0.005)],
                           train_df$eid[enet_preds < quantile(enet_preds, 0.0005)],
                           train_df$eid[enet_preds < quantile(enet_preds, 0.0001)])

control_test_wrong_list <- list(valid_df$eid[test_preds < quantile(enet_preds, 0.01)],
                                valid_df$eid[test_preds < quantile(enet_preds, 0.005)],
                                valid_df$eid[test_preds < quantile(enet_preds, 0.0005)],
                                valid_df$eid[test_preds < quantile(enet_preds, 0.0001)])


saveRDS(mod, paste0("out_batch/wrong.model.", author, ".RDS"))

saveRDS(wrong_list, paste0("out_batch/case_groups.train.", author, ".RDS"))
saveRDS(test_wrong_list, paste0("out_batch/case_groups.test.", author, ".RDS"))

saveRDS(control_wrong_list, paste0("out_batch/control_groups.train.", author, ".RDS"))
saveRDS(control_test_wrong_list, paste0("out_batch/control_groups.test.", author, ".RDS"))
