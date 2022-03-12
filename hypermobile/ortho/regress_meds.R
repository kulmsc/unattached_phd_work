
meds <- readRDS("meds.RDS")
meds <- meds[,colSums(meds) > 50]

covars <- read.table("use_covar", stringsAsFactors=F, header=T)
meds <- meds[meds$eid %in% covars[,1],]
meds <- meds[order(meds$eid)[rank(covars[,1])],]

disease_eid <- read.table("eid_gwas_3", stringsAsFactors=F)

covars$pheno <- 0
covars$pheno[covars[,1] %in% disease_eid[,1]] <- 1
covars <- covars[,3:ncol(covars)]
meds <- meds[,-ncol(meds)]

all_fish_p <- rep(NA, ncol(meds))
all_reg_p <- rep(NA, ncol(meds))

for(i in 1:ncol(meds)){
  
  mat <- c(sum(meds[,i] == 1 & covars$pheno == 1),
           sum(meds[,i] == 1 & covars$pheno == 0),
           sum(meds[,i] == 0 & covars$pheno == 1),
           sum(meds[,i] == 0 & covars$pheno == 0))
  all_fish_p[i] <- fisher.test(matrix(mat, nrow = 2))$p.value

  covars$med <-meds[,i]
  mod <- glm(pheno ~ ., data = covars, family = "binomial")
  all_reg_p[i] <- summary(mod)$coef[nrow(summary(mod)$coef), 4]
}
