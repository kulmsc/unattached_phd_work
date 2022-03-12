library(stringr)

#args = commandArgs(trailingOnly=TRUE)
#index <- args[1]
index = 3

hla <- read.table("~/athena/ukbiobank/hla/ukb26867.csv.gz", sep=',', header = T, stringsAsFactors = F)
hla_desc <- read.table("~/athena/COVID/common_files/ukb_hla_v2.txt", stringsAsFactors=F)

total_phen <- read.table("use_pheno", stringsAsFactors=F, header = F)
total_covar <- read.table("~/athena/xwas/gwas/covars", stringsAsFactors=F, header = T)
total_covar <- total_covar[total_covar[,1] %in% total_phen[,1],]
total_phen <- total_phen[total_phen[,1] %in% total_covar[,1],]
total_phen <- total_phen[order(total_phen[,1])[rank(total_covar[,1])],]

total_df <- cbind(pheno = total_phen[,3], total_covar[,-c(1,2)])

total_hla_val <- matrix(0, nrow = nrow(total_df), ncol = 362)
for(i in 1:nrow(total_df)){
  curr_row <- which(hla[,1] == total_phen[i,1])
  total_hla_val[i,] <- as.numeric(str_split(hla[curr_row,2], ",", simplify=T))
}

total_df <- total_df[!is.na(total_hla_val[,1]),]
total_phen <- total_phen[!is.na(total_hla_val[,1]),]
total_hla_val <- total_hla_val[!is.na(total_hla_val[,1]),]

#uniq_vals <- apply(total_hla_val, 2, function(x) length(unique(x)))
nonzero_vals <- apply(total_hla_val, 2, function(x) sum(x!=0))
#total_hla_val <- total_hla_val[, uniq_vals > 10]
total_hla_val <- total_hla_val[,nonzero_vals/nrow(total_hla_val) > 0.01]
#hla_desc <- hla_desc[,uniq_vals > 10, drop = F]
hla_desc <- hla_desc[,nonzero_vals/nrow(total_hla_val) > 0.01]

gwas_res <- matrix(0, nrow = ncol(total_hla_val), ncol = 4)
i <- 1
for(i in 1:nrow(gwas_res)){
  total_df$geno <- total_hla_val[,i]
  if(sum(total_df$geno[total_df$pheno == 2] != 0) < 5){
    gwas_res[i,] <- rep(NA, 4)
  } else {
    mod <- glm((pheno-1) ~ ., family = "binomial", data = total_df)
    gwas_res[i,] <- unname(summary(mod)$coefficients[nrow(summary(mod)$coefficients),])
  }
}

rownames(gwas_res) <- hla_desc[1,]
write.table(gwas_res, paste0("hla_gwas_res_", index, ".txt"), sep = '\t', row.names = T, col.names = F)



