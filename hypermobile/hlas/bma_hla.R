library(stringr)
library(BMA)

args = commandArgs(trailingOnly=TRUE)
index <- args[1]
cov_type <- args[2]
ethnic <- args[3]

hla <- read.table("~/athena/ukbiobank/hla/ukb26867.csv.gz", sep=',', header = T, stringsAsFactors = F)
hla_desc <- read.table("../common_files/ukb_hla_v2.txt", stringsAsFactors=F)

total_phen <- read.table("use_pheno", stringsAsFactors=F, header = F)
total_covar <- read.table("use_covar", stringsAsFactors=F)
total_sex <- read.table("use_fam", stringsAsFactors=F)
total_df <- cbind(pheno = total_phen[,3], total_covar[,-1], sex = total_sex[,5])

total_hla_val <- matrix(0, nrow = nrow(total_df), ncol = 362)
for(i in 1:nrow(total_df)){
  curr_row <- which(hla[,1] == total_sex[i,1])
  total_hla_val[i,] <- as.numeric(str_split(hla[curr_row,2], ",", simplify=T))
}

res <- read.table(paste0("september_results/hla_gwas_res_", index, "_", cov_type, "_", ethnic, ".txt"), stringsAsFactors=F)
best_names <- res[res[,5] < 0.05, 1]
if(length(best_names) > 10){
  best_names <- res[res[,5] < sort(res[,5])[10], 1]
}

if(length(best_names) > 2){

  sig_hla_val <- total_hla_val[,hla_desc %in% best_names]
  bayes_mod <- bic.glm(y = total_phen[,3]-1, x = sig_hla_val, glm.family = "binomial")
  post <- bayes_mod$postmean
  names(post) <- c("intercept", hla_desc[hla_desc %in% best_names])

  write.table(post, paste0("september_results/bayes_res_", index, "_", cov_type, "_", ethnic, ".txt"), 
  row.names = T, col.names = F, quote = F, sep = '\t')

} else {
  write.table("empty", paste0("september_results/bayes_res_", index, "_", cov_type, "_", ethnic, ".txt"),
  row.names = T, col.names = F, quote = F, sep = '\t')

}
