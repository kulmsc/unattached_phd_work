library(SKAT)
library(data.table)

df <- readRDS("vt_pheno_mat.RDS")
df$waist[is.na(df$waist)] <- mean(df$waist, na.rm = T)
df$hip_ratio[is.na(df$hip_ratio)] <- mean(df$hip_ratio, na.rm = T)
df$BMI[is.na(df$BMI)] <- mean(df$BMI, na.rm = T)
df$tobacco[is.na(df$tobacco)] <- 0


list_burden_pval <- list()
list_single_pval <- list()
i <- 1
for(f in grep("all", list.files("../get_geno/", "raw.gz"), value = T)){
  print(f)
  geno <- as.data.frame(fread(paste0("../get_geno/", f), stringsAsFactors=F, header = T))
  geno <- geno[geno[,1] %in% df$eid,]
  geno <- geno[order(geno[,1])[rank(df$eid)],]
  geno <- geno[,7:ncol(geno), drop=F]



  genes <- read.table(paste0("../get_geno/all_gene.", strsplit(f, ".", fixed = T)[[1]][2]), stringsAsFactors=F)
  id <- read.table(paste0("../get_geno/all_id.", strsplit(f, ".", fixed = T)[[1]][2]), stringsAsFactors=F)
  if(sum(duplicated(id[,1])) > 0){
    genes <- genes[!duplicated(id[,1]),,drop=F]
  }
  use_genes <- names(table(genes)[table(genes) > 10])



  #null_mod_inter <- SKAT_Null_Model(df$pheno ~ 1, out_type = "D")
  null_mod_other <- with(df, SKAT_Null_Model(pheno ~ age + sex + brit + euro + asian + african + BMI + waist + hip_ratio  , out_type = "D"))
  #null_mod_comor <- with(df, SKAT_Null_Model(pheno ~ age + sex + brit + euro + asian + african + BMI + waist + hip_ratio + asthma +  cad + ckd + copd + depress + dm + hld + htn +hyper + hypo + hf + stroke + nic + osa + pvd + icm + tobacco + freq_etoh, out_type = "D"))

   all_burden_pvals <- data.frame(matrix(0, nrow = length(use_genes), ncol = 3))
   colnames(all_burden_pvals) <- c("inter", "other", "comor")
   for(j in 1:length(use_genes)){
     print(use_genes[j])
     subgeno <- as.matrix(geno[,genes[,1] == use_genes[j]])
     subgeno <- subgeno[,colSums(subgeno, na.rm = T) > 10]
     #all_burden_pvals[j,1] <- SKATBinary_Robust(subgeno, null_mod_inter, max_maf = 0.05)$p.value
     all_burden_pvals[j,2] <- SKATBinary_Robust(subgeno, null_mod_other, max_maf = 0.05)$p.value
     #all_burden_pvals[j,3] <- SKATBinary_Robust(subgeno, null_mod_comor, max_maf = 0.05)$p.value
   }
   all_burden_pvals$gene <- use_genes


   subgeno <- geno[, (colSums(geno, na.rm = T) / nrow(geno)*2) > 0.005]
   subgenenames <- genes[ (colSums(geno, na.rm = T) / nrow(geno)*2) > 0.005,]
   all_single_pval <- rep(NA, ncol(subgeno))
   for(k in 1:ncol(subgeno)){
     df$gene <- subgeno[,k]
     mod <- glm(pheno ~ age + sex + brit + euro + asian + african + BMI + waist + hip_ratio + gene, data = df, family = "binomial")
     all_single_pval[k] <- summary(mod)$coef[nrow(summary(mod)$coef), 4]
   }
   names(all_single_pval) <- colnames(subgeno)
   names(all_single_pval) <- paste0(names(all_single_pval), "_", subgenenames)


  list_burden_pval[[i]] <- all_burden_pvals
  list_single_pval[[i]] <- all_single_pval
  i <- i + 1
}


final_burden_pval <- do.call("rbind", list_burden_pval)
final_single_pval <- unlist(list_single_pval)


saveRDS(list(final_burden_pval, final_single_pval), "results/pvc_all_pval.RDS")
