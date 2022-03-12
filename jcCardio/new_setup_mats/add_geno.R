library(data.table)

af_df <- readRDS("vt_pheno_mat.RDS")


opt = "clinvar_loftee"
#opt = "loftee_all"

j <- 1
all_geno <- list()
for(i in 1:22){
  print(i)
  if(file.exists(paste0("../get_geno/curr.", i, ".raw.gz"))){
    if(opt == "clinvar" | opt == "loftee_all"){
      addon = "/all_"
      otheron = "_all_"
    } else {
      addon = "/"
      otheron = "_"
    }

    geno <- as.data.frame(fread(paste0("../get_geno", addon, "curr.", i, ".raw.gz"), stringsAsFactors=F, header = T))
    geno <- geno[geno[,1] %in% af_df$eid,]
    geno <- geno[order(geno[,1])[rank(af_df$eid)],]
    geno <- geno[,7:ncol(geno), drop=F]

    genes <- read.table(paste0("../get_geno", addon, "gene.", i), stringsAsFactors=F)
    id <- read.table(paste0("../get_geno", addon, "id.", i), stringsAsFactors=F)
    clin_var <- as.data.frame(fread(paste0("../get_geno/clin", otheron, "data.", i), header = F, sep = "\t"))
    clin_cldn <- as.data.frame(fread(paste0("../get_geno/clin", otheron, "cldn.", i), header = F, sep = "\t"))

    if(nrow(clin_var) > ncol(geno)){
      clin_id <- read.table(paste0("../get_geno/clin_id.", i), stringsAsFactors=F)
      geno_id <- unlist(lapply(strsplit(colnames(geno), "_"), function(x) x[1]))

      clin_var <- clin_var[clin_id[,1] %in% geno_id, 1, drop = F]
      clin_cldn <- clin_cldn[clin_id[,1] %in% geno_id, 1, drop = F]
    }

    if(nrow(clin_var) < ncol(geno)){
      print("ERROR")
    }


    if(sum(duplicated(id[,1])) > 0){
      genes <- genes[!duplicated(id[,1]),,drop=F]
    }


    #subset the total genotype, either all of the rare variants in the genes or the loftee HC variants by clinvar
    if(opt == "clinvar_loftee" | opt == "clinvar"){
      clin_bool <- clin_var[,1] == "Pathogenic" | clin_var[,1] == "Likely_pathogenic" | clin_var[,1] == "Pathogenic/Likely_pathogenic"
      clin_bool[is.na(clin_bool)] <- F
      if(any(clin_bool)){
        geno <- geno[,clin_bool,drop=F]
        genes <- genes[clin_bool,,drop = F]
        #id <- id[clin_bool,,drop=F]

        colnames(geno) <- paste0(colnames(geno), "_",  genes[,1])
        all_geno[[j]] <- geno
        j <- j + 1
      } 
 
    } else if(opt == "loftee"){
      colnames(geno) <- paste0(colnames(geno), "_",  genes[,1])
      all_geno[[j]] <-geno
      j <- j + 1


    } else if(opt == "loftee_all"){
      hc <- read.table(paste0("../get_geno/loftee_hc_all_data.", i), stringsAsFactors=F)

      geno_pos <- unlist(lapply(strsplit(colnames(geno), ":", fixed = T), function(x) x[2]))
      hc_pos <- unlist(lapply(strsplit(hc[,1], ":", fixed = T), function(x) x[2]))
      genes <- genes[geno_pos %in% hc_pos,,drop=F]
      geno <- geno[,geno_pos %in% hc_pos]
      colnames(geno) <- paste0(colnames(geno), "_",  genes[,1])
      all_geno[[j]] <-geno
      j <- j + 1
    }





  }
}


geno <- do.call("cbind", all_geno)
colnames(geno) <- paste0("X1_", colnames(geno))

for(f in list.files(pattern = "pheno_mat")){
  pheno <- readRDS(f)
  pheno <- cbind(pheno, geno)
  
  header <- strsplit(f, "_", fixed=T)[[1]][1]
  saveRDS(pheno, paste0(header, "_ready.", opt, ".RDS"))

}


