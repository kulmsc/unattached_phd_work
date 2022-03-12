library(vroom)
library(stringr)
library(bigsnpr)

options(warn=2)

#impute <- as.data.frame(vroom("common_files/impute_rsids", col_names = F))
impute <- read.table("common_files/impute_rsids", header = F, stringsAsFactors = F)
colnames(impute) <- c("LOC", "SNP", "POS", "A1", "A2", "MAF", "AX", "INFO")

#First remove SNPs where the base is longer than one allele
#And make sure the SNP is either A, C, G, T
#Remove ambigous SNPs, those where the bases are actually pairs
impute <- impute[nchar(impute$A1) == 1 & nchar(impute$A2) == 1,]
impute <- impute[impute$A1 %in% c("A", "C", "G", "T") & impute$A2 %in% c("A", "C", "G", "T"),]
impute <- impute[!((impute$A1 == "A" & impute$A2 == "T") |
             (impute$A1 == "T" & impute$A2 == "A") |
             (impute$A1 == "G" & impute$A2 == "C") |
             (impute$A1 == "C" & impute$A2 == "G")),]


adjust_ss <- function(author){

  #Start by reading in the summary statistics file and a parameters file that states the column
  #names of the desired columns
  #ss <- as.data.frame(vroom(paste0(author, "/raw_", tolower(author), ".ss.gz")))
  ss <- read.table(paste0(author, "/raw_", tolower(author), ".ss.gz"), stringsAsFactors=F, header = T)
  params <- read.table(paste0(author, "/parameters"), stringsAsFactors=F)

  system(paste0("echo clean notes > ", author, "/clean.log"))
  system(paste0("echo after read: ", nrow(ss), " >> ", author, "/clean.log"))

  #Check if there are missing columns, and fill them in just so things don't go haywire below
  missing_cols <- params[!(params[,1] %in% colnames(ss)),1]
  if(length(missing_cols) > 0){
    if("NO_CHR" %in% missing_cols){
      ss$NO_CHR <- 0
    }
    if("NO_POS" %in% missing_cols){
      ss$NO_POS <- 0
    }
    if("NO_SE" %in% missing_cols){
      ss$NO_SE <- 0
    }
  }


  #Construct a common ss object, with common column names based on the manually curated params
  ss <-  ss[,c(which(colnames(ss) == params[1,1]),
               which(colnames(ss) == params[2,1]),
               which(colnames(ss) == params[3,1]),
               which(colnames(ss) == params[4,1]),
               which(colnames(ss) == params[5,1]),
               which(colnames(ss) == params[6,1]),
               which(colnames(ss) == params[7,1]),
               which(colnames(ss) == params[8,1]))]
  colnames(ss) <- c("CHR", "BP", "RSID", "A1", "A2", "SE", "BETA", "P")
  ss$A1 <- toupper(ss$A1)
  ss$A2 <- toupper(ss$A2)
  ss <- ss[!is.na(ss$RSID),]
  ss <- ss[!is.na(ss$BETA),]
  ss <- ss[!is.na(ss$P),]

  system(paste0("echo after remove NA: ", nrow(ss), " >> ", author, "/clean.log"))

  #Remove duplicate SNP IDs
  dup_ids <- ss$RSID[duplicated(ss$RSID)]
  dup_ids <- c(dup_ids, duplicated(ss[,c('CHR','BP')]) | duplicated(ss[,c('CHR','BP')],fromLast = TRUE))
  ss <- ss[!(ss$RSID %in% dup_ids),]

  system(paste0("echo after remove duplicated IDs: ", nrow(ss), " >> ", author, "/clean.log"))

  #First remove SNPs where the base is longer than one allele
  #And make sure the SNP is either A, C, G, T
  #Remove ambigous SNPs, those where the bases are actually pairs
  ss <- ss[nchar(ss$A1) == 1 & nchar(ss$A2) == 1,]
  ss <- ss[ss$A1 %in% c("A", "C", "G", "T") & ss$A2 %in% c("A", "C", "G", "T"),]
  system(paste0("echo after remove non ACGT SNPs: ", nrow(ss), " >> ", author, "/clean.log"))

  ss <- ss[!((ss$A1 == "A" & ss$A2 == "T") |
             (ss$A1 == "T" & ss$A2 == "A") |
             (ss$A1 == "G" & ss$A2 == "C") |
             (ss$A1 == "C" & ss$A2 == "G")),]

  system(paste0("echo after remove ambigous SNPs: ", nrow(ss), " >> ", author, "/clean.log"))

  #Make sure the imputation reference and the ss objects have the same number of SNP IDs
  ss$RSID <- tolower(ss$RSID)
  ss <- ss[ss$RSID %in% impute$SNP,]
  sub_impute <- impute[impute$SNP %in% ss$RSID,]

  system(paste0("echo after remove SNPs not in the UKBB imputation: ", nrow(ss), " >> ", author, "/clean.log"))

  if(nrow(ss) != nrow(sub_impute)){
    print("SIZE MISMATCH")
  }

  #Extract the chromosomes from the imputation reference object
  #Sometimes the chromosome is not listed for an object
  #But since the rows seem to be ordered we can impute the chromosome based on the surrounding chromosomes
  chr <- sapply(strsplit(sub_impute[,1], ":"), function(x) x[[1]])
  for(i in as.character(1:22)){
    if(i %in% chr){
      min_ind <- min(which(chr==i))
      max_ind <- max(which(chr==i))
      chr[min_ind:max_ind] <- i
    }
  }

  #Finish cleaning up the chr, because "X" could be there instead of 23
  chr[!(chr %in% as.character(1:22))] <- "23"
  chr <- as.numeric(chr)

  #Sort the imputation reference and the ss objects so they are the same order
  #Then assign the imputation chr and pos to the ss object
  ss <- ss[order(ss$RSID)[rank(sub_impute$SNP)],]
  ss <- ss[order(chr, sub_impute$POS),]
  sub_impute <- sub_impute[order(chr, sub_impute$POS),]
  chr <- chr[order(chr)]

  #make the ss object chromosome and position the same as the imputed object
  ss$CHR <- chr
  ss$BP <- sub_impute$POS


  #If there is a name in the parameters that states CHANGE_BOTH assuming that it's not BETA but OR
  #This would just be a simple exponentiation
  if(nrow(params) > 8){
    if(params[9,1] == "CHANGE_BOTH"){
      #will assume that the effect name is for odds ratio
      #following the walds ratio tests can just switch things over assuming a normal distribution
      #This may not be exactly what was completed in the published GWAS, but it should be a good approximation
      #And few methods seem to use the SE column anyway
      ss$BETA <- log(ss$BETA)
      ss$SE <- abs(ss$BETA/qnorm(ss$P))
    }
  }


  #use the bignspr function to flip and reverse the summmary statistics
  colnames(ss) <- c("chr", "pos", "rsid", "a0", "a1", "beta_set", "beta", "p")
  save_impute_cols <- colnames(sub_impute)
  colnames(sub_impute) <- c("loc", "rsid", "pos", "a0", "a1", "maf", "ax", "info")
  sub_impute$chr <- ss$chr
  ss <- snp_match(ss, sub_impute)
  ss <- ss[,c(1,2,5,3,4,6,7,8)]
  colnames(ss) <- c("CHR", "BP", "RSID", "A1", "A2", "SE", "BETA", "P")
  colnames(sub_impute) <- save_impute_cols
  
  system(paste0("echo after removing SNPs that do not flip: ", nrow(ss), " >> ", author, "/clean.log"))


  #Finally add in effective sample size, as defined by PLINK
  notes <- read.table(paste0(author, "/notes"), sep = '\t', stringsAsFactors=F)
  ss$ESS <- round(4/(1/as.numeric(notes[4,1]) + 1/as.numeric(notes[5,1])))

  #According to LDPRED-2 it is a good idea to compare to a validation set following:
  #sd_ss <- 2/(ss$SE * sqrt(ss$ESS))

  #For the validation set I will use summary statistics from the FinnGen Biobank
  #First I have to go through several steps to line everything up
  #finn_notes <- read.table("../finngen_ss/sample_size", stringsAsFactors=F)
  #if(author %in% finn_notes[,4]){
  #  finn_ess <- 4/((1/finn_notes[finn_notes[,4] == author,2]) + (1/finn_notes[finn_notes[,4] == author,3]))
  #  finn_pheno <- finn_notes[finn_notes[,4] == author,1]

  #  finn_ss <- read.table(paste0("../finngen_ss/summary_stats_finngen_r3_", finn_pheno, ".gz"), header =F, stringsAsFactors=F, sep = '\t')  
  #  finn_ss <- finn_ss[finn_ss[,5] %in% ss$RSID,]
  #  dup_rs <- finn_ss[duplicated(finn_ss[,5]), 5]
  #  finn_ss <- finn_ss[!finn_ss[,5] %in% dup_rs,]

  #  finn_rs <- finn_ss[,5]
  #  finn_sd <- finn_ss[,9]
  
  #  finn_sd <- c(finn_sd, rep(0, (nrow(ss) - length(finn_rs))))
  #  finn_rs <- c(finn_rs, ss$RSID[!(ss$RSID %in% finn_rs)])
  #  finn_sd <- finn_sd[order(finn_rs)[rank(ss$RSID)]]

  #  finn_sd <- 2/(finn_sd * sqrt(finn_ess))

  #  good_rs <- ss$RSID[sd_ss < 2 * finn_sd & sd_ss > 0.5 * finn_sd]
  #  good_rs <- c(good_rs, ss$RSID[finn_sd == Inf])
  #  ss <- ss[ss$RSID %in% good_rs,]
  #}
  #system(paste0("echo after removing SNPs that do not align with FinnGen: ", nrow(ss), " >> ", author, "/clean.log"))


  #Now write the answer
  write.table(ss, paste0(author, "/clean_", tolower(author), ".txt"), row.names = F, col.names = T, sep = '\t', quote = F)
  system(paste0("gzip ", author, "/clean_", tolower(author), ".txt"))
}


#all_authors <- read.table("common_files/list_authors", stringsAsFactors=F)
#for(author in all_authors[,1]){
#author <- "european"
for(author in c("african", "eastasian", "hispanic", "total")){
  print(author)
  adjust_ss(author)
}


