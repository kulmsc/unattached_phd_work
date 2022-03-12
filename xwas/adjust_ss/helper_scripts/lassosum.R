library(lassosum)

print("starting r script")
#args are: d, chr, low_author
args <- commandArgs(trailingOnly=TRUE)

print("reading ss")
ss <- read.table(paste0("temp_files/ss.", args[3], ".", args[2]), stringsAsFactors=F, header=T)
ss <- ss[!(ss$BP %in% ss$BP[duplicated(ss$BP)]),]
ref.bfile <- paste0("geno_files/", args[3], ".", args[2])
print("done reading ss")

###########

pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]

pheno <- read.table(paste0("../get_pheno/pheno_defs/diag.", args[3], ".txt.gz"), stringsAsFactors=F, header=F)
output <- data.frame(eid = pheno_eids, eid2 = pheno_eids, pheno = apply(pheno, 1, function(x) 1 %in% x[1:4]) * 1 + 1)

train <- read.table("~/athena/doc_score/qc/cv_files/train_eid.0.5.txt", stringsAsFactors=F)
output <- output[output[,1] %in% train[,1],]
ess <- 4/((1/sum(output$pheno == 2)) + (1/sum(output$pheno == 1)))

###########

ss <- ss[!is.na(ss$P) & !is.na(ss$BETA),]
cor <- p2cor(p = ss$P, n = ess, sign=ss$BETA)

out <- lassosum.pipeline(cor=cor, chr=ss$CHR, pos=ss$BP, 
                         A1=ss$A1, A2=ss$A2, # A2 is not required but advised
                         ref.bfile=ref.bfile, LDblocks = "EUR.hg19", sample = 5000,
                         s = c(0.1, 0.4, 0.8), lambda = c(0.002, 0.004, 0.007, 0.1))

new_ss <- ss[ss$BP %in% out$sumstats[,2],]

k <- 1
for(i in 1:3){
  for(j in 1:4){
    new_ss$BETA <- out$beta[[i]][,j]
    write.table(new_ss, paste0("~/athena/xwas/mod_sets/", args[3], "/", args[3], ".", args[2],  ".lassosum.", k ,".ss"), row.names = F, col.names = T, sep = '\t', quote = F)
    k <- k + 1
  }
} 
