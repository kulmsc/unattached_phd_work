library(vroom)
options(warn=2)

pheno <- as.data.frame(vroom("~/athena/ukbiobank/phenotypes/ukb26867.csv.gz", delim = ","))

pheno <- pheno[,c(1,200,2081:2120)]

#job_vals <- unique(c(as.matrix(pheno[,-1])))
job_vals <- table(c(as.matrix(pheno[,-1])))
job_vals <- names(job_vals)[job_vals > 500]
job_vals <- job_vals[!is.na(job_vals)]
job_vals <- job_vals[job_vals != 0]

poss_jobs <- matrix(0, nrow = nrow(pheno), ncol = length(job_vals))

for(i in 1:nrow(pheno)){
  print(i)
  poss_jobs[i,job_vals %in% pheno[i,-1]] <- 1
}
poss_jobs <- as.data.frame(poss_jobs)
colnames(poss_jobs) <- job_vals
poss_jobs$eid <- pheno$eid

saveRDS(poss_jobs, "common_jobs.RDS")
saveRDS(pheno, "all_jobs.RDS")
