

survey_data <- read.table("imptoved_features.txt", stringsAsFactors=F, header=T)
census <- read.table("/home/kulmsc/athena/lior_project/polysocial_score_predictor/census.txt", stringsAsFactors=F, header=T)
job_data <- readRDS("common_jobs.RDS")
colnames(survey_data) <- paste0("SURVEY_", colnames(survey_data))
colnames(census) <- paste0("CENSUS_", colnames(census))
colnames(job_data) <- paste0("JOB_", colnames(job_data))

cancer_codes <- read.table("combo_output/diag.cancer.txt.gz", stringsAsFactors=F, header=F)
noncancer_codes <- read.table("combo_output/diag.noncancer.txt.gz", stringsAsFactors=F, header=F)
icd10_codes <- read.table("combo_output/inci.icd10.txt.gz", stringsAsFactors=F, header=F)
icd9_codes <- read.table("combo_output/inci.icd9.txt.gz", stringsAsFactors=F, header=F)
meds_codes <- read.table("combo_output/diag.meds.txt.gz", stringsAsFactors=F, header=F)
opcs_codes <- read.table("combo_output/inci.opcs.txt.gz", stringsAsFactors=F, header=F)

cancer_names <- read.table("raw_output/diag.cancer.bentham.txt.gz", stringsAsFactors=F)
noncancer_names <- read.table("raw_output/diag.noncancer.bentham.txt.gz", stringsAsFactors=F)
icd10_names <- read.table("raw_output/diag.icd10.bentham.txt.gz", stringsAsFactors=F)
icd9_names <- read.table("raw_output/diag.icd9.bentham.txt.gz", stringsAsFactors=F)
meds_names <- read.table("raw_output/diag.meds.bentham.txt.gz", stringsAsFactors=F)
opcs_names <- read.table("raw_output/diag.opcs.bentham.txt.gz", stringsAsFactors=F)

colnames(cancer_codes) <- paste0("CANCER_", cancer_names[,1])
colnames(noncancer_codes) <- paste0("NONCANCER_", noncancer_names[,1])
colnames(icd10_codes) <- paste0("ICD10_", icd10_names[,1])
colnames(icd9_codes) <- paste0("ICD9_", icd9_names[,1])
colnames(meds_codes) <- paste0("MEDS_", meds_names[,1])
colnames(opcs_codes) <- paste0("OPCS_", opcs_names[,1])

codes_eid <- read.table("big_eid.txt", stringsAsFactors=F)
codes_eid <- codes_eid[1:nrow(icd10_codes),,drop=F]

cancer_codes$eid <- codes_eid[,1]
noncancer_codes$eid <- codes_eid[,1]
icd10_codes$eid <- codes_eid[,1]
icd9_codes$eid <- codes_eid[,1]
meds_codes$eid <- codes_eid[,1]
opcs_codes$eid <- codes_eid[,1]

colnames(census)[colnames(census) == "CENSUS_eid"] <- "eid"
colnames(survey_data)[colnames(survey_data) == "SURVEY_eid"] <- "eid"
colnames(job_data)[colnames(job_data) == "JOB_eid"] <- "eid"

big_data <- merge(survey_data, census, by = "eid")
big_data <- merge(big_data, job_data, by = "eid")
big_data <- merge(big_data, cancer_codes, by = "eid")
big_data <- merge(big_data, noncancer_codes, by = "eid")
big_data <- merge(big_data, icd10_codes, by = "eid")
big_data <- merge(big_data, icd9_codes, by = "eid")
big_data <- merge(big_data, meds_codes, by = "eid")
big_data <- merge(big_data, opcs_codes, by = "eid")

saveRDS(big_data, "big_data.RDS")
