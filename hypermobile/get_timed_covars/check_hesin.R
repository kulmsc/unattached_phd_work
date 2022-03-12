options(warn=2)

hesin_diag <- read.table("~/athena/ukbiobank/hesin/hesin_diag.txt", stringsAsFactors=F, header=T, sep="\t")
hesin_oper <- read.table("~/athena/ukbiobank/hesin/hesin_oper.txt", stringsAsFactors=F, header=T, sep="\t")
hesin <- read.table("~/athena/ukbiobank/hesin/hesin.txt", stringsAsFactors=F, header=T, sep="\t")

hesin$epistart <- as.Date(hesin$epistart, "%d/%m/%Y")
hesin$epiend <- as.Date(hesin$epiend, "%d/%m/%Y")
hesin$admidate <- as.Date(hesin$admidate, "%d/%m/%Y")
hesin$disdate <- as.Date(hesin$disdate, "%d/%m/%Y")

primary_eid <- read.table("primary_eid", stringsAsFactors=F)

hesin_diag <- hesin_diag[hesin_diag[,1] %in% primary_eid[,1],]
hesin_oper <- hesin_oper[hesin_oper[,1] %in% primary_eid[,1],]
hesin <- hesin[hesin[,1] %in% primary_eid[,1],]
#hesin_diag$short_icd10 <- substr(hesin_diag$diag_icd10, 1, 3)


primary_icd <- read.table("kr_primary_icd", stringsAsFactors=F)
primary_opcs <- strsplit(primary_icd[2,2], "\\|", fixed = T)[[1]]
primary_icd <- strsplit(primary_icd[1,2], "\\|", fixed = T)[[1]]
risk_icd <- read.table("kr_risk_icd", stringsAsFactors=F)
risk_icd_list <- list()
for(i in 1:nrow(risk_icd)){
  risk_icd_list[[i]] <- strsplit(risk_icd[i,1], "\\|", fixed = T)[[1]]
}


get_date <- function(eid, icd_code, hesin_file, hesin_index){ #hesin_diag = 6, hesin_oper = 8
  if("XXX" %in% icd_code){
    return(NA)
  }
  sub_diag <- hesin_file[hesin_file[,1] == eid & grepl(paste(icd_code, collapse = "|"), hesin_file[,hesin_index]),]
  if(nrow(sub_diag) > 0){
    sub_hesin <- hesin[hesin[,1] == eid & hesin$ins_index %in% sub_diag$ins_index,]
    return(min(c(sub_hesin$epistart, sub_hesin$admidate, sub_hesin$disdate), na.rm=T))
  } else {
    return(NA)
  }
}


###############################

covar_counter <- data.frame(matrix(0, nrow = nrow(primary_eid), ncol =  nrow(risk_icd)))
for(i in 1:nrow(primary_eid)){
  
  if(!"XXX" %in% primary_icd & !"XXX" %in% primary_opcs){
    comp_date <- min(c(get_date(primary_eid[i,1], primary_icd, hesin_diag, 7),
                   get_date(primary_eid[i,1], primary_opcs, hesin_oper, 8)), na.rm=T)
  } else if(!"XXX" %in% primary_icd){
    comp_date <- get_date(primary_eid[i,1], primary_icd, hesin_diag, 7)
  } else if(!"XXX" %in% primary_opcs){
    comp_date <- get_date(primary_eid[i,1], primary_opcs, hesin_oper, 8)
  }

  for(j in 1:nrow(risk_icd)){
    try_date <- get_date(primary_eid[i,1], risk_icd_list[[j]], hesin_diag, 7)
    if(!is.na(try_date)){
      print(try_date)
      if(try_date < comp_date){
        print(j)
        covar_counter[i,j] <- 1
      }
    }
  }

}


for(i in 1:nrow(risk_icd)){
  write.table(primary_eid[covar_counter[,i] == 1,], paste0("kr_covar.", risk_icd[i,2], ".eid"), row.names = F, col.names = F, quote = F)
}

