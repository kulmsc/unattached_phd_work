
get_missing_p <- function(x){
  return(chisq.test(matrix(x, nrow = 2))$p.value)
}

####################################################

female_miss <- read.table("female.vmiss", stringsAsFactors=F, header=T, comment.char="%")
male_miss <- read.table("male.vmiss", stringsAsFactors=F, header=T, comment.char="%")

test_df <- data.frame(mm = male_miss$MISSING_CT, mo = male_miss$OBS_CT, fm = female_miss$MISSING_CT, fo = female_miss$OBS_CT)
test_df[,2] <- test_df[,2] - test_df[,1]
test_df[,4] <- test_df[,4] - test_df[,3]

missing_p <- apply(test_df, 1, get_missing_p)
missing_p[is.nan(missing_p)] <- 1

######################################################

female_freq <- read.table("female.afreq", stringsAsFactors=F, header=T, comment.char="%")
male_freq <- read.table("male.afreq", stringsAsFactors=F, header=T, comment.char="%")

test_df <- data.frame(mm = round(male_freq$ALT_FREQS * male_freq$OBS_CT),
                      mo = male_freq$OBS_CT - round(male_freq$ALT_FREQS * male_freq$OBS_CT),
                      fm = round(female_freq$ALT_FREQS * female_freq$OBS_CT),
                      fo = female_freq$OBS_CT - round(female_freq$ALT_FREQS * female_freq$OBS_CT))

freq_p <- apply(test_df, 1, get_missing_p)
freq_p[is.nan(freq_p)] <- 0

bad_rsid <- male_freq[freq_p < 5e-8 | missing_p < 1e-50, 2]
write.table(bad_rsid, "bad_rsid", row.names = F, col.names = F, quote = F)
