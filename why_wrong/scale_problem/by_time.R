
#batch individuals by their age and time in the study
#note that the extreme values are only good for plotting

#have already confirmed the general nature of the problem graphically
#so now just go through and guess
run_time <- function(tort){
time_group <- list()
if(tort == "train"){
  timeline <- readRDS("../comp_outcomes/timeline_res/all_eid.time.RDS")
} else {
  timeline <- readRDS("../comp_outcomes/timeline_res/test_eid.time.RDS")
}
timeline$last_icd_diag <- as.Date(timeline$last_icd_diag)
timeline$last_diag <- as.Date(timeline$last_diag)

#note that you can have NA time from last diag, and positive diagnoses
time_from_diag <- as.Date("2021-01-01") - timeline$last_diag
time_group[[1]] <- timeline$eid[time_from_diag > quantile(time_from_diag, 0.95, na.rm=T) & !is.na(time_from_diag)]
time_group[[2]] <- timeline$eid[time_from_diag > quantile(time_from_diag, 0.99, na.rm=T) & !is.na(time_from_diag)]
time_group[[3]] <- timeline$eid[time_from_diag > quantile(time_from_diag, 0.999, na.rm=T) & !is.na(time_from_diag)]
time_group[[4]] <- timeline$eid[is.na(time_from_diag)]
time_group[[5]] <- c(time_group[[4]], time_group[[2]])

time_group[[6]] <- timeline$eid[timeline$total_visit <= 1]
time_group[[7]] <- timeline$eid[timeline$total_visit <= 0]
time_group[[8]] <- timeline$eid[timeline$total_visit > quantile(timeline$total_visit, 0.95)]
time_group[[9]] <- timeline$eid[timeline$total_visit > quantile(timeline$total_visit, 0.99)]
time_group[[10]] <- timeline$eid[timeline$total_visit > quantile(timeline$total_visit, 0.999)]
time_group[[11]] <- c(time_group[[6]], time_group[[8]])
time_group[[12]] <- c(time_group[[7]], time_group[[10]])

#all death
timeline$death_date <- as.Date(timeline$death_date)
time_from_death <- as.Date("2021-01-01") - timeline$death_date
time_group[[13]] <- timeline$eid[timeline$is_dead == 1]
time_group[[14]] <- timeline$eid[time_from_death > 2000 & timeline$is_dead == 1]
time_group[[15]] <- timeline$eid[time_from_death > 4000 & timeline$is_dead == 1]

#young
time_group[[16]] <- timeline$eid[timeline$age < quantile(timeline$age, 0.05)]
time_group[[17]] <- timeline$eid[timeline$age < quantile(timeline$age, 0.01)]
time_group[[18]] <- timeline$eid[timeline$age < quantile(timeline$age, 0.001)]

if(tort == "train"){
saveRDS(time_group, "out_batch/timeline_groups.train.RDS")
} else {
saveRDS(time_group, "out_batch/timeline_groups.test.RDS")
}


}


run_time("train")
run_time("test")
