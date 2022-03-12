

author <- commandArgs(trailingOnly=TRUE)[1]

all_auc <- rep(NA, 18)
all_or <- rep(NA, 18)

for(i in 1:18){
  all_res <- readRDS(paste0("acc_results/all_indiv.time.", i, ".", author, ".RDS"))
  better_res <- readRDS(paste0("acc_results/better_indiv.time.", i, ".", author, ".RDS"))

  all_auc[i] <- better_res[["score"]][["auc"]][2] - all_res[["score"]][["auc"]][2]
  all_or[i] <- better_res[["score"]][["or"]][5,2] - all_res[["score"]][["or"]][5,2]
}

i <- which.max(all_or)
print(i)
print(all_auc)

all_res <- readRDS(paste0("acc_results/all_indiv.time.", i, ".", author, ".RDS"))
better_res <- readRDS(paste0("acc_results/better_indiv.time.", i, ".", author, ".RDS"))

print("score diff")
print(better_res[["score"]][["auc"]][2] - all_res[["score"]][["auc"]][2])

print("score perc gain")
print((better_res[["score"]][["auc"]][2] - all_res[["score"]][["auc"]][2])/all_res[["score"]][["auc"]][2])

print("score imp diff")
print((better_res[["score"]][["auc"]][2] - better_res[["base"]][["auc"]][2]) - (all_res[["score"]][["auc"]][2] - all_res[["base"]][["auc"]][2]))

print("")
print("")
print(paste("all - base: ", all_res[["base"]][["auc"]][2]))
print(paste("all - score: ", all_res[["score"]][["auc"]][2]))
print(paste("better - base: ", better_res[["base"]][["auc"]][2]))
print(paste("better - score: ", better_res[["score"]][["auc"]][2]))

print("")
print("")
print(paste("OR better:",better_res[["score"]][["or"]][5,2]))
print(paste("OR all:",all_res[["score"]][["or"]][5,2]))
