library(pROC)
library(glmnet)
library(doParallel)
registerDoParallel(4)


predict_pheno_all <- function(train_group, test_group, spec_ind, group_ind, input_df){
  train_xdata <- big_df[train_df$eid %in% train_group,]
  test_xdata <- big_df[train_df$eid %in% test_group,]

  if(file.exists(paste0("minlambs/lam.", ncol(input_df), ".phenoall.RDS"))){
    minlam <- readRDS(paste0("minlambs/lam.", ncol(input_df), ".phenoall.RDS"))
  } else {
    cvmod <- cv.glmnet(y = input_df$pheno[train_df$eid %in% train_group], x = as.matrix(train_xdata), alpha = 1, parallel = TRUE)
    saveRDS(cvmod$lambda.min, paste0("minlambs/lam.", ncol(input_df), ".phenoall.RDS"))
    minlam <- cvmod$lambda.min
  }
  mod <- glmnet(y = input_df$pheno[train_df$eid %in% train_group], x = as.matrix(train_xdata), lambda = minlam, alpha = 1)

  train_preds <- predict(mod, as.matrix(train_xdata))
  test_preds <- predict(mod, as.matrix(test_xdata))


  if(group_ind == 3){
  remove_eids <- list(c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.99)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.99)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.995)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.995)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.9995)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.9995)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.9999)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.9999)]))

  } else if(group_ind == 4){
  remove_eids <- list(c(input_df$eid[input_df$eid %in% train_group][train_preds < quantile(train_preds, 0.01)],
                        input_df$eid[input_df$eid %in% test_group][test_preds < quantile(test_preds, 0.01)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds < quantile(train_preds, 0.005)],
                        input_df$eid[input_df$eid %in% test_group][test_preds < quantile(test_preds, 0.005)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds < quantile(train_preds, 0.0005)],
                        input_df$eid[input_df$eid %in% test_group][test_preds < quantile(test_preds, 0.0005)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds < quantile(train_preds, 0.0001)],
                        input_df$eid[input_df$eid %in% test_group][test_preds < quantile(test_preds, 0.0001)]))
  }

  return(remove_eids[[spec_ind]])
}



predict_residual <- function(train_group, test_group, spec_ind, group_ind, input_df){
  train_xdata <- big_df[train_df$eid %in% train_group,]
  test_xdata <- big_df[train_df$eid %in% test_group,]

  if("sex" %in% colnames(train_df)){
    base_mod <- glm(pheno ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, train_df[train_df$eid %in% train_group,], family = "binomial")
  } else {
    base_mod <- glm(pheno ~ age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, train_df[train_df$eid %in% train_group,], family = "binomial")
  }
  resids <- base_mod$residuals

  if(file.exists(paste0("minlambs/lam.", ncol(input_df), ".residual.RDS"))){
    minlam <- readRDS(paste0("minlambs/lam.", ncol(input_df), ".residual.RDS"))
  } else {
    cvmod <- cv.glmnet(y = resids, x = as.matrix(train_xdata), alpha = 1, parallel = TRUE)
    saveRDS(cvmod$lambda.min, paste0("minlambs/lam.", ncol(input_df), ".residual.RDS"))
    minlam <- cvmod$lambda.min
  }

  #lambda - 0.4
  mod <- glmnet(y = resids, x = as.matrix(train_xdata), lambda = minlam, alpha = 1)

  train_preds <- predict(mod, as.matrix(train_xdata))
  test_preds <- predict(mod, as.matrix(test_xdata))


  remove_eids <- list(c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.99)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.99)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.995)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.995)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.9995)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.9995)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.9999)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.9999)]))


  return(remove_eids[[spec_ind]])
}

#NEED TO FIGURE THIS OUT!!!
predict_pheno_extreme <- function(train_group, test_group, spec_ind, group_ind, input_df){
  train_xdata <- big_df[train_df$eid %in% train_group,]
  test_xdata <- big_df[train_df$eid %in% test_group,]
  train_input <- input_df[input_df$eid %in% train_group,]
  test_input <- input_df[input_df$eid %in% test_group,]

  mod_func <- function(case_eid, control_eid){

  tail_eid <- c(case_eid, control_eid)
  cvmod <- cv.glmnet(y = input_df$pheno[train_df$eid %in% train_group & train_df$eid %in% tail_eid], x = as.matrix(train_xdata[train_group %in% tail_eid,]), alpha = 1)
  mod <- glmnet(y = input_df$pheno[train_df$eid %in% train_group & train_df$eid %in% tail_eid], x = as.matrix(train_xdata[train_group %in% tail_eid,]), lambda = cvmod$lambda.min, alpha = 1)

  train_preds <- predict(mod, as.matrix(train_xdata))
  test_preds <- predict(mod, as.matrix(test_xdata))


  if(group_ind == 5){ #case
  remove_eids <- list(c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.99)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.99)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.995)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.995)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds > quantile(train_preds, 0.9995)],
                        input_df$eid[input_df$eid %in% test_group][test_preds > quantile(test_preds, 0.9995)]))

  } else if(group_ind == 6){ #control
  remove_eids <- list(c(input_df$eid[input_df$eid %in% train_group][train_preds < quantile(train_preds, 0.01)],
                        input_df$eid[input_df$eid %in% test_group][test_preds < quantile(test_preds, 0.01)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds < quantile(train_preds, 0.005)],
                        input_df$eid[input_df$eid %in% test_group][test_preds < quantile(test_preds, 0.005)]),
                      c(input_df$eid[input_df$eid %in% train_group][train_preds < quantile(train_preds, 0.0005)],
                        input_df$eid[input_df$eid %in% test_group][test_preds < quantile(test_preds, 0.0005)]))
  }

  return(remove_eids)
  }


  top_case <- train_input$eid[train_input$pheno == 1][train_input$score[train_input$pheno == 1] > sort(train_input$score[train_input$pheno == 1], decreasing = T)[100]]
  top_control <- train_input$eid[train_input$pheno == 0][train_input$score[train_input$pheno == 0] > sort(train_input$score[train_input$pheno == 0], decreasing = T)[100]]
  bottom_case <- train_input$eid[train_input$pheno == 1][train_input$score[train_input$pheno == 1] < sort(train_input$score[train_input$pheno == 1], decreasing = F)[100]]
  bottom_control <- train_input$eid[train_input$pheno == 0][train_input$score[train_input$pheno == 0] < sort(train_input$score[train_input$pheno == 0], decreasing = F)[100]]

  #top_case <- train_xdata$eid[train_xdata$pheno == 1][train_xdata$score[train_xdata$pheno == 1] > sort(train_xdata$score[train_xdata$pheno == 1], decreasing = T)[100]]
  #top_control <- train_xdata$eid[train_xdata$pheno == 0][train_xdata$score[train_xdata$pheno == 0] > sort(train_xdata$score[train_xdata$pheno == 0], decreasing = T)[100]]
  #bottom_case <- train_xdata$eid[train_xdata$pheno == 1][train_xdata$score[train_xdata$pheno == 1] < sort(train_xdata$score[train_xdata$pheno == 1], decreasing = F)[100]]
  #bottom_control <- train_xdata$eid[train_xdata$pheno == 0][train_xdata$score[train_xdata$pheno == 0] < sort(train_xdata$score[train_xdata$pheno == 0], decreasing = F)[100]]
  
  remove_eids <- c(mod_func(top_case, top_control),
                 mod_func(bottom_case, bottom_control),
                 mod_func(top_control, bottom_case))
  
  return(remove_eids[[spec_ind]])
}


