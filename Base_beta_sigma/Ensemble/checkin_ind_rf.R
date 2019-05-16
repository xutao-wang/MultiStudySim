lst_0 = readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/lst.rds')
test_0 = readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/test_0.rds')

lst_3 = readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/lst_3.rds')
test_3 = readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/test_3.rds')

lst = get_input_dat(lst_sum)[[1]]
dat_train = get_input_dat(lst_sum)[[2]]
dat_test_ori = readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/test_0.rds')
dat_test = lapply(dat_test_ori,'[[',1)

ensemble_learning <- function(lst,dat_train,dat_test,nk_train){
  
  # Elnet
  con <- caret::trainControl(method = 'cv',number=5,allowParallel = T)
  fit_elnet <- lapply(lst,function(lst) caret::train(y ~ ., data=lst, method = "glmnet",trControl=con,preProc = c("center", "scale")))
  #fitM_elnet <- caret::train(y ~ ., data=dat_train, method = "glmnet",trControl=con,preProc = c("center", "scale"))
  
  #pred_train_M_elnet <- predict(fitM_elnet,as.matrix(dat_train[,-which(names(dat_train) %in% c('y'))]))
  
  #predM_elnet <- lapply(dat_test,function(dat_test) predict(fitM_elnet,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  pred_train_elnet <- lapply(fit_elnet, function(fit_elnet) predict(fit_elnet,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_train_ref <- pred_train_elnet
  
  pred_test_elnet <- lapply(1:length(dat_test),function(x) predict(fit_elnet,dat_test[[x]]  ))
  
  # Neural Net
  fit_nnet <- lapply(lst, function(lst) monmlp::monmlp.fit(x = as.matrix(lst[,names(lst)!='y']), y = as.matrix(lst$y), hidden1 = 2, hidden2 = 2,
                                                           bag = TRUE, iter.max = 500, iter.stopped = 10))
  #fitM_nnet <- monmlp::monmlp.fit(x = as.matrix(dat_train[,names(dat_train)!='y']), y = as.matrix(dat_train$y), hidden1 = 2, hidden2 = 2,
  #                           bag = TRUE, iter.max = 500, iter.stopped = 10)
  #pred_train_M_nnet <- monmlp.predict(as.matrix(dat_train[,-which(names(dat_train) %in% c('y'))]),fitM_nnet)
  
  #predM_nnet <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fitM_nnet))
  pred_train_nnet <- lapply(fit_nnet, function(x) monmlp.predict(as.matrix(dat_train[,-which(names(dat_train) %in% c('y'))]),x))
  pred_test_nnet <- lapply(1:length(dat_test),function(x)  
    lapply(1:length(fit_nnet), function(y)  monmlp.predict(as.matrix(dat_test[[x]][,names(dat_test[[x]])!='y']),fit_nnet[[y]])) )
  
  # Boosting
  fit_boost <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "gbm",trControl=con,preProc = c("center", "scale")))
  
  #fitM_boost <- caret::train(y~.,data=dat_train, method = "gbm",trControl=con,preProc = c("center", "scale"))
  
  #pred_train_M_boost <- predict(fitM_boost,as.matrix(dat_train[,-which(names(dat_train) %in% c('y'))]))
  
  #predM_boost <- lapply(dat_test,function(dat_test) predict(fitM_boost,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  pred_train_boost <- lapply(fit_boost, function(fit_boost) predict(fit_boost,dat_train))
  pred_test_boost <- lapply(1:length(dat_test),function(x) predict(fit_boost,dat_test[[x]]  ))
  
  # RF
  fit_rf <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "ranger",trControl=con,preProc = c("center", "scale")))
  #fitM_rf <- caret::train(y~., data=dat_train,method = "ranger",trControl=con,preProc = c("center", "scale"))
  
  #pred_train_M_rf <- predict(fitM_rf,as.matrix(dat_train[,-which(names(dat_train) %in% c('y'))]))
  
  #predM_rf <- lapply(dat_test,function(dat_test) predict(fitM_rf,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  pred_train_rf <- lapply(fit_rf, function(fit_rf) predict(fit_rf,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test_rf <- lapply(1:length(dat_test),function(x) predict(fit_rf,dat_test[[x]][,-which(names(dat_test[[x]]) %in% c('y'))]  ))
  
  # Merge Multi-study Together
  pred_train <- c(pred_train_elnet,pred_train_nnet,pred_train_boost,pred_train_rf)
  
  pred_test_elnet_mod <- lapply(1:length(pred_test_elnet), function(x) do.call(cbind,pred_test_elnet[[x]]))
  pred_test_nnet_mod <- lapply(1:length(pred_test_nnet), function(x) do.call(cbind,pred_test_nnet[[x]]))
  pred_test_boost_mod <- lapply(1:length(pred_test_boost), function(x) do.call(cbind,pred_test_boost[[x]]))
  pred_test_rf_mod <- lapply(1:length(pred_test_rf), function(x) do.call(cbind,pred_test_rf[[x]]))
  pred_test <- mapply(cbind, pred_test_elnet_mod, pred_test_nnet_mod, pred_test_boost_mod, pred_test_rf_mod,SIMPLIFY=FALSE)
  
  
  nl <- stack_weight_combine(pred_train,dat_train)
  stacking_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(pred_test[[x]] %*%nl,dat_test[[x]]$y))
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  
  # lapply(1:length(pred_test), function(y) apply(pred_test[[y]], 1,mean))
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) 
    Metrics::rmse(lapply(1:length(pred_test), function(y) apply(pred_test[[y]], 1,mean))[[x]],dat_test[[x]]$y))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  
  nl_int <- stack_weight_int(pred_train,dat_train)
  stacking_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(pred_test[[x]] %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]),dat_test[[x]]$y))
  stacking_int_rmse_mean <- mean(unlist(stacking_int_rmse_list))
  
  nl_wz <- stack_weight_wz_ensemble(pred_train_elnet,pred_train_nnet,pred_train_boost,pred_train_rf,dat_train,nk_train,pred_train_ref)
  stacking_wz_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(pred_test[[x]] %*%nl_wz,dat_test[[x]]$y))
  stacking_wz_rmse_mean <- mean(unlist(stacking_wz_rmse_list))
  
  nl_elnet_int <- elnet_weight_int(pred_train, dat_train)
  stacking_elnet_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(pred_test[[x]] %*%nl_elnet_int[-1]+nl_elnet_int[1],dat_test[[x]]$y))
  stacking_elnet_int_rmse_mean <- mean(unlist(stacking_elnet_int_rmse_list))
  
  
  return(list(stacking=stacking_rmse_mean,stakcing_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
              elnet_int = stacking_elnet_int_rmse_mean,simp_avg=simp_avg_rmse_mean
  ))
  
  
  
}
### Beta Sigma = 0 ####
library(grid)

pred_stacking <- lapply(1:length(dat_test), function(x) as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)


pred_stacking_wz <- lapply(1:length(dat_test), function(x) 
  as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)


pred_stacking_int <- lapply(1:length(dat_test), function(x) 
  as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))

pred_elnet_int <- lapply(1:length(dat_test), function(x) 
  as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_elnet_int[-1]+nl_elnet_int[1])

predM <- lapply(dat_test,function(dat_test) predict(fitM,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))

pred_simp_avg <- lapply(1:length(dat_test), function(i) 
  sapply(1:nrow(dat_test[[1]]), function(y) mean(sapply(pred_test[[i]], function(x) x[y]))))

### Elastic Net
fit_single_1 <-  caret::train(y ~ ., data=lst[[1]], method = "glmnet",trControl=con,preProc = c("center", "scale"))
pred_single_1 <- lapply(dat_test,function(dat_test) predict(fit_single_1,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))

### Neural Net
fit_single_1 <- monmlp::monmlp.fit(x = as.matrix(lst[[1]][,names(lst[[1]])!='y']), y = as.matrix(lst[[1]]$y), hidden1 = 2, hidden2 = 2,
                                   bag = TRUE, iter.max = 500, iter.stopped = 10)
pred_single_1 <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fit_single_1))


## Boosting
fit_single_1 <- caret::train(y~., data=lst[[1]],method = "gbm",trControl=con,preProc = c("center", "scale"))
pred_single_1 <- lapply(dat_test,function(dat_test) predict(fit_single_1,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))


### Random Forest
fit_single_1 <- caret::train(y~., data=lst[[1]],method = "ranger",trControl=con,preProc = c("center", "scale"))
pred_single_1 <- lapply(dat_test,function(dat_test) predict(fit_single_1,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))


dat_test_sig0_rf <- rbind(data.frame(pred = unlist(pred_stacking),X_value = unlist(pred_single_1), Weight = rep('Stacking',500)),
                             data.frame(pred = unlist(pred_stacking_wz),X_value = unlist(pred_single_1), Weight = rep('Stacking_wz',500)),
                             data.frame(pred = unlist(pred_stacking_int),X_value = unlist(pred_single_1), Weight = rep('Stacking_int',500)),
                             data.frame(pred = unlist(pred_simp_avg),X_value = unlist(pred_single_1), Weight = rep('Simp_Avg',500)),
                             data.frame(pred = unlist(predM),X_value = unlist(pred_single_1), Weight = rep('Merge',500)),
                             data.frame(pred = do.call(rbind,dat_test)$y,X_value = unlist(pred_single_1), Weight = rep('Actual',500))
)
change_data_structure <- function(dat){
  dat$Actual <- rep(dat[which(dat$Weight %in% 'Actual'),]$pred,nrow(dat)/nrow(dat[which(dat$Weight %in% 'Actual'),]))
  dat[which(dat$Weight %in% 'Actual'),]$pred <- dat$X_value[1:nrow(dat[which(dat$Weight %in% 'Actual'),])]
  
  levels(dat$Weight) <- sub("^Actual$", "Single_Study", levels(dat$Weight))
  dat <- dat[,-which(names(dat) %in% 'X_value')]
  colnames(dat) <- c('pred','Weight','X_value')
  return(dat)
}

#write.csv(dat_test_sig2_1,'~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/dat_test_sig2_1.csv')
#dat_test_sig3_1 <- read.csv('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/dat_test_sig3_1.csv',header=T)
compare_plot <- function(dat,size,learner){
  dat <- change_data_structure(dat)
  dat <- dat[-which(dat$Weight %in% c('Stacking_wz','Stacking_int')),]
  p_plot <- ggplot(dat, aes(x=pred, y=X_value, shape=Weight, color=Weight)) +
    geom_point(alpha=0.5) + 
    guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) + 
    #geom_smooth(se=FALSE, method="loess") +
    #geom_line(data = dat, aes(x = pred, y = X_value)) +
    ggtitle(paste('Beta Sigma =',size,',',learner)) + geom_abline(intercept = 0,slope = 1) +
    xlab('Prediction') + ylab('Actual Value')
  return(p_plot)
  
}

compare_plot_smooth <- function(dat,size,learner){
  dat <- change_data_structure(dat)
  dat <- dat[-which(dat$Weight %in% c('Stacking_wz','Stacking_int')),]
  p_plot <- ggplot(dat, aes(x=pred, y=X_value, shape=Weight, color=Weight)) +
    #geom_point(alpha=0.5) + 
    guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) + 
    geom_smooth(se=FALSE, method="loess") +
    #geom_line(data = dat, aes(x = pred, y = X_value)) +
    ggtitle(paste('Beta Sigma =',size,',',learner)) + geom_abline(intercept = 0,slope = 1) +
    xlab('Prediction') + ylab('Actual Value')
  return(p_plot)
  
}


compare_plot(dat_test_sig0_rf,0,'RF')



grid.arrange(compare_plot(dat_test_sig0_elnet,0,'Elnet'), compare_plot(dat_test_sig0_Nnet,0,'Nnet'), 
             compare_plot(dat_test_sig0_boost,0,'Boosting'), compare_plot(dat_test_sig0_rf,0,'RF'),
             ncol=2)

grid.arrange(compare_plot_smooth(dat_test_sig0_elnet,0,'Elnet'), compare_plot_smooth(dat_test_sig0_Nnet,0,'Nnet'), 
             compare_plot_smooth(dat_test_sig0_boost,0,'Boosting'), compare_plot_smooth(dat_test_sig0_rf,0,'RF'),
             ncol=2)

grid.arrange(compare_plot(dat_test_sig0.75_elnet,0.75,'Elnet'), compare_plot(dat_test_sig0.75_Nnet,0.75,'Nnet'), 
             compare_plot(dat_test_sig0.75_boost,0.75,'Boosting'), compare_plot(dat_test_sig0.75_rf,0.75,'RF'),
             ncol=2)

grid.arrange(compare_plot(dat_test_sig2_elnet,2,'Elnet'), compare_plot(dat_test_sig2_Nnet,2,'Nnet'), 
             compare_plot(dat_test_sig2_boost,2,'Boosting'), compare_plot(dat_test_sig2_rf,2,'RF'),
             ncol=2)

grid.arrange(compare_plot(dat_test_sig3_elnet,3,'Elnet'), compare_plot(dat_test_sig3_Nnet,3,'Nnet'), 
             compare_plot(dat_test_sig3_boost,3,'Boosting'), compare_plot(dat_test_sig3_rf,3,'RF'),
             ncol=2)

grid.arrange(compare_plot_smooth(dat_test_sig3_elnet,3,'Elnet'), compare_plot_smooth(dat_test_sig3_Nnet,3,'Nnet'), 
             compare_plot_smooth(dat_test_sig3_boost,3,'Boosting'), compare_plot_smooth(dat_test_sig3_rf,3,'RF'),
             ncol=2)

save.image('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/Single_Study_pred.RData')

load('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/Single_Study_pred.RData')
