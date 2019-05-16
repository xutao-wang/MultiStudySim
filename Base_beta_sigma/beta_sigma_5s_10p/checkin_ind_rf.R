lst_sum = readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/lst.rds')
lst = get_input_dat(lst_sum)[[1]]
dat_train = get_input_dat(lst_sum)[[2]]
dat_test_ori = readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/test_0.rds')
dat_test = lapply(dat_test_ori,'[[',1)


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
    geom_point(alpha=0.05) + 
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
