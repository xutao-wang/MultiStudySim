library(grid)
get_input_dat <- function(list){
  maindat <- lapply(list,'[[',1)
  dat_train <- do.call(rbind,maindat)
  return(list(maindat_list=maindat,dat_train=dat_train))
}
lst_sum_1 = readRDS('~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/lst.rds')
lst = get_input_dat(lst_sum)[[1]]
dat_train = get_input_dat(lst_sum)[[2]]
test_1 = readRDS('~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/test_1.rds')
#test_1 = lapply(dat_test_ori,'[[',1)

unlist(lapply(lst_sum,'[[',3))

elnet_pred_value <- function(lst_ori,dat_test){
  lst <- get_input_dat(lst_ori)[[1]]
  dat_train <- get_input_dat(lst_ori)[[2]]
  con <- caret::trainControl(method = 'cv',number=5,allowParallel = T)
  fit <- lapply(lst,function(lst) caret::train(y ~ X + X:X1, data=lst, method = "glmnet",trControl=con,preProc = c("center", "scale")))
  #insamp_rmse_list <- lapply(1:length(fit),function(x) 
  #  Metrics::rmse(predict(fit[[x]],lst[[x]][,-which(names(lst[[x]]) %in% c('y'))]),lst[[x]]$y))
  
  fitM <- caret::train(y ~ X + X:X1, data=dat_train, method = "glmnet",trControl=con,preProc = c("center", "scale"))
  
  pred_train <- lapply(fit, function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  
  pred_test <- lapply(1:length(dat_test),function(x) predict(fit,dat_test[[x]]  ))
  
  nl <- stack_weight_combine(pred_train,dat_train)
  pred_stacking <- lapply(1:length(dat_test), function(x) as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)
  
  nl_int <- stack_weight_int(pred_train,dat_train)
  pred_stacking_int <- lapply(1:length(dat_test), function(x) 
    as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))
  
  nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
  pred_stacking_wz <- lapply(1:length(dat_test), function(x) 
    as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  pred_simp_avg <- lapply(1:length(dat_test), function(i) 
    sapply(1:nrow(dat_test[[1]]), function(y) mean(sapply(pred_test[[i]], function(x) x[y]))))
  
  index <- sapply(lst_ori,'[[',3)
  
  # From positive side
  fit_single_1 <-  caret::train(y ~ X + X:X1, data=lst[[as.numeric(which(index==2)[1])]], method = "glmnet",trControl=con,preProc = c("center", "scale"))
  pred_single_1 <- lapply(dat_test,function(dat_test) predict(fit_single_1,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  # From negative side
  fit_single_2 <-  caret::train(y ~ X + X:X1, data=lst[[as.numeric(which(index==1)[1])]], method = "glmnet",trControl=con,preProc = c("center", "scale"))
  pred_single_2 <- lapply(dat_test,function(dat_test) predict(fit_single_2,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  dat <- rbind(data.frame(pred = unlist(pred_stacking),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking',500)),
               data.frame(pred = unlist(pred_stacking_wz),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_wz',500)),
               data.frame(pred = unlist(pred_stacking_int),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_int',500)),
               data.frame(pred = unlist(pred_simp_avg),X_value = do.call(rbind,dat_test)$X, Weight = rep('Simp_Avg',500)),
               data.frame(pred = unlist(predM),X_value = do.call(rbind,dat_test)$X, Weight = rep('Merge',500)),
               data.frame(pred = unlist(pred_single_1),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 1',500)),
               data.frame(pred = unlist(pred_single_2),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 2',500)),
               data.frame(pred = do.call(rbind,dat_test)$y,X_value = do.call(rbind,dat_test)$X, Weight = rep('Actual',500))
  )
  return(dat)
}
elnet_pred_0.1 <- elnet_pred_value(lst_sum_0.1,lapply(test_0.1,'[[',1))
elnet_pred_0.5 <- elnet_pred_value(lst_sum_0.5,lapply(test_0.5,'[[',1))
elnet_pred_1 <- elnet_pred_value(lst_sum_1,lapply(test_1,'[[',1))
elnet_pred_4 <- elnet_pred_value(lst_sum_4,lapply(test_4,'[[',1))

nnet_pred_value <- function(lst_ori,dat_test){
  
  lst <- get_input_dat(lst_ori)[[1]]
  dat_train <- get_input_dat(lst_ori)[[2]]
  
  fit <- lapply(lst, function(lst) monmlp::monmlp.fit(x = as.matrix(lst[,names(lst)!='y']), y = as.matrix(lst$y), hidden1 = 2, hidden2 = 2,
                                                      bag = TRUE, iter.max = 500, iter.stopped = 10))
  
  fitM <- monmlp::monmlp.fit(x = as.matrix(dat_train[,names(dat_train)!='y']), y = as.matrix(dat_train$y), hidden1 = 2, hidden2 = 2,
                             bag = TRUE, iter.max = 500, iter.stopped = 10)
  
  pred_train <- lapply(fit, function(x) monmlp.predict(as.matrix(dat_train[,-which(names(dat_train) %in% c('y'))]),x))
  pred_test <- lapply(1:length(dat_test),function(x)  
    lapply(1:length(fit), function(y)  monmlp.predict(as.matrix(dat_test[[x]][,names(dat_test[[x]])!='y']),fit[[y]])) )
  
  nl <- stack_weight_combine(pred_train,dat_train)
  pred_stacking <- lapply(1:length(dat_test), function(x) as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)
  
  nl_int <- stack_weight_int(pred_train,dat_train)
  pred_stacking_int <- lapply(1:length(dat_test), function(x) 
    as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))
  
  nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
  pred_stacking_wz <- lapply(1:length(dat_test), function(x) 
    as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)
  
  predM <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fitM))
  pred_simp_avg <- lapply(1:length(dat_test), function(i) 
    sapply(1:nrow(dat_test[[1]]), function(y) mean(sapply(pred_test[[i]], function(x) x[y]))))
  
  index <- sapply(lst_ori,'[[',3)
  
  fit_single_1 <- monmlp::monmlp.fit(x = as.matrix(lst[[as.numeric(which(index == 2)[1])]][,names(lst[[1]])!='y']), y = as.matrix(lst[[as.numeric(which(index == 2)[1])]]$y), hidden1 = 2, hidden2 = 2,
                                     bag = TRUE, iter.max = 500, iter.stopped = 10)
  pred_single_1 <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fit_single_1))
  
  fit_single_2 <- monmlp::monmlp.fit(x = as.matrix(lst[[as.numeric(which(index == 1)[1])]][,names(lst[[2]])!='y']), y = as.matrix(lst[[as.numeric(which(index == 1)[1])]]$y), hidden1 = 2, hidden2 = 2,
                                     bag = TRUE, iter.max = 500, iter.stopped = 10)
  pred_single_2 <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fit_single_2))
  
  dat <- rbind(data.frame(pred = unlist(pred_stacking),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking',500)),
               data.frame(pred = unlist(pred_stacking_wz),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_wz',500)),
               data.frame(pred = unlist(pred_stacking_int),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_int',500)),
               data.frame(pred = unlist(pred_simp_avg),X_value = do.call(rbind,dat_test)$X, Weight = rep('Simp_Avg',500)),
               data.frame(pred = unlist(predM),X_value = do.call(rbind,dat_test)$X, Weight = rep('Merge',500)),
               data.frame(pred = unlist(pred_single_1),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 1',500)),
               data.frame(pred = unlist(pred_single_2),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 2',500)),
               data.frame(pred = do.call(rbind,dat_test)$y,X_value = do.call(rbind,dat_test)$X, Weight = rep('Actual',500))
  )
  return(dat)
  
}

nnet_pred_0.1 <- nnet_pred_value(lst_sum_0.1,lapply(test_0.1,'[[',1))
nnet_pred_0.5 <- nnet_pred_value(lst_sum_0.5,lapply(test_0.5,'[[',1))
nnet_pred_1 <- nnet_pred_value(lst_sum_1,lapply(test_1,'[[',1))
nnet_pred_4 <- nnet_pred_value(lst_sum_4,lapply(test_4,'[[',1))


boost_pred_value <- function(lst_ori,dat_test){
  lst <- get_input_dat(lst_ori)[[1]]
  dat_train <- get_input_dat(lst_ori)[[2]]
  con <- caret::trainControl(method = 'cv',number=5,allowParallel = T)
  fit <- lapply(lst,function(lst) caret::train(y ~ ., data=lst, method = "gbm",trControl=con,preProc = c("center", "scale")))
  #insamp_rmse_list <- lapply(1:length(fit),function(x) 
  #  Metrics::rmse(predict(fit[[x]],lst[[x]][,-which(names(lst[[x]]) %in% c('y'))]),lst[[x]]$y))
  
  fitM <- caret::train(y ~ ., data=dat_train, method = "gbm",trControl=con,preProc = c("center", "scale"))
  
  pred_train <- lapply(fit, function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  
  pred_test <- lapply(1:length(dat_test),function(x) predict(fit,dat_test[[x]]  ))
  
  nl <- stack_weight_combine(pred_train,dat_train)
  pred_stacking <- lapply(1:length(dat_test), function(x) as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)
  
  nl_int <- stack_weight_int(pred_train,dat_train)
  pred_stacking_int <- lapply(1:length(dat_test), function(x) 
    as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))
  
  nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
  pred_stacking_wz <- lapply(1:length(dat_test), function(x) 
    as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  pred_simp_avg <- lapply(1:length(dat_test), function(i) 
    sapply(1:nrow(dat_test[[1]]), function(y) mean(sapply(pred_test[[i]], function(x) x[y]))))
  
  index <- sapply(lst_ori,'[[',3)
  
  # From positive side
  fit_single_1 <-  caret::train(y ~ ., data=lst[[as.numeric(which(index==2)[1])]], method = "gbm",trControl=con,preProc = c("center", "scale"))
  pred_single_1 <- lapply(dat_test,function(dat_test) predict(fit_single_1,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  # From negative side
  fit_single_2 <-  caret::train(y ~ ., data=lst[[as.numeric(which(index==1)[1])]], method = "gbm",trControl=con,preProc = c("center", "scale"))
  pred_single_2 <- lapply(dat_test,function(dat_test) predict(fit_single_2,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  dat <- rbind(data.frame(pred = unlist(pred_stacking),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking',500)),
               data.frame(pred = unlist(pred_stacking_wz),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_wz',500)),
               data.frame(pred = unlist(pred_stacking_int),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_int',500)),
               data.frame(pred = unlist(pred_simp_avg),X_value = do.call(rbind,dat_test)$X, Weight = rep('Simp_Avg',500)),
               data.frame(pred = unlist(predM),X_value = do.call(rbind,dat_test)$X, Weight = rep('Merge',500)),
               data.frame(pred = unlist(pred_single_1),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 1',500)),
               data.frame(pred = unlist(pred_single_2),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 2',500)),
               data.frame(pred = do.call(rbind,dat_test)$y,X_value = do.call(rbind,dat_test)$X, Weight = rep('Actual',500))
  )
  return(dat)
}

boost_pred_0.1 <- boost_pred_value(lst_sum_0.1,lapply(test_0.1,'[[',1))
boost_pred_0.5 <- boost_pred_value(lst_sum_0.5,lapply(test_0.5,'[[',1))
boost_pred_1 <- boost_pred_value(lst_sum_1,lapply(test_1,'[[',1))
boost_pred_4 <- boost_pred_value(lst_sum_4,lapply(test_4,'[[',1))


treebag_pred_value <- function(lst_ori,dat_test){
  lst <- get_input_dat(lst_ori)[[1]]
  dat_train <- get_input_dat(lst_ori)[[2]]
  con <- caret::trainControl(method = 'cv',number=5,allowParallel = T)
  fit <- lapply(lst,function(lst) caret::train(y ~ ., data=lst, method = "treebag",trControl=con,preProc = c("center", "scale")))

  fitM <- caret::train(y ~ ., data=dat_train, method = "treebag",trControl=con,preProc = c("center", "scale"))
  
  pred_train <- lapply(fit, function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  
  pred_test <- lapply(1:length(dat_test),function(x) predict(fit,dat_test[[x]]  ))
  
  nl <- stack_weight_combine(pred_train,dat_train)
  pred_stacking <- lapply(1:length(dat_test), function(x) as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)
  
  nl_int <- stack_weight_int(pred_train,dat_train)
  pred_stacking_int <- lapply(1:length(dat_test), function(x) 
    as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))
  
  nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
  pred_stacking_wz <- lapply(1:length(dat_test), function(x) 
    as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  pred_simp_avg <- lapply(1:length(dat_test), function(i) 
    sapply(1:nrow(dat_test[[1]]), function(y) mean(sapply(pred_test[[i]], function(x) x[y]))))
  
  index <- sapply(lst_ori,'[[',3)
  
  # From positive side
  fit_single_1 <-  caret::train(y ~ ., data=lst[[as.numeric(which(index==2)[1])]], method = "treebag",trControl=con,preProc = c("center", "scale"))
  pred_single_1 <- lapply(dat_test,function(dat_test) predict(fit_single_1,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  # From negative side
  fit_single_2 <-  caret::train(y ~ ., data=lst[[as.numeric(which(index==1)[1])]], method = "treebag",trControl=con,preProc = c("center", "scale"))
  pred_single_2 <- lapply(dat_test,function(dat_test) predict(fit_single_2,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  
  dat <- rbind(data.frame(pred = unlist(pred_stacking),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking',500)),
               data.frame(pred = unlist(pred_stacking_wz),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_wz',500)),
               data.frame(pred = unlist(pred_stacking_int),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_int',500)),
               data.frame(pred = unlist(pred_simp_avg),X_value = do.call(rbind,dat_test)$X, Weight = rep('Simp_Avg',500)),
               data.frame(pred = unlist(predM),X_value = do.call(rbind,dat_test)$X, Weight = rep('Merge',500)),
               data.frame(pred = unlist(pred_single_1),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 1',500)),
               data.frame(pred = unlist(pred_single_2),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 2',500)),
               data.frame(pred = do.call(rbind,dat_test)$y,X_value = do.call(rbind,dat_test)$X, Weight = rep('Actual',500))
  )
  return(dat)
}

treebag_pred_0.1 <- treebag_pred_value(lst_sum_0.1,lapply(test_0.1,'[[',1))
treebag_pred_0.5 <- treebag_pred_value(lst_sum_0.5,lapply(test_0.5,'[[',1))
treebag_pred_1 <- treebag_pred_value(lst_sum_1,lapply(test_1,'[[',1))
treebag_pred_4 <- treebag_pred_value(lst_sum_4,lapply(test_4,'[[',1))

Metrics::rmse(treebag_pred_4[which(treebag_pred_1$Weight%in% 'Stacking_int'),]$pred,
              treebag_pred_4[which(treebag_pred_1$Weight%in% 'Actual'),]$pred)

Metrics::rmse(treebag_pred_4[which(treebag_pred_1$Weight%in% 'Merge'),]$pred,
              treebag_pred_4[which(treebag_pred_1$Weight%in% 'Actual'),]$pred)
compare_plot <- function(dat,size,learner){
  dat <- dat[-which(dat$Weight %in% c('Single 1','Single 2')),]
  p_plot <- ggplot(dat, aes(x=X_value, y=pred, shape=Weight, color=Weight)) +
    geom_point(alpha=0.1) + guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) +
    ggtitle(paste(learner,'predictors vs. X for different weighiting methods X Sigma = ', size))
  return(p_plot)
  
}
grid.arrange(compare_plot(elnet_pred_0.1,0.1,'Elnet'),compare_plot(nnet_pred_0.1,0.1,'Nnet'),
             compare_plot(boost_pred_0.1,0.1,'Boosting'),compare_plot(treebag_pred_0.1,0.1,'tree'),
             ncol=2,top=textGrob('Single covariate X vs. Y',gp=gpar(fontsize=20,font=3)))

grid.arrange(compare_plot(elnet_pred_0.5,0.5,'Elnet'),compare_plot(nnet_pred_0.5,0.5,'Nnet'),
             compare_plot(boost_pred_0.5,0.5,'Boosting'),compare_plot(treebag_pred_0.5,0.5,'tree'),
             ncol=2,top=textGrob('Single covariate X vs. Y',gp=gpar(fontsize=20,font=3)))

grid.arrange(compare_plot(elnet_pred_1,1,'Elnet'),compare_plot(nnet_pred_1,1,'Nnet'),
             compare_plot(boost_pred_1,1,'Boosting'),compare_plot(treebag_pred_1,1,'tree'),
             ncol=2,top=textGrob('Single covariate X vs. Y',gp=gpar(fontsize=20,font=3)))

grid.arrange(compare_plot(elnet_pred_4,4,'Elnet'),compare_plot(nnet_pred_4,4,'Nnet'),
             compare_plot(boost_pred_4,4,'Boosting'),compare_plot(treebag_pred_4,4,'tree'),
             ncol=2,top=textGrob('Single covariate X vs. Y',gp=gpar(fontsize=20,font=3)))


get_single_dat <- function(dat){
  dat <- dat[-which(dat$Weight %in% c('Stacking_wz','Stacking_int')),]
  return(dat)
}

compare_single_plot <- function(dat,size,learner){
  dat <- get_single_dat(dat)
  p_plot <- ggplot(dat, aes(x=X_value, y=pred, shape=Weight, color=Weight)) +
    geom_point(alpha=0.1) + guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) +
    ggtitle(paste(learner,'predictors vs. Single Predictor X Sigma = ', size)) +
    xlab('X value') + ylab('Predictors')
  return(p_plot)
  
}
grid.arrange(compare_single_plot(elnet_pred_0.1,0.1,'Elnet'),compare_single_plot(nnet_pred_0.1,0.1,'Nnet'),
             compare_single_plot(boost_pred_0.1,0.1,'Boosting'),compare_single_plot(treebag_pred_0.1,0.1,'tree'),
             ncol=2,top=textGrob('Single covariate X vs. Y',gp=gpar(fontsize=20,font=3)))

grid.arrange(compare_single_plot(elnet_pred_0.5,0.5,'Elnet'),compare_single_plot(nnet_pred_0.5,0.5,'Nnet'),
             compare_single_plot(boost_pred_0.5,0.5,'Boosting'),compare_single_plot(treebag_pred_0.5,0.5,'tree'),
             ncol=2,top=textGrob('Single covariate X vs. Y',gp=gpar(fontsize=20,font=3)))

grid.arrange(compare_single_plot(elnet_pred_1,1,'Elnet'),compare_single_plot(nnet_pred_1,1,'Nnet'),
             compare_single_plot(boost_pred_1,1,'Boosting'),compare_single_plot(treebag_pred_1,1,'tree'),
             ncol=2,top=textGrob('Single covariate X vs. Y',gp=gpar(fontsize=20,font=3)))

grid.arrange(compare_single_plot(elnet_pred_4,4,'Elnet'),compare_single_plot(nnet_pred_4,4,'Nnet'),
             compare_single_plot(boost_pred_4,4,'Boosting'),compare_single_plot(treebag_pred_4,4,'tree'),
             ncol=2,top=textGrob('Single covariate X vs. Y',gp=gpar(fontsize=20,font=3)))


#compare_single_plot(compare_dat_boost,1,'Boosting','Single 1')
grid.arrange(compare_single_plot(compare_dat_elnet,1,'Elnet','Single 1'),compare_single_plot(compare_dat_Nnet,1,'Nnet','Single 1'),
             compare_single_plot(compare_dat_boost,1,'Boosting','Single 1'),compare_single_plot(compare_dat_tree,1,'tree','Single 1'),
             ncol=2)
grid.arrange(compare_single_plot(compare_dat_elnet,1,'Elnet','Single 2'),compare_single_plot(compare_dat_Nnet,1,'Nnet','Single 2'),
             compare_single_plot(compare_dat_boost,1,'Boosting','Single 2'),compare_single_plot(compare_dat_tree,1,'tree','Single 2'),
             ncol=2)









pred_stacking <- lapply(1:length(dat_test), function(x) as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)


pred_stacking_wz <- lapply(1:length(dat_test), function(x) 
  as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)


pred_stacking_int <- lapply(1:length(dat_test), function(x) 
  as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))

predM <- lapply(dat_test,function(dat_test) predict(fitM,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))

pred_simp_avg <- lapply(1:length(dat_test), function(i) 
  sapply(1:nrow(dat_test[[1]]), function(y) mean(sapply(pred_test[[i]], function(x) x[y]))))

### Elastic Net
fit_single_1 <-  caret::train(y ~ X + X:X1, data=lst[[1]], method = "glmnet",trControl=con,preProc = c("center", "scale"))
pred_single_1 <- lapply(dat_test,function(dat_test) predict(fit_single_1,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))

fit_single_2 <-  caret::train(y ~ X + X:X1, data=lst[[2]], method = "glmnet",trControl=con,preProc = c("center", "scale"))
pred_single_2 <- lapply(dat_test,function(dat_test) predict(fit_single_2,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))

### Neural Net
fit_single_1 <- monmlp::monmlp.fit(x = as.matrix(lst[[1]][,names(lst[[1]])!='y']), y = as.matrix(lst[[1]]$y), hidden1 = 2, hidden2 = 2,
                                 bag = TRUE, iter.max = 500, iter.stopped = 10)
pred_single_1 <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fit_single_1))

fit_single_2 <- monmlp::monmlp.fit(x = as.matrix(lst[[2]][,names(lst[[1]])!='y']), y = as.matrix(lst[[2]]$y), hidden1 = 2, hidden2 = 2,
                                   bag = TRUE, iter.max = 500, iter.stopped = 10)
pred_single_2 <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fit_single_2))

## Boosting
fit_single_1 <- caret::train(y~., data=lst[[1]],method = "gbm",trControl=con,preProc = c("center", "scale"))
pred_single_1 <- lapply(dat_test,function(dat_test) predict(fit_single_1,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))

fit_single_2 <- caret::train(y~., data=lst[[2]],method = "gbm",trControl=con,preProc = c("center", "scale"))
pred_single_2 <- lapply(dat_test,function(dat_test) predict(fit_single_2,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))

## Tree Bag
fit_single_1 <- caret::train(y~., data=lst[[1]],method = "treebag",trControl=con,preProc = c("center", "scale"))
pred_single_1 <- lapply(dat_test,function(dat_test) predict(fit_single_1,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))

fit_single_2 <- caret::train(y~., data=lst[[2]],method = "treebag",trControl=con,preProc = c("center", "scale"))
pred_single_2 <- lapply(dat_test,function(dat_test) predict(fit_single_2,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))

#write.csv(dat_test_sig2_1,'~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/dat_test_sig2_1.csv')
#dat_test_sig3_1 <- read.csv('~/Desktop/master_project/base_case/Base_beta_sigma/beta_sigma_5s_10p/dat_test_sig3_1.csv',header=T)

compare_dat_tree <- rbind(data.frame(pred = unlist(pred_stacking),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking',500)),
                          data.frame(pred = unlist(pred_stacking_wz),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_wz',500)),
                          data.frame(pred = unlist(pred_stacking_int),X_value = do.call(rbind,dat_test)$X, Weight = rep('Stacking_int',500)),
                          data.frame(pred = unlist(pred_simp_avg),X_value = do.call(rbind,dat_test)$X, Weight = rep('Simp_Avg',500)),
                          data.frame(pred = unlist(predM),X_value = do.call(rbind,dat_test)$X, Weight = rep('Merge',500)),
                          data.frame(pred = unlist(pred_single_1),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 1',500)),
                          data.frame(pred = unlist(pred_single_2),X_value = do.call(rbind,dat_test)$X, Weight = rep('Single 2',500)),
                          data.frame(pred = do.call(rbind,dat_test)$y,X_value = do.call(rbind,dat_test)$X, Weight = rep('Actual',500))
                          )

write.csv(compare_dat_tree,'~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/compare_dat_tree.csv')

compare_dat_elnet <- read.csv('~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/compare_dat_elnet.csv')
compare_dat_boost <- read.csv('~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/compare_dat_boost.csv')
compare_dat_Nnet <- read.csv('~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/compare_dat_Nnet.csv')
compare_dat_tree <- read.csv('~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/compare_dat_tree.csv')
compare_plot <- function(dat,size,learner){
  dat <- dat[-which(dat$Weight %in% c('Single 1','Single 2','Stacking','Stacking_wz')),]
  p_plot <- ggplot(dat, aes(x=X_value, y=pred, shape=Weight, color=Weight)) +
    geom_point() +
    ggtitle(paste(learner,'Sigma = ', size)) +
    theme(axis.text.x = element_text(angle = 0,size=18),legend.text = element_text(size = 18), 
          legend.title = element_text(size = 20),axis.title.x = element_text(size = 18), 
          axis.title.y = element_text(size = 18), plot.title = element_text(size = 20))
  
  return(p_plot)
  
}

compare_plot(compare_dat_tree,1,'Treebag')

grid.arrange(compare_plot(compare_dat_elnet,1,'Elnet'),compare_plot(compare_dat_Nnet,1,'Nnet'),
             compare_plot(compare_dat_boost,1,'Boosting'),compare_plot(compare_dat_tree,1,'Treebag'),
             ncol=2)



get_single_dat <- function(dat,num){
  compare_single_dat <- rbind(data.frame(pred = dat[which(dat$Weight %in% 'Stacking'),]$pred,Single_pred = dat[which(dat$Weight %in% num),]$pred, Weight = rep('Stacking',500)),
                                data.frame(pred = dat[which(dat$Weight %in% 'Simp_Avg'),]$pred,Single_pred = dat[which(dat$Weight %in% num),]$pred, Weight = rep('Simp_Avg',500)),
                                data.frame(pred = dat[which(dat$Weight %in% 'Merge'),]$pred,Single_pred = dat[which(dat$Weight %in% num),]$pred, Weight = rep('Merge',500)),
                                data.frame(pred = dat[which(dat$Weight %in% 'Actual'),]$pred,Single_pred = dat[which(dat$Weight %in% num),]$pred, Weight = rep('Actual',500))
  )
  return(compare_single_dat)
}
compare_single_plot <- function(dat,size,learner,num){
  dat <- get_single_dat(dat,num)
  p_plot <- ggplot(dat, aes(x=Single_pred, y=pred, shape=Weight, color=Weight)) +
    geom_point(alpha=0.01) + guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) +
    ggtitle(paste(learner,'predictors vs. Single Predictor X Sigma = ', size)) +
    xlab('Single Predictor') + ylab('Predictors')
  return(p_plot)

}

#compare_single_plot(compare_dat_boost,1,'Boosting','Single 1')
grid.arrange(compare_single_plot(compare_dat_elnet,1,'Elnet','Single 1'),compare_single_plot(compare_dat_Nnet,1,'Nnet','Single 1'),
             compare_single_plot(compare_dat_boost,1,'Boosting','Single 1'),compare_single_plot(compare_dat_tree,1,'tree','Single 1'),
             ncol=2)
grid.arrange(compare_single_plot(compare_dat_elnet,1,'Elnet','Single 2'),compare_single_plot(compare_dat_Nnet,1,'Nnet','Single 2'),
             compare_single_plot(compare_dat_boost,1,'Boosting','Single 2'),compare_single_plot(compare_dat_tree,1,'tree','Single 2'),
             ncol=2)
save.image('~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/Single_pred.RData')
