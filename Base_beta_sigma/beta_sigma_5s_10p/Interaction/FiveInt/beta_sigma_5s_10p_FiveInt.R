library(MASS)
library(MCMCpack)
library(caret)
library(monmlp)

#load('~/Desktop/master_project/multisimstudy_test.Rdata')
#### Input:
# Parameter of Interest
# k- No. of studies
# n_k - sample size of each study of k: a vector indicating each length
# p - No. of covariates, assume equal firstly
# mu - k*p matrix of covariance mean
# SIGMA - p*p covariance matrix
# link function - x*beta in our case
# m - No. of mixing component
# Beta - m*p matrix of coefficient means
# sigma - m*p matrix of coefficient variance
# pi - m*1 vector of mxing probability
# sigma_error - noise term 
#### Output:
# List of k
# Each list with y_k,X_k,beta_k,C_k
MultiStudySim <- function(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta){
  x_vec_list <- lapply(1:k,function(x) mvrnorm(n = nk[x], mu_x[x,], SIG)) 
  
  Ck <- sample(1:m,k,replace = T,prob=pii)
  beta_vec_list <- lapply(1:k,function(z) apply(matrix(c(1:p),nrow=p,ncol=1),1,function(x) 
    rnorm(1,mu_beta[Ck[z],x],sigma_beta[Ck[z],x])))
  output_Y_list <- lapply(1:k,function(x) 
    x_vec_list[[x]]%*%beta_vec_list[[x]] + x_vec_list[[x]][,1]*x_vec_list[[x]][,2]+  x_vec_list[[x]][,1]*x_vec_list[[x]][,3] +  x_vec_list[[x]][,1]*x_vec_list[[x]][,4] + 
      x_vec_list[[x]][,1]*x_vec_list[[x]][,5] +  x_vec_list[[x]][,1]*x_vec_list[[x]][,6] + rnorm(nrow(x_vec_list[[x]]),0,1))
  
  final_data <- lapply(1:k,function(x) list(SimulatedOutput=data.frame(x_vec_list[[x]],y=output_Y_list[[x]],row.names = c(1:nrow(output_Y_list[[x]]))),
                                            BetaValue=beta_vec_list[[x]],Ck=Ck[x]))
  
  names(final_data) <- paste0('Study',c(1:k))
  return(final_data)
}

##### Simple Size Avg ####
simp_avg <- function(list1,dat_test){
  result <- rep(0,nrow(dat_test))
  for (i in 1:nrow(dat_test)){
    result[i] <- mean(sapply(list1, function(x) x[i]))
  }
  return(result)
}

###### Sample Size Average #####
#simp_size_avg <- function(list_final,list_ori,dat_train,dat_test){
#  result <- rep(0,nrow(dat_test))
#  for (i in 1:nrow(dat_test)){
#    result[i] <- sum(unlist(lapply(list_ori, nrow))/nrow(dat_train)* sapply(list_final, function(x) x[i]))
#  }
#  return(result)
#}

stack_weight_combine <- function(pred_train,dat_train){
  stack_train_x <- as.matrix.data.frame(do.call(cbind,pred_train))  
  #  stack_test_x <- as.matrix.data.frame(do.call(cbind,pred_test)) 
  nl <- nnls::nnls(stack_train_x, dat_train$y) 
  #  pred <- stack_test_x%*%coef(nl)
  #  pred_norm <- stack_test_x%*%(coef(nl)/sum(coef(nl)))
  return (coef(nl))
  
}

stack_weight_int <- function(pred_train,dat_train){
  stack_train_x <- as.matrix.data.frame(cbind(1,-1,do.call(cbind,pred_train)))
  nl <- nnls::nnls(stack_train_x, dat_train$y) 
  return (coef(nl))
}

stack_weight_wz <- function(pred_train,dat_train,nk_train){
  pred_train_wz <- data.frame(do.call(cbind,pred_train))
  index <- lapply(colnames(pred_train_wz),function(x) grep(paste0(x,'.'),rownames(pred_train_wz),fixed = T))
  for (i in 1:length(index)){
    pred_train_wz[c(index[[i]]),i] = 0
  }
  nl <- nnls::nnls(as.matrix.data.frame(pred_train_wz), dat_train$y)
  nl_wz <- sapply(1:length(coef(nl)), function(i) coef(nl)[i]*(1-nk_train/sum(nk_train))[i])
  return(nl_wz)
}

elnet_weight_int <- function(pred_train,dat_train){
  train.data <- as.matrix.data.frame(cbind(do.call(cbind,pred_train),dat_train$y))  
  colnames(train.data) <- c(1:(ncol(train.data)-1),'y')
  model <- train(
    y ~., data = train.data, method = "glmnet",
    trControl = trainControl("cv", number = 5),
    tuneLength = 10
  )
  nl_elnet_int <- coef(model$finalModel, model$bestTune$lambda)
  return(nl_elnet_int)
  
}

lasso_weight_noint <- function(pred_train,dat_train){
  x1 <- as.matrix.data.frame(do.call(cbind,pred_train))  
  y <- dat_train$y
  mod <- cv.glmnet(x1,y,alpha=1,intercept = F)
  nl_lasso_noint <- coef(mod)
  return(nl_lasso_noint)
}

ridge_weight_noint <- function(pred_train,dat_train){
  x1 <- as.matrix.data.frame(do.call(cbind,pred_train))  
  y <- dat_train$y
  mod <- cv.glmnet(x1,y,alpha=0,intercept = F)
  nl_ridge_noint <- coef(mod)
  return(nl_ridge_noint)
}

elnet_class_all <- function(lst,dat_train,dat_test,nk_train){
  # List: generate from multiStudySim
  # dat_train: one large dataframe (unlist List)
  # dat_testL List of another new data
  
  con <- caret::trainControl(method = 'cv',number=5,allowParallel = T)
  fit <- lapply(lst,function(lst) caret::train(y ~ ., data=lst, method = "glmnet",trControl=con,preProc = c("center", "scale")))
  #insamp_rmse_list <- lapply(1:length(fit),function(x) 
  #  Metrics::rmse(predict(fit[[x]],lst[[x]][,-which(names(lst[[x]]) %in% c('y'))]),lst[[x]]$y))
  
  fitM <- caret::train(y ~ ., data=dat_train, method = "glmnet",trControl=con,preProc = c("center", "scale"))
  
  pred_train <- lapply(fit, function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test <- lapply(1:length(dat_test),function(x) predict(fit,dat_test[[x]]  ))
  
  nl <- stack_weight_combine(pred_train,dat_train)
  stacking_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl,dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  nl_int <- stack_weight_int(pred_train,dat_train)
  stacking_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]),dat_test[[x]]$y))
  stacking_int_rmse_mean <- mean(unlist(stacking_int_rmse_list))
  
  nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
  stacking_wz_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz,dat_test[[x]]$y))
  stacking_wz_rmse_mean <- mean(unlist(stacking_wz_rmse_list))
  
  nl_elnet_int <- elnet_weight_int(pred_train, dat_train)
  stacking_elnet_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_elnet_int[-1]+nl_elnet_int[1],dat_test[[x]]$y))
  stacking_elnet_int_rmse_mean <- mean(unlist(stacking_elnet_int_rmse_list))
  
  nl_lasso_noint <- lasso_weight_noint(pred_train,dat_train)
  stacking_lasso_noint_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_lasso_noint[-1],dat_test[[x]]$y))
  stacking_lasso_noint_rmse_mean <- mean(unlist(stacking_lasso_noint_rmse_list))
  
  nl_ridge_noint <- ridge_weight_noint(pred_train,dat_train)
  stacking_ridge_noint_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_ridge_noint[-1],dat_test[[x]]$y))
  stacking_ridge_noint_rmse_mean <- mean(unlist(stacking_ridge_noint_rmse_list))
  
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
              elnet_int = stacking_elnet_int_rmse_mean, lasso_noint = stacking_lasso_noint_rmse_mean, ridge_noint = stacking_ridge_noint_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean))
}

nnet_class_all <- function(lst,dat_train,dat_test,nk_train){
  
  fit <- lapply(lst, function(lst) monmlp::monmlp.fit(x = as.matrix(lst[,names(lst)!='y']), y = as.matrix(lst$y), hidden1 = 2, hidden2 = 2,
                                                      bag = TRUE, iter.max = 500, iter.stopped = 10))
  
  fitM <- monmlp::monmlp.fit(x = as.matrix(dat_train[,names(dat_train)!='y']), y = as.matrix(dat_train$y), hidden1 = 2, hidden2 = 2,
                             bag = TRUE, iter.max = 500, iter.stopped = 10)
  
  pred_train <- lapply(fit, function(x) monmlp.predict(as.matrix(dat_train[,-which(names(dat_train) %in% c('y'))]),x))
  pred_test <- lapply(1:length(dat_test),function(x)  
    lapply(1:length(fit), function(y)  monmlp.predict(as.matrix(dat_test[[x]][,names(dat_test[[x]])!='y']),fit[[y]])) )
  
  nl <- stack_weight_combine(pred_train,dat_train)
  stacking_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl,dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  nl_int <- stack_weight_int(pred_train,dat_train)
  stacking_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]),dat_test[[x]]$y))
  stacking_int_rmse_mean <- mean(unlist(stacking_int_rmse_list))
  
  nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
  stacking_wz_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz,dat_test[[x]]$y))
  stacking_wz_rmse_mean <- mean(unlist(stacking_wz_rmse_list))
  
  nl_elnet_int <- elnet_weight_int(pred_train, dat_train)
  stacking_elnet_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_elnet_int[-1]+nl_elnet_int[1],dat_test[[x]]$y))
  stacking_elnet_int_rmse_mean <- mean(unlist(stacking_elnet_int_rmse_list))
  
  nl_lasso_noint <- lasso_weight_noint(pred_train,dat_train)
  stacking_lasso_noint_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_lasso_noint[-1],dat_test[[x]]$y))
  stacking_lasso_noint_rmse_mean <- mean(unlist(stacking_lasso_noint_rmse_list))
  
  nl_ridge_noint <- ridge_weight_noint(pred_train,dat_train)
  stacking_ridge_noint_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_ridge_noint[-1],dat_test[[x]]$y))
  stacking_ridge_noint_rmse_mean <- mean(unlist(stacking_ridge_noint_rmse_list))
  
  predM <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fitM))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  #insamp_rmse_mean <- mean(unlist(insamp_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
              elnet_int = stacking_elnet_int_rmse_mean, lasso_noint = stacking_lasso_noint_rmse_mean, ridge_noint = stacking_ridge_noint_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean))  
}

boost_class_all <- function(lst,dat_train,dat_test,nk_train){
  con <- caret::trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "gbm",trControl=con,preProc = c("center", "scale")))
  
  fitM <- caret::train(y~.,data=dat_train, method = "gbm",trControl=con,preProc = c("center", "scale"))
  
  pred_train <- lapply(fit, function(fit) predict(fit,dat_train))
  pred_test <- lapply(1:length(dat_test),function(x) predict(fit,dat_test[[x]]  ))
  
  
  nl <- stack_weight_combine(pred_train,dat_train)
  stacking_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl,dat_test[[x]]$y))
  
  nl_int <- stack_weight_int(pred_train,dat_train)
  stacking_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]),dat_test[[x]]$y))
  stacking_int_rmse_mean <- mean(unlist(stacking_int_rmse_list))
  
  nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
  stacking_wz_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz,dat_test[[x]]$y))
  stacking_wz_rmse_mean <- mean(unlist(stacking_wz_rmse_list))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  nl_elnet_int <- elnet_weight_int(pred_train, dat_train)
  stacking_elnet_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_elnet_int[-1]+nl_elnet_int[1],dat_test[[x]]$y))
  stacking_elnet_int_rmse_mean <- mean(unlist(stacking_elnet_int_rmse_list))
  
  nl_lasso_noint <- lasso_weight_noint(pred_train,dat_train)
  stacking_lasso_noint_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_lasso_noint[-1],dat_test[[x]]$y))
  stacking_lasso_noint_rmse_mean <- mean(unlist(stacking_lasso_noint_rmse_list))
  
  nl_ridge_noint <- ridge_weight_noint(pred_train,dat_train)
  stacking_ridge_noint_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_ridge_noint[-1],dat_test[[x]]$y))
  stacking_ridge_noint_rmse_mean <- mean(unlist(stacking_ridge_noint_rmse_list))
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,dat_test)) 
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  #insamp_rmse_mean <- mean(unlist(insamp_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
              elnet_int = stacking_elnet_int_rmse_mean, lasso_noint = stacking_lasso_noint_rmse_mean, ridge_noint = stacking_ridge_noint_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean))
}

rf_class_all <- function(lst,dat_train,dat_test,nk_train){
  con <- caret::trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "ranger",trControl=con,preProc = c("center", "scale")))
  #insamp_rmse_list <- lapply(1:length(fit),function(x) 
  #  Metrics::rmse(predict(fit[[x]],lst[[x]][,-which(names(lst[[x]]) %in% c('y'))]),lst[[x]]$y))
  
  fitM <- caret::train(y~., data=dat_train,method = "ranger",trControl=con,preProc = c("center", "scale"))
  
  pred_train <- lapply(fit, function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test <- lapply(1:length(dat_test),function(x) predict(fit,dat_test[[x]][,-which(names(dat_test[[x]]) %in% c('y'))]  ))
  
  nl <- stack_weight_combine(pred_train,dat_train)
  
  stacking_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl,dat_test[[x]]$y))
  
  nl_int <- stack_weight_int(pred_train,dat_train)
  stacking_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]),dat_test[[x]]$y))
  stacking_int_rmse_mean <- mean(unlist(stacking_int_rmse_list))
  
  nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
  stacking_wz_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz,dat_test[[x]]$y))
  stacking_wz_rmse_mean <- mean(unlist(stacking_wz_rmse_list))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  nl_elnet_int <- elnet_weight_int(pred_train, dat_train)
  stacking_elnet_int_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_elnet_int[-1]+nl_elnet_int[1],dat_test[[x]]$y))
  stacking_elnet_int_rmse_mean <- mean(unlist(stacking_elnet_int_rmse_list))
  
  nl_lasso_noint <- lasso_weight_noint(pred_train,dat_train)
  stacking_lasso_noint_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_lasso_noint[-1],dat_test[[x]]$y))
  stacking_lasso_noint_rmse_mean <- mean(unlist(stacking_lasso_noint_rmse_list))
  
  nl_ridge_noint <- ridge_weight_noint(pred_train,dat_train)
  stacking_ridge_noint_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_ridge_noint[-1],dat_test[[x]]$y))
  stacking_ridge_noint_rmse_mean <- mean(unlist(stacking_ridge_noint_rmse_list))
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  
 return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
             elnet_int = stacking_elnet_int_rmse_mean, lasso_noint = stacking_lasso_noint_rmse_mean, ridge_noint = stacking_ridge_noint_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean))  
}

get_input_dat <- function(list){
  maindat <- lapply(list,'[[',1)
  dat_train <- do.call(rbind,maindat)
  return(list(maindat_list=maindat,dat_train=dat_train))
}
get_rmse <- function(list,dat_train,dat_test,nk_train){
  elnet <- elnet_class_all(list,dat_train,dat_test,nk_train)
  nnet <- nnet_class_all(list,dat_train,dat_test,nk_train)
  boost <- boost_class_all(list,dat_train,dat_test,nk_train)
  rf <- rf_class_all(list,dat_train,dat_test,nk_train)
  #treebag <- treebag_class_all(list,dat_train,dat_test,nk_train)
  output <- list(elnet=elnet,nnet=nnet,boost=boost,rf=rf)
  
  return(output)
}

get_boxplot <- function(output){
  tt_elnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',1)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'elnet')
  tt_nnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',2)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'Nnet')
  tt_boost <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',3)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'boost')
  tt_rf <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',4)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'rf')
  colnames(tt_elnet)=colnames(tt_nnet)=colnames(tt_boost)=colnames(tt_rf)
  plotting <- data.frame(rbind(tt_elnet,tt_nnet,tt_boost,tt_rf))
  colnames(plotting)[3] <- 'Learners'
  plotting$measurement <- as.numeric(plotting$measurement)
  plotting$measurement <- log(plotting$measurement)
  re <- ggplot(data = plotting, aes(x=Learners, y=measurement,fill=condition)) + geom_boxplot()+ xlab("Learners") + ylab("log(RMSE)")
  return(re)
}
## Base case ###
k <- 5
nk <- rep(500,k)
nk_100 <- rep(100,k)
nk_train <- nk_100
#nk_150 <- rep(150,10)
#nk_200 <- rep(200,10)
#nk_250 <- rep(250,10)
#nk_300 <- rep(300,10)
#nk_350 <- rep(350,10)
#nk_400 <- rep(400,10)
#nk_450 <- rep(450,10)
p <- 10

SIG <- diag(x=1,nrow = p,ncol = p)
# SIG diagonal
# same mu x
# same mu beta
# 0 sigma beta
m <- 4
pii <- matrix(rep(0.25,4),nrow=m,ncol=1,byrow=T)
mu_x <- matrix(rep(0,k*p),nrow=k,ncol=p,byrow=T)

# Case X.1 : Similar mu_beta small sigma_beta
get_sep_unif <- function(a1,b1,a2,b2){
  k = sample(1:2,p,replace = T)
  re = rep(0,length(k))
  for (i in 1:length(k)){
    if (k[i] > 1){
      re[i] = runif(1,a1,b1)
    }else{re[i]= runif(1,a2,b2)}
  }
  return(re)
}
mu_beta_1.1 <- matrix(rep(get_sep_unif(2,3,-3,-2),m),nrow = m,ncol=p,byrow=T);mu_beta_1.1
# Make beta absolutely the same in every cluster
sigma_beta_0 <- matrix(rep(0,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_0.25 <- matrix(rep(0.25,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_0.5 <- matrix(rep(0.5,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_0.75 <- matrix(rep(0.75,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_1 <- matrix(rep(1,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_1.5 <- matrix(rep(1.5,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_2 <- matrix(rep(2,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_2.5 <- matrix(rep(2.5,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_3 <- matrix(rep(3,m*p),nrow = m,ncol=p,byrow = T)

test_0 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0)),recursive=F)
test_0.25 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0.25)),recursive=F)
test_0.5 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0.5)),recursive=F)
test_0.75 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0.75)),recursive=F)
test_1 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1)),recursive=F)
test_1.5 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.5)),recursive=F)
test_2 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_2)),recursive=F)
test_2.5 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_2.5)),recursive=F)
test_3 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_3)),recursive=F)

dat_input_0 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0))
dat_input_0.25 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0.25))
dat_input_0.5 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0.5))
dat_input_0.75 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0.75))

dat_input_1 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1))
dat_input_1.5 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.5))
dat_input_2 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_2))
dat_input_2.5 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_2.5))
dat_input_3 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_3))


#### Some Test ####
lst = get_input_dat(MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0))[[1]]
dat_train = get_input_dat(MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0))[[2]]
dat_test = lapply(test_0,'[[',1)
#test = elnet_class_all(lst,dat_train,dat_test,nk_train)
#testnnet = nnet_class_all(lst,dat_train,dat_test,nk_train)
#testboost = boost_class_all(lst,dat_train,dat_test,nk_train)
#testrf = rf_class_all(lst,dat_train,dat_test,nk_train)



## Do parallel ##
library(snow)
library(snowfall)
sfInit( parallel=TRUE, cpus=18 )

sfLibrary(MASS)
sfLibrary(MCMCpack)
sfLibrary(caret)
sfLibrary(monmlp)
sfLibrary(glmnet)
sfExportAll( )

case_0 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_0[[x]])[[1]],
                                                              get_input_dat(dat_input_0[[x]])[[2]],
                                                              lapply(test_0, '[[',1),nk_train))
saveRDS(case_0, file = "~/MultiSim/result/beta_sigma_5s_10p/Interaction/FiveInt/case_0.rds")

case_0.25 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_0.25[[x]])[[1]],
                                                                 get_input_dat(dat_input_0.25[[x]])[[2]],
                                                                 lapply(test_0.25, '[[',1),nk_train))
saveRDS(case_0.25, file = "~/MultiSim/result/beta_sigma_5s_10p/Interaction/FiveInt/case_0.25.rds")

case_0.5 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_0.5[[x]])[[1]],
                                                                get_input_dat(dat_input_0.5[[x]])[[2]],
                                                                lapply(test_0.5, '[[',1),nk_train))
saveRDS(case_0.5, file = "~/MultiSim/result/beta_sigma_5s_10p/Interaction/FiveInt/case_0.5.rds")

case_0.75 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_0.75[[x]])[[1]],
                                                                 get_input_dat(dat_input_0.75[[x]])[[2]],
                                                                 lapply(test_0.75, '[[',1),nk_train))
saveRDS(case_0.75, file = "~/MultiSim/result/beta_sigma_5s_10p/Interaction/FiveInt/case_0.75.rds")

case_1 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_1[[x]])[[1]],
                                                              get_input_dat(dat_input_1[[x]])[[2]],
                                                              lapply(test_1, '[[',1),nk_train))
saveRDS(case_1, file = "~/MultiSim/result/beta_sigma_5s_10p/Interaction/FiveInt/case_1.rds")

case_1.5 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_1.5[[x]])[[1]],
                                                                get_input_dat(dat_input_1.5[[x]])[[2]],
                                                                lapply(test_1.5, '[[',1),nk_train))
saveRDS(case_1.5, file = "~/MultiSim/result/beta_sigma_5s_10p/Interaction/FiveInt/case_1.5.rds")

case_2 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_2[[x]])[[1]],
                                                              get_input_dat(dat_input_2[[x]])[[2]],
                                                              lapply(test_2, '[[',1),nk_train))
saveRDS(case_2, file = "~/MultiSim/result/beta_sigma_5s_10p/Interaction/FiveInt/case_2.rds")

case_2.5 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_2.5[[x]])[[1]],
                                                                get_input_dat(dat_input_2.5[[x]])[[2]],
                                                                lapply(test_2.5, '[[',1),nk_train))
saveRDS(case_2.5, file = "~/MultiSim/result/beta_sigma_5s_10p/Interaction/FiveInt/case_2.5.rds")

case_3 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_3[[x]])[[1]],
                                                              get_input_dat(dat_input_3[[x]])[[2]],
                                                              lapply(test_3, '[[',1),nk_train))
saveRDS(case_3, file = "~/MultiSim/result/beta_sigma_5s_10p/Interaction/FiveInt/case_3.rds")



