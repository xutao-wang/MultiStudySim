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
    x_vec_list[[x]]%*%beta_vec_list[[x]] + rnorm(nrow(x_vec_list[[x]]),0,1))
  
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
simp_size_avg <- function(list_final,list_ori,dat_train,dat_test){
  result <- rep(0,nrow(dat_test))
  for (i in 1:nrow(dat_test)){
    result[i] <- sum(unlist(lapply(list_ori, nrow))/nrow(dat_train)* sapply(list_final, function(x) x[i]))
  }
  return(result)
}

stack_weight_combine <- function(pred_train,dat_train){
  stack_train_x <- as.matrix.data.frame(do.call(cbind,pred_train))  
  #  stack_test_x <- as.matrix.data.frame(do.call(cbind,pred_test)) 
  nl <- nnls::nnls(stack_train_x, dat_train$y) 
  #  pred <- stack_test_x%*%coef(nl)
  #  pred_norm <- stack_test_x%*%(coef(nl)/sum(coef(nl)))
  return (coef(nl))
  
}

elnet_class_all <- function(lst,dat_train,dat_test){
  # List: generate from multiStudySim
  # dat_train: one large dataframe (unlist List)
  # dat_testL List of another new data
  
  con <- caret::trainControl(method = 'cv',number=5,allowParallel = T)
  fit <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "glmnet",trControl=con,preProc = c("center", "scale")))
  insamp_rmse_list <- lapply(1:length(fit),function(x) 
    Metrics::rmse(predict(fit[[x]],lst[[x]][,-which(names(lst[[x]]) %in% c('y'))]),lst[[x]]$y))
  
  fitM <- caret::train(y~., data=dat_train, method = "glmnet",trControl=con,preProc = c("center", "scale"))
  list_fit <- list(MultiFit=fit,MergeFit=fitM)
  
  pred_train <- lapply(list_fit[[1]], function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test <- lapply(1:length(dat_test),function(x) predict(list_fit[[1]],dat_test[[x]][,-which(names(dat_test[[x]]) %in% c('y'))]  ))
  
  nl <- stack_weight_combine(pred_train,dat_train)
  nl_norm <- nl/sum(nl)
  stacking_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl,dat_test[[x]]$y))
  stacking_norm_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_norm,dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  
  predM <- lapply(dat_test,function(dat_test) predict(list_fit[[2]],as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  stacking_norm_rmse_mean <- mean(unlist(stacking_norm_rmse_list))
  insamp_rmse_mean <- mean(unlist(insamp_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_norm=stacking_norm_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,insamp_avg=insamp_rmse_mean,
              nl_coef=nl,nl_coef_norm=nl_norm))
}

#### Some Test ####
#a =get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.1))
#b = lapply(test_1.1, '[[',1)

#test = elnet_class_all(a[[1]],a[[2]],b)
#testnnet = nnet_class_all(a[[1]],a[[2]],b)
#testboost = boost_class_all(a[[1]],a[[2]],b)
#testrf = rf_class_all(a[[1]],a[[2]],b)


nnet_class_all <- function(lst,dat_train,dat_test){
  
  fit <- lapply(lst, function(lst) monmlp::monmlp.fit(x = as.matrix(lst[,names(lst)!='y']), y = as.matrix(lst$y), hidden1 = 2, hidden2 = 2,
                                                      bag = TRUE, iter.max = 500, iter.stopped = 10))
  insamp_rmse_list <- lapply(1:length(fit),function(x) 
    Metrics::rmse(monmlp.predict(as.matrix(lst[[x]][,names(lst[[x]])!='y']), fit[[x]]),lst[[x]]$y))
  
  fitM <- monmlp::monmlp.fit(x = as.matrix(dat_train[,names(dat_train)!='y']), y = as.matrix(dat_train$y), hidden1 = 2, hidden2 = 2,
                             bag = TRUE, iter.max = 500, iter.stopped = 10)
  
  pred_train <- lapply(fit, function(x) monmlp.predict(as.matrix(dat_train[,-which(names(dat_train) %in% c('y'))]),x))
  pred_test <- lapply(1:length(dat_test),function(x)  
    lapply(1:length(fit), function(y)  monmlp.predict(as.matrix(dat_test[[x]][,names(dat_test[[x]])!='y']),fit[[y]])) )
  
  nl <- stack_weight_combine(pred_train,dat_train)
  nl_norm <- nl/sum(nl)
  stacking_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl,dat_test[[x]]$y))
  stacking_norm_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_norm,dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  predM <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fitM))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  stacking_norm_rmse_mean <- mean(unlist(stacking_norm_rmse_list))
  insamp_rmse_mean <- mean(unlist(insamp_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_norm=stacking_norm_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,insamp_avg=insamp_rmse_mean,
              nl_coef=nl,nl_coef_norm=nl_norm))  
}

boost_class_all <- function(lst,dat_train,dat_test){
  con <- caret::trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "gbm",trControl=con,preProc = c("center", "scale")))
  insamp_rmse_list <- lapply(1:length(fit),function(x) 
    Metrics::rmse(predict(fit[[x]],lst[[x]][,-which(names(lst[[x]]) %in% c('y'))]),lst[[x]]$y))
  
  fitM <- caret::train(y~.,data=dat_train, method = "gbm",trControl=con,preProc = c("center", "scale"))
  list_fit <- list(MultiFit=fit,MergeFit=fitM)
  
  pred_train <- lapply(list_fit[[1]], function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test <- lapply(1:length(dat_test),function(x) predict(list_fit[[1]],dat_test[[x]][,-which(names(dat_test[[x]]) %in% c('y'))]  ))
  
  
  nl <- stack_weight_combine(pred_train,dat_train)
  nl_norm <- nl/sum(nl)
  stacking_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl,dat_test[[x]]$y))
  stacking_norm_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_norm,dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  
  predM <- lapply(dat_test,function(dat_test) predict(list_fit[[2]],as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  stacking_norm_rmse_mean <- mean(unlist(stacking_norm_rmse_list))
  insamp_rmse_mean <- mean(unlist(insamp_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_norm=stacking_norm_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,insamp_avg=insamp_rmse_mean,
              nl_coef=nl,nl_coef_norm=nl_norm))
}

rf_class_all <- function(lst,dat_train,dat_test){
  con <- caret::trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "ranger",trControl=con,preProc = c("center", "scale")))
  insamp_rmse_list <- lapply(1:length(fit),function(x) 
    Metrics::rmse(predict(fit[[x]],lst[[x]][,-which(names(lst[[x]]) %in% c('y'))]),lst[[x]]$y))
  
  fitM <- caret::train(y~., data=dat_train,method = "ranger",trControl=con,preProc = c("center", "scale"))
  list_fit <- list(MultiFit=fit,MergeFit=fitM)
  
  pred_train <- lapply(list_fit[[1]], function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test <- lapply(1:length(dat_test),function(x) predict(list_fit[[1]],dat_test[[x]][,-which(names(dat_test[[x]]) %in% c('y'))]  ))
  
  nl <- stack_weight_combine(pred_train,dat_train)
  nl_norm <- nl/sum(nl)
  stacking_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl,dat_test[[x]]$y))
  stacking_norm_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_norm,dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  
  predM <- lapply(dat_test,function(dat_test) predict(list_fit[[2]],as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  
  stacking_norm_rmse_mean <- mean(unlist(stacking_norm_rmse_list))
  insamp_rmse_mean <- mean(unlist(insamp_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_norm=stacking_norm_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,insamp_avg=insamp_rmse_mean,
              nl_coef=nl,nl_coef_norm=nl_norm))  
}
get_input_dat <- function(list){
  maindat <- lapply(list,'[[',1)
  dat_train <- do.call(rbind,maindat)
  return(list(maindat_list=maindat,dat_train=dat_train))
}
get_rmse <- function(list,dat_train,dat_test){
  elnet <- elnet_class_all(list,dat_train,dat_test)
  nnet <- nnet_class_all(list,dat_train,dat_test)
  boost <- boost_class_all(list,dat_train,dat_test)
  rf <- rf_class_all(list,dat_train,dat_test)
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
k <- 10
nk <- rep(500,10)
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
# Make beta absolutrly the same in every cluster

sigma_beta_1.1 <- matrix(rep(0,m*p),nrow = m,ncol=p,byrow = T)

#train_1.1 <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.1)
test_1.1 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.1)),recursive=F)

dat_input <- lapply(1:50,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.1))
## Do parallel ##
library(snow)
library(snowfall)
sfInit( parallel=TRUE, cpus=12 )

sfLibrary(MASS)
sfLibrary(MCMCpack)
sfLibrary(caret)
sfLibrary(readxl)
sfLibrary(monmlp)
sfExportAll( )

case_base <- sfLapply(1:length(dat_input),function(x) get_rmse(get_input_dat(dat_input[[x]])[[1]],
                                              get_input_dat(dat_input[[x]])[[2]],
                                              lapply(test_1.1, '[[',1)))
case_base <- readRDS('~/Desktop/master_project/base_case/case_base.RDS')
p1.1 <- get_boxplot(case_base)

p1.1 <- p1.1 + ggtitle("Similar mu_beta, small sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF','yellow'))
p1.1

get_boxplot_2 <- function(output){
  #tt_elnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',1)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'elnet')
  #tt_nnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',2)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'Nnet')
  tt_boost <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',3)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'boost')
  tt_rf <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',4)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'rf')
  colnames(tt_boost)=colnames(tt_rf)
  plotting <- data.frame(rbind(tt_boost,tt_rf))
  colnames(plotting)[3] <- 'Learners'
  plotting$measurement <- as.numeric(plotting$measurement)
  plotting$measurement <- log(plotting$measurement)
  re <- ggplot(data = plotting, aes(x=Learners, y=measurement,fill=condition)) + geom_boxplot()+ xlab("Learners") + ylab("log(RMSE)")
  return(re)
}
get_boxplot_3 <- function(output){
  tt_elnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',1)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'elnet')
  tt_nnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',2)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'Nnet')
  #tt_boost <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',3)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'boost')
  #tt_rf <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',4)))[,-c(5:7)], condition, measurement, stacking:merge, factor_key=TRUE),'rf')
  colnames(tt_elnet)=colnames(tt_nnet)
  plotting <- data.frame(rbind(tt_elnet,tt_nnet))
  colnames(plotting)[3] <- 'Learners'
  plotting$measurement <- as.numeric(plotting$measurement)
  plotting$measurement <- log(plotting$measurement)
  re <- ggplot(data = plotting, aes(x=Learners, y=measurement,fill=condition)) + geom_boxplot()+ xlab("Learners") + ylab("log(RMSE)")
  return(re)
}

p1 <- get_boxplot(case_base)
p1 <- p1 + ggtitle("Same mu_beta, 0 sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF','yellow'))

p2 <- get_boxplot_2(case_base)
p2 <- p2 + ggtitle("Same mu_beta, 0 sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF','yellow'))

p3 <- get_boxplot_3(case_base)
p3 <- p3 + ggtitle("Same mu_beta, 0 sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF','yellow'))

library(grid)
library(gridExtra)
grid.arrange(p1,p2,p3,ncol=2,top=textGrob("Same mu_x",gp=gpar(fontsize=20,font=3)))

saveRDS(case_base, file = "~/MultiSim/result/case_base.rds")
