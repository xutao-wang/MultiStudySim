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
nk_100 <- rep(100,10)
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
# Make beta absolutrly the same in every cluster
sigma_beta_0 <- matrix(rep(0,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_1 <- matrix(rep(1,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_3 <- matrix(rep(3,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_5 <- matrix(rep(5,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_7 <- matrix(rep(7,m*p),nrow = m,ncol=p,byrow = T)
sigma_beta_9 <- matrix(rep(9,m*p),nrow = m,ncol=p,byrow = T)

test_1 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1)),recursive=F)
test_3 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_3)),recursive=F)
test_5 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_5)),recursive=F)
test_7 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_7)),recursive=F)
test_9 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_9)),recursive=F)

dat_input_1 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1))
dat_input_3 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_3))
dat_input_5 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_5))
dat_input_7 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_7))
dat_input_9 <- lapply(1:50,function(x) MultiStudySim(k,nk_100,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_9))
## Do parallel ##
library(snow)
library(snowfall)
sfInit( parallel=TRUE, cpus=8 )

sfLibrary(MASS)
sfLibrary(MCMCpack)
sfLibrary(caret)
sfLibrary(monmlp)
sfExportAll( )

case_1 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_1[[x]])[[1]],
                                                                  get_input_dat(dat_input_1[[x]])[[2]],
                                                                  lapply(test_1, '[[',1)))
case_3 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_3[[x]])[[1]],
                                                                  get_input_dat(dat_input_3[[x]])[[2]],
                                                                  lapply(test_3, '[[',1)))
case_5 <- sfLapply(1:length(dat_input_5),function(x) get_rmse(get_input_dat(dat_input_5[[x]])[[1]],
                                                                  get_input_dat(dat_input_5[[x]])[[2]],
                                                                  lapply(test_5, '[[',1)))
case_7 <- sfLapply(1:length(dat_input_7),function(x) get_rmse(get_input_dat(dat_input_7[[x]])[[1]],
                                                                  get_input_dat(dat_input_7[[x]])[[2]],
                                                                  lapply(test_7, '[[',1)))
case_9 <- sfLapply(1:length(dat_input_9),function(x) get_rmse(get_input_dat(dat_input_9[[x]])[[1]],
                                                                  get_input_dat(dat_input_9[[x]])[[2]],
                                                                  lapply(test_9, '[[',1)))

#saveRDS(case_1, file = "~/MultiSim/result/case_1.rds")
#saveRDS(case_3, file = "~/MultiSim/result/case_3.rds")
#saveRDS(case_5, file = "~/MultiSim/result/case_5.rds")
#saveRDS(case_7, file = "~/MultiSim/result/case_7.rds")
#saveRDS(case_9, file = "~/MultiSim/result/case_9.rds")
case_0 <- readRDS('~/Desktop/master_project/base_case/base_sample_size/case_100.RDS')
case_1 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/case_1.RDS')
case_3 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/case_3.RDS')
case_5 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/case_5.RDS')
case_7 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/case_7.RDS')
case_9 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/case_9.RDS')

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  #datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  #ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  #datac$ci <- datac$se * ciMult
  
  return(datac)
}

get_long_data <- function(output,size){
  tt_elnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',1)))[,-c(5:7)], weighting, measurement, stacking:merge, factor_key=TRUE),'elnet')
  tt_nnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',2)))[,-c(5:7)], weighting, measurement, stacking:merge, factor_key=TRUE),'Nnet')
  tt_boost <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',3)))[,-c(5:7)], weighting, measurement, stacking:merge, factor_key=TRUE),'boost')
  tt_rf <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',4)))[,-c(5:7)], weighting, measurement, stacking:merge, factor_key=TRUE),'rf')
  colnames(tt_elnet)=colnames(tt_nnet)=colnames(tt_boost)=colnames(tt_rf)
  plotting <- data.frame(rbind(tt_elnet,tt_nnet,tt_boost,tt_rf))
  colnames(plotting)[3] <- 'Learners'
  plotting$measurement <- as.numeric(plotting$measurement)
  plotting$measurement <- log(plotting$measurement)
  plotting <- cbind(plotting,size)
  names(plotting)[4] <- 'Size'
  #re <- ggplot(data = plotting, aes(x=Learners, y=measurement,fill=condition)) + geom_boxplot()+ xlab("Learners") + ylab("log(RMSE)")
  return(plotting)
}

dat_0 <- get_long_data(case_0,0)
dat_1 <- get_long_data(case_1,1)
dat_3 <- get_long_data(case_3,3)
dat_5 <- get_long_data(case_5,5)
dat_7 <- get_long_data(case_7,7)
dat_9 <- get_long_data(case_9,9)

dat_total <- rbind(dat_0,dat_1,dat_3,dat_5,dat_7,dat_9)

get_learner_plotting <- function(dat_total,title,xlab){
  
  dat_elnet <- dat_total[which(dat_total$Learners %in% 'elnet'),names(dat_total)!='Learners']
  dat_elnet_plot <- summarySE(dat_elnet, measurevar="measurement", groupvars=c("weighting","Size"))
  pd <- position_dodge(1)
  p_elnet <- ggplot(dat_elnet_plot, aes(x=Size, y=measurement, color=weighting)) + ylab("log(RMSE)") +
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.3,position = pd) +
    geom_line(position = pd) + xlab(xlab) +
    geom_point(position = pd) + ggtitle("Elnet") +
    scale_x_continuous(breaks = dat_elnet_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  
  
  dat_nnet <- dat_total[which(dat_total$Learners %in% 'Nnet'),names(dat_total)!='Learners']
  dat_nnet_plot <- summarySE(dat_nnet, measurevar="measurement", groupvars=c("weighting","Size"))
  p_nnet <- ggplot(dat_nnet_plot, aes(x=Size, y=measurement, color=weighting)) + ylab("log(RMSE)") +
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.3,position = pd) +
    geom_line(position = pd) + xlab(xlab) +
    geom_point(position = pd) + ggtitle("Nnet") +
    scale_x_continuous(breaks = dat_nnet_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  
  
  dat_boost <- dat_total[which(dat_total$Learners %in% 'boost'),names(dat_total)!='Learners']
  dat_boost_plot <- summarySE(dat_boost, measurevar="measurement", groupvars=c("weighting","Size"))
  p_boost <- ggplot(dat_boost_plot, aes(x=Size, y=measurement, color=weighting)) + ylab("log(RMSE)") +
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.3,position = pd) +
    geom_line(position = pd) + xlab(xlab) +
    geom_point(position = pd) + ggtitle("Boosting") +
    scale_x_continuous(breaks = dat_boost_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  
  
  dat_rf <- dat_total[which(dat_total$Learners %in% 'rf'),names(dat_total)!='Learners']
  dat_rf_plot <- summarySE(dat_rf, measurevar="measurement", groupvars=c("weighting","Size"))
  p_rf <- ggplot(dat_rf_plot, aes(x=Size, y=measurement, color=weighting)) + ylab("log(RMSE)") +
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.3,position = pd) +
    geom_line(position = pd) + xlab(xlab) +
    geom_point(position = pd) + ggtitle("RF") +
    scale_x_continuous(breaks = dat_rf_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  
  
  library(grid)
  return(grid.arrange(p_elnet,p_nnet,p_boost,p_rf,ncol=2,top=textGrob(title,gp=gpar(fontsize=20,font=3))))
  
}

get_learner_plotting(dat_total,'Variance of Beta','Beta Sigma')

par(mfrow=c(3,2))
hist(unlist(lapply(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0),'[[',2)),density=T,main='Beta Sigma = 0',xlab='Beta Value')
hist(unlist(lapply(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1),'[[',2)),main='Beta Sigma = 1',xlab='Beta Value')
hist(unlist(lapply(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_3),'[[',2)),main='Beta Sigma = 3',xlab='Beta Value')
hist(unlist(lapply(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_5),'[[',2)),main='Beta Sigma = 5',xlab='Beta Value')
hist(unlist(lapply(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_7),'[[',2)),main='Beta Sigma = 7',xlab='Beta Value')
hist(unlist(lapply(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_9),'[[',2)),main='Beta Sigma = 9',xlab='Beta Value')
dev.off()

get_plot_ratio <- function(dat_total,learner){
  dat_plot <- summarySE(dat_total[which(dat_total$Learners %in% learner),names(dat_total)!='Learners'], measurevar="measurement", groupvars=c("weighting","Size"))
  dat_plot_ratio <- dat_elnet_plot[-which(dat_plot$weighting %in% 'stacking'),]
  base <- dat_plot$measurement[1:6]
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'stacking_norm'),]$measurement <- dat_plot[which(dat_plot$weighting %in% 'stacking_norm'),]$measurement/base
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'simp_avg'),]$measurement <- dat_plot[which(dat_elnet_plot$weighting %in% 'simp_avg'),]$measurement/base
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'merge'),]$measurement <- dat_plot[which(dat_elnet_plot$weighting %in% 'merge'),]$measurement/base
  pd <- position_dodge(0)
  p_ratio <- ggplot(dat_plot_ratio, aes(x=Size, y=measurement, color=weighting)) +
    #geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.3,position = pd) +
    geom_line(position = pd) + xlab('Beta Sigma') + ylab("log(RMSE)") +
    geom_point(position = pd) + ggtitle(learner) +
    scale_x_continuous(breaks = dat_plot_ratio$Size) + theme(axis.text.x = element_text(angle = 0))
  return(p_ratio)
}
grid.arrange(get_plot_ratio(dat_total,'elnet'),get_plot_ratio(dat_total,'Nnet'),get_plot_ratio(dat_total,'boost'),get_plot_ratio(dat_total,'rf'),
             ncol=2,top=textGrob('Beta Sigma Sampe Size 100',gp=gpar(fontsize=20,font=3)))
