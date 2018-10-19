library(MASS)
library(MCMCpack)
library(caret)
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

############## Demo--Input ##############
k <- 10
nk <- rep(100,k) # No. of samples in each study should be larger than 1
p <- 3
m <- 4
mu_x <- matrix(runif(k*p,1,20),nrow=k,ncol=p,byrow=T)
SIG <- MCMCpack::riwish(v = p+1, S = diag(1, p,p))
pii <- matrix(runif(m,0,1),nrow=m,ncol=1,byrow=T)

mu_beta <- matrix(runif(m*p,1,20),nrow = m,ncol=p,byrow=T)
sigma_beta <- matrix(runif(m*p,1,1.1),nrow = m,ncol=p,byrow = T) # Use uniform distribution

final <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta) #lapply(final,"[[", 3)


########### Simple Average #########
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

###### Stacked Weight with non-negaive least square #####
stack_weight_con <- function(pred_train,pred_test,dat_train){
  
  stack_train_x <- as.matrix.data.frame(do.call(cbind,pred_train))  
  stack_test_x <- as.matrix.data.frame(do.call(cbind,pred_test)) 
  nl <- nnls::nnls(stack_train_x, dat_train$y) 
  pred <- stack_test_x%*%coef(nl)
  return(pred)
  
}
stack_weight_norm <- function(pred_train,pred_test,dat_train){
  
  stack_train_x <- as.matrix.data.frame(do.call(cbind,pred_train))  
  stack_test_x <- as.matrix.data.frame(do.call(cbind,pred_test)) 
  nl <- nnls::nnls(stack_train_x, dat_train$y) 
  pred_norm <- stack_test_x%*%(coef(nl)/sum(coef(nl)))
  return(pred_norm)
  
}

get_input_dat <- function(list){
  maindat <- lapply(list,'[[',1)
  dat_train <- do.call(rbind,maindat)
  return(list(maindat_list=maindat,dat_train=dat_train))
}

elnet_class_all <- function(list,dat_train,dat_test){
  # List: generate from multiStudySim
  # dat_train: one large dataframe (unlist List)
  # dat_testL List of another new data

  con <- caret::trainControl(method = 'cv',number=5,allowParallel = T)
  fit <- lapply(list,function(list) caret::train(y~., data=list, method = "glmnet",trControl=con,preProc = c("center", "scale")))
  fitM <- caret::train(y~., data=dat_train, method = "glmnet",trControl=con,preProc = c("center", "scale"))
  list_fit <- list(MultiFit=fit,MergeFit=fitM)
  
  pred_train <- lapply(list_fit[[1]], function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test <- lapply(1:length(dat_test),function(x) predict(list_fit[[1]],dat_test[[x]][,-which(names(dat_test[[x]]) %in% c('y'))]  ))
  stack_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(stack_weight_con(pred_train,pred_test[[x]],dat_train),dat_test[[x]]$y))
  stack_norm_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(stack_weight_norm(pred_train,pred_test[[x]],dat_train),dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  
  predM <- lapply(dat_test,function(dat_test) predict(list_fit[[2]],as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))

  stacking_rmse_mean <- mean(unlist(stack_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  
  stacking_norm_rmse_mean <- mean(unlist(stack_norm_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_norm=stacking_norm_rmse_mean,simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean))
}


#a = get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.1))
#dat_test = lapply(test_1.1, '[[',1)
#dat_train = a[[2]]

nnet_class_all <- function(list,dat_train,dat_test){
  con <- trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(list, function(list) caret::train(y~., data=list, method = "nnet",trace=TRUE, maxit=1000, linout = 1,trControl=con,preProc = c("center", "scale")))
  fitM <- caret::train(y~., data=dat_train, method = "nnet",trace=TRUE, maxit=1000, linout = 1,trControl=con,preProc = c("center", "scale"))
  list_fit <- list(MultiFit=fit,MergeFit=fitM)
  
  pred_train <- lapply(list_fit[[1]], function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test <- lapply(1:length(dat_test),function(x) predict(list_fit[[1]],dat_test[[x]][,-which(names(dat_test[[x]]) %in% c('y'))]  ))
  stack_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(stack_weight_con(pred_train,pred_test[[x]],dat_train),dat_test[[x]]$y))
  stack_norm_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(stack_weight_norm(pred_train,pred_test[[x]],dat_train),dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  predM <- lapply(dat_test,function(dat_test) predict(list_fit[[2]],as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stack_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  stacking_norm_rmse_mean <- mean(unlist(stack_norm_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_norm=stacking_norm_rmse_mean,simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean))
  
}

boost_class_all <- function(list,dat_train,dat_test){
  con <- trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(list,function(list) caret::train(y~., data=list, method = "gbm",trControl=con,preProc = c("center", "scale")))
  fitM <- caret::train(y~.,data=dat_train, method = "gbm",trControl=con,preProc = c("center", "scale"))
  list_fit <- list(MultiFit=fit,MergeFit=fitM)
  
  pred_train <- lapply(list_fit[[1]], function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test <- lapply(1:length(dat_test),function(x) predict(list_fit[[1]],dat_test[[x]][,-which(names(dat_test[[x]]) %in% c('y'))]  ))
  stack_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(stack_weight_con(pred_train,pred_test[[x]],dat_train),dat_test[[x]]$y))
  stack_norm_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(stack_weight_norm(pred_train,pred_test[[x]],dat_train),dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  
  predM <- lapply(dat_test,function(dat_test) predict(list_fit[[2]],as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stack_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  stacking_norm_rmse_mean <- mean(unlist(stack_norm_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_norm=stacking_norm_rmse_mean,simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean))
}

rf_class_all <- function(list,dat_train,dat_test){
  con <- trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(list,function(list) caret::train(y~., data=list, method = "ranger",trControl=con,preProc = c("center", "scale")))
  fitM <- caret::train(y~., data=dat_train,method = "ranger",trControl=con,preProc = c("center", "scale"))
  list_fit <- list(MultiFit=fit,MergeFit=fitM)
  
  pred_train <- lapply(list_fit[[1]], function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
  pred_test <- lapply(1:length(dat_test),function(x) predict(list_fit[[1]],dat_test[[x]][,-which(names(dat_test[[x]]) %in% c('y'))]  ))
  stack_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(stack_weight_con(pred_train,pred_test[[x]],dat_train),dat_test[[x]]$y))
  stack_norm_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(stack_weight_norm(pred_train,pred_test[[x]],dat_train),dat_test[[x]]$y))
  
  simp_avg_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(simp_avg(pred_test[[x]],dat_test[[x]]),dat_test[[x]]$y))
  
  
  predM <- lapply(dat_test,function(dat_test) predict(list_fit[[2]],as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stack_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  
  stacking_norm_rmse_mean <- mean(unlist(stack_norm_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_norm=stacking_norm_rmse_mean,simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean))
}



get_rmse <- function(list,dat_train,dat_test){
  elnet <- elnet_class_all(list,dat_train,dat_test)
  nnet <- nnet_class_all(list,dat_train,dat_test)
  boost <- boost_class_all(list,dat_train,dat_test)
  rf <- rf_class_all(list,dat_train,dat_test)
  output <- list(elnet=elnet,nnet=nnet,boost=boost,rf=rf)
  return(output)
}

##### Make Boxplot #####
library(ggplot2)
library(grid)
library(gridExtra)
get_boxplot <- function(output){
  tt_elnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',1))), condition, measurement, stacking:merge, factor_key=TRUE),'elnet')
  tt_nnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',2))), condition, measurement, stacking:merge, factor_key=TRUE),'Nnet')
  tt_boost <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',3))), condition, measurement, stacking:merge, factor_key=TRUE),'boost')
  tt_rf <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',4))), condition, measurement, stacking:merge, factor_key=TRUE),'rf')
  colnames(tt_elnet)=colnames(tt_nnet)=colnames(tt_boost)=colnames(tt_rf)
  plotting <- data.frame(rbind(tt_elnet,tt_nnet,tt_boost,tt_rf))
  colnames(plotting)[3] <- 'Learners'
  plotting$measurement <- as.numeric(plotting$measurement)
  re <- ggplot(data = plotting, aes(x=Learners, y=measurement,fill=condition)) + geom_boxplot()+ xlab("Learners") + ylab("RMSE")
  return(re)
}

####### Test with different cases #####
# Case 1: large variation in covariates(X) mean small variance and cov
#Case 2: Small variation in covariates(X) mean small variance and cov

k <- 10
nk <- rep(500,10)
p <- 10
set.seed(123)
SIG <- riwish(v = p+1, S = diag(1, p,p))

m <- 4
set.seed(123)
pii <- matrix(rep(0.25,4),nrow=m,ncol=1,byrow=T);pii
set.seed(123)
mu_x <- matrix(runif(k*p,1,20),nrow=k,ncol=p,byrow=T);mu_x
set.seed(123)
mu_x_2 <- matrix(runif(k*p,5,10),nrow=k,ncol=p,byrow=T);mu_x_2

# Case X.1 : Similar mu_beta small sigma_beta
set.seed(123)
mu_beta_1.1 <- matrix(runif(m*p,2,5),nrow = m,ncol=p,byrow=T)
set.seed(123)
sigma_beta_1.1 <- matrix(runif(m*p,0,1),nrow = m,ncol=p,byrow = T) # Use uniform distribution

#train_1.1 <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.1)
test_1.1 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.1)),recursive=F)

#train_2.1 <- MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.1,sigma_beta_1.1)
test_2.1 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.1,sigma_beta_1.1)),recursive=F)

# Case X.2 : Different mu_beta small sigma_beta
set.seed(123)
mu_beta_1.2 <- matrix(runif(m*p,2,20),nrow = m,ncol=p,byrow=T)
set.seed(123)
sigma_beta_1.2 <- matrix(runif(m*p,0,1),nrow = m,ncol=p,byrow = T) # Use uniform distribution
#train_1.2 <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.2,sigma_beta_1.2)
test_1.2 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.2,sigma_beta_1.2)),recursive=F)

#train_2.2 <- MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.2,sigma_beta_1.2)
test_2.2 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.2,sigma_beta_1.2)),recursive=F)

# Case X.3 : Similar mu_beta large sigma_beta
set.seed(123)
mu_beta_1.3 <- matrix(runif(m*p,2,5),nrow = m,ncol=p,byrow=T)
set.seed(123)
sigma_beta_1.3 <- matrix(runif(m*p,5,10),nrow = m,ncol=p,byrow = T) # Use uniform distribution

#train_1.3 <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.3,sigma_beta_1.3)
test_1.3 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.3,sigma_beta_1.3)),recursive=F)

#train_2.3 <- MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.3,sigma_beta_1.3)
test_2.3 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.3,sigma_beta_1.3)),recursive=F)

# Case X.4: Different mu_beta large sigma_beta
set.seed(123)
mu_beta_1.4 <- matrix(runif(m*p,2,20),nrow = m,ncol=p,byrow=T)
set.seed(123)
sigma_beta_1.4 <- matrix(runif(m*p,5,10),nrow = m,ncol=p,byrow = T) # Use uniform distribution

#train_1.4 <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.4,sigma_beta_1.4)
test_1.4 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.4,sigma_beta_1.4)),recursive=F)

#train_2.4 <- MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.4,sigma_beta_1.4)
test_2.4 <- unlist(lapply(1:10,function(x) MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.4,sigma_beta_1.4)),recursive=F)

###### Running the Code #####
## Do paralleli ###
library(snow)
library(snowfall)
sfInit( parallel=TRUE, cpus=10 )

sfLibrary(MASS)
sfLibrary(MCMCpack)
sfLibrary(caret)
sfLibrary(readxl)
sfExportAll( )
sfClusterEval(ls())  
case1.1 <- sfLapply(1:50,function(x) get_rmse(get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.1))[[1]],
                                            get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_1.1))[[2]],
                                            lapply(test_1.1, '[[',1)))
p1.1 <- get_boxplot(case1.1)
p1.1 <- p1.1 + ggtitle("Similar mu_beta, small sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF'))

######
case2.1 <- sfLapply(1:50,function(x) get_rmse(get_input_dat(MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.1,sigma_beta_1.1))[[1]],
                                            get_input_dat(MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.1,sigma_beta_1.1))[[2]],
                                            lapply(test_2.1, '[[',1)))
p2.1 <- get_boxplot(case2.1)
p2.1 <- p2.1 + ggtitle("Similar mu_beta, small sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF'))

#####
case1.2 <- sfLapply(1:50,function(x) get_rmse(get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.2,sigma_beta_1.2))[[1]],
                                            get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.2,sigma_beta_1.2))[[2]],
                                            lapply(test_1.2, '[[',1)))
p1.2 <- get_boxplot(case1.2)
p1.2 <- p1.2 + ggtitle("Different mu_beta, small sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF'))

#####
case2.2 <- sfLapply(1:50,function(x) get_rmse(get_input_dat(MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.2,sigma_beta_1.2))[[1]],
                                            get_input_dat(MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.2,sigma_beta_1.2))[[2]],
                                            lapply(test_2.2, '[[',1)))
p2.2 <- get_boxplot(case2.2)
p2.2 <- p2.2 + ggtitle("Different mu_beta, small sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF'))

#####
case1.3 <- sfLapply(1:50,function(x) get_rmse(get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.3,sigma_beta_1.3))[[1]],
                                            get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.3,sigma_beta_1.3))[[2]],
                                            lapply(test_1.3, '[[',1)))
p1.3 <- get_boxplot(case1.3)
p1.3 <- p1.3 + ggtitle("Similar mu_beta large sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF'))

#####
#case2.3 <- sfLapply(1:50,function(x) get_rmse(get_input_dat(MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.3,sigma_beta_1.3))[[1]],
#                                            get_input_dat(MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.3,sigma_beta_1.3))[[2]],
#                                            lapply(test_2.3, '[[',1)))
case2.3 <- readRDS('~/Desktop/master_project/case2.3.rds')
p2.3 <- get_boxplot(case2.3)
p2.3 <- p2.3 + ggtitle("Similar mu_beta large sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF'))

######
#case1.4 <- sfLapply(1:50,function(x) get_rmse(get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.4,sigma_beta_1.4))[[1]],
#                                            get_input_dat(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.4,sigma_beta_1.4))[[2]],
#                                            lapply(test_1.4, '[[',1)))
case1.4 <- readRDS('~/Desktop/master_project/case1.4.rds')

p1.4 <- get_boxplot(case1.4)
p1.4 <- p1.4 + ggtitle("Different mu_beta large sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF'))

#####
#case2.4 <- sfLapply(1:50,function(x) get_rmse(get_input_dat(MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.4,sigma_beta_1.4))[[1]],
#                                            get_input_dat(MultiStudySim(k,nk,p,m,mu_x_2,SIG,pii,mu_beta_1.4,sigma_beta_1.4))[[2]],
#                                            lapply(test_2.4, '[[',1)))
case2.4 <- readRDS('~/Desktop/master_project/case2.4.rds')
p2.4 <- get_boxplot(case2.4)
p2.4 <- p2.4 + ggtitle("Different mu_beta large sigma_beta") + scale_fill_manual(values=c("#F8766D",'#C77CFF' ,'#00BA38', '#619CFF'))



grid.arrange(p1.1,p1.2,p1.3,p1.4,ncol=2,top=textGrob("Different mu_x",gp=gpar(fontsize=20,font=3)))
grid.arrange(p2.1,p2.2,p2.3,p2.4,ncol=2,top=textGrob("Similar mu_x",gp=gpar(fontsize=20,font=3)))



load('~/Desktop/master_project/multisimstudy_test3.0.Rdata')
save.image('~/Desktop/master_project/multisimstudy_test3.0.Rdata')

sfStop() 
showConnections()



