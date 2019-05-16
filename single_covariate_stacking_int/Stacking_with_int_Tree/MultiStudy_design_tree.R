library(MASS)
library(MCMCpack)
library(caret)
library(monmlp)

#### Input:
# Parameter of Interest
# k- No. of studies
# n_k - sample size of each study of k: a vector indicating each length
# Assume p=1 - No. of covariates, assume equal firstly
# mu_x - vectors covariate mean
# sigma_x - vectors covariate variance
# link function - x*beta in our case
# m - No. of mixing component 
# pi - m*1 vector of mxing probability
# beta_vec - vector, 3 elements
#### Output:
# List of k
# Each list with y_k,X_k,beta_k,C_k

MultiStudySim_2 <- function(k,nk,m,mu_x,sigma_x,pii,beta_vec){
  # Sampling random number from cluster (2 in total)
  Ck <- sample(1:m,size=k,replace = T, prob = pii)
  
  # Sampling X from mixture distribution
  x_vec_list <- lapply(1:k,function(x) rnorm(n = nk[x], mu_x[Ck[x]], sigma_x) ) 
  x_vec_ind <- lapply(1:k,function(x) ifelse(x_vec_list[[x]]>0,1,0))
  x_list <- lapply(1:length(x_vec_list), function(y) data.frame(X=x_vec_list[[y]],X1=x_vec_ind[[y]]))
  # Convert beta vector into list
  beta_vec_list <- lapply(1:k,function(z) beta_vec)
  
  # Determine whether X is greater than 0
  #Indicator <- lapply(1:length(x_vec_list), function(z) ifelse(x_vec_list[[z]]>0,1,0)) 
  output_Y_list <- lapply(1:k,function(x) 
    beta_vec_list[[x]][1]+beta_vec_list[[x]][2]*x_list[[x]][,1]+beta_vec_list[[x]][3]*x_list[[x]][,1]*x_list[[x]][,2])
  
  final_data <- lapply(1:k,function(x) list(SimulatedOutput=data.frame(X = x_list[[x]][,1],X1 = x_list[[x]][,2] ,y=output_Y_list[[x]],row.names = c(1:length(output_Y_list[[x]]))),
                                            BetaValue=beta_vec,Ck=Ck[[x]]))
  
  names(final_data) <- paste0('Study',c(1:k))
  return(final_data)
}

k <- 10
nk <- rep(500,k)
nk_100 <- rep(100,k)
nk_train <- nk_100
beta_vec <- c(1,2,-3)
mu_x <- c(-2,2)
sigma_x_0.1 <- 0.1
sigma_x_0.5 <- 0.5
sigma_x_1 <- 1
sigma_x_4 <- 4

m <- 2 # Two clusters
pii <- matrix(rep(0.5,m),nrow=m,ncol=1,byrow=T)

#ttt = MultiStudySim_2(k,nk_100,m,mu_x,sigma_x_0.1,pii,beta_vec)
#MultiStudySim_2(k,nk_100,m,mu_x,sigma_x,pii,beta_vec)
set.seed(11)
test_0.1 <- unlist(lapply(1:10,function(x) MultiStudySim_2(k,nk,m,mu_x,sigma_x_0.1,pii,beta_vec)),recursive=F)

set.seed(12)
test_0.5 <- unlist(lapply(1:10,function(x) MultiStudySim_2(k,nk,m,mu_x,sigma_x_0.5,pii,beta_vec)),recursive=F)

set.seed(15)
test_1 <- unlist(lapply(1:10,function(x) MultiStudySim_2(k,nk,m,mu_x,sigma_x_1,pii,beta_vec)),recursive=F)

set.seed(13)
test_4 <- unlist(lapply(1:10,function(x) MultiStudySim_2(k,nk,m,mu_x,sigma_x_4,pii,beta_vec)),recursive=F)

dat_input_0.1 <- lapply(1:50,function(x) MultiStudySim_2(k,nk_100,m,mu_x,sigma_x_0.1,pii,beta_vec))
dat_input_0.5 <- lapply(1:50,function(x) MultiStudySim_2(k,nk_100,m,mu_x,sigma_x_0.5,pii,beta_vec))
dat_input_1 <- lapply(1:50,function(x) MultiStudySim_2(k,nk_100,m,mu_x,sigma_x_1,pii,beta_vec))
dat_input_4 <- lapply(1:50,function(x) MultiStudySim_2(k,nk_100,m,mu_x,sigma_x_4,pii,beta_vec))

##### Simple Size Avg ####
simp_avg <- function(list1,dat_test){
  result <- rep(0,nrow(dat_test))
  for (i in 1:nrow(dat_test)){
    result[i] <- mean(sapply(list1, function(x) x[i]))
  }
  return(result)
}
#pred_simp_avg <- lapply(1:length(dat_test), function(i) 
#  sapply(1:nrow(dat_test[[1]]), function(y) mean(sapply(pred_test[[i]], function(x) x[y]))))

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
  #stack_train_x <- as.matrix.data.frame(cbind(-1,1,do.call(cbind,pred_train)))
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


elnet_class_all <- function(lst,dat_train,dat_test,nk_train){
  # List: generate from multiStudySim
  # dat_train: one large dataframe (unlist List)
  # dat_testL List of another new data
  
  con <- caret::trainControl(method = 'cv',number=5,allowParallel = T)
  fit <- lapply(lst,function(lst) caret::train(y ~ X + X:X1, data=lst, method = "glmnet",trControl=con,preProc = c("center", "scale")))
  #insamp_rmse_list <- lapply(1:length(fit),function(x) 
  #  Metrics::rmse(predict(fit[[x]],lst[[x]][,-which(names(lst[[x]]) %in% c('y'))]),lst[[x]]$y))
  
  fitM <- caret::train(y ~ X + X:X1, data=dat_train, method = "glmnet",trControl=con,preProc = c("center", "scale"))
  
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
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,
              nl_coef=nl,nl_coef_wz=nl_wz))
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
  
  predM <- lapply(dat_test,function(dat_test) monmlp.predict(as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))]),fitM))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  #insamp_rmse_mean <- mean(unlist(insamp_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,
              nl_coef=nl,nl_coef_wz=nl_wz))  
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
  
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,dat_test)) 
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
   #insamp_rmse_mean <- mean(unlist(insamp_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,
              nl_coef=nl,nl_coef_wz=nl_wz))
}

lst_sum <- MultiStudySim_2(k,rep(100,10),m,mu_x,sigma_x_1,pii,beta_vec)
lst = get_input_dat(lst_sum)[[1]]
dat_train = get_input_dat(lst_sum)[[2]]
dat_test = lapply(test_1,'[[',1)

rf_class_all <- function(lst,dat_train,dat_test,nk_train){
  con <- caret::trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "ranger",mtry=1,trControl=con,preProc = c("center", "scale")))
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
  
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,as.matrix(dat_test[,-which(names(dat_test) %in% c('y'))])))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))
  
   #insamp_rmse_mean <- mean(unlist(insamp_rmse_list))
  return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,
              nl_coef=nl,nl_coef_wz=nl_wz))  
}
treebag_class_all <- function(lst,dat_train,dat_test,nk_train){
  con <- caret::trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "treebag",trControl=con,preProc = c("center", "scale")))

  fitM <- caret::train(y~., data=dat_train,method = "treebag",trControl=con,preProc = c("center", "scale"))
  
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
  
  
  predM <- lapply(dat_test,function(dat_test) predict(fitM,dat_test))
  merge_rmse_list <- lapply(1:length(dat_test), function(x) Metrics::rmse(predM[[x]],dat_test[[x]]$y))
  
  stacking_rmse_mean <- mean(unlist(stacking_rmse_list))
  simp_avg_rmse_mean <- mean(unlist(simp_avg_rmse_list))
  merge_rmse_mean <- mean(unlist(merge_rmse_list))

  return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,stacking_int=stacking_int_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,
              nl_coef=nl,nl_coef_wz=nl_wz))  
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
  #rf <- rf_class_all(list,dat_train,dat_test,nk_train)
  treebag <- treebag_class_all(list,dat_train,dat_test,nk_train)
  output <- list(elnet=elnet,nnet=nnet,boost=boost,treebag=treebag)
  
  return(output)
}

#### Some Test ####
#set.seed(999)
#lst = MultiStudySim_2(k,nk_100,m,mu_x,sigma_x_1,pii,beta_vec)
#set.seed(555)
#test_1 <- unlist(lapply(1:10,function(x) MultiStudySim_2(k,nk,m,mu_x,sigma_x_1,pii,beta_vec)),recursive=F)
set.seed(998)
# set sample size to 500 in traning sets
lst_sum_0.1 <- MultiStudySim_2(k,nk,m,mu_x,sigma_x_0.1,pii,beta_vec)

saveRDS(lst_sum_0.1,'~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/lst_sum_0.1.rds')
saveRDS(test_0.1,'~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/test_0.1.rds')

set.seed(990)
lst_sum_0.5 <- MultiStudySim_2(k,nk,m,mu_x,sigma_x_0.5,pii,beta_vec)
saveRDS(lst_sum_0.5,'~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/lst_sum_0.5.rds')
saveRDS(test_0.5,'~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/test_0.5.rds')

set.seed(92)
# set sample size to 500 in traning sets
lst_sum_1 <- MultiStudySim_2(k,nk,m,mu_x,sigma_x_1,pii,beta_vec)


set.seed(9)
lst_sum_4 <- MultiStudySim_2(k,nk,m,mu_x,sigma_x_4,pii,beta_vec)
saveRDS(lst_sum_4,'~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/lst_sum_4.rds')
saveRDS(test_4,'~/Desktop/master_project/base_case/single_covariate_stacking_int/Stacking_with_int_Tree/test_4.rds')

lst_sum <- MultiStudySim_2(k,nk_100,m,mu_x,sigma_x_1,pii,beta_vec)
lst = get_input_dat(lst_sum)[[1]]
dat_train = get_input_dat(lst_sum)[[2]]
dat_test = lapply(test_1,'[[',1)

#plot(dat_train$X,dat_train$y)
#test = elnet_class_all(lst,dat_train,dat_test,nk_train) # work

#testnnet = nnet_class_all(lst,dat_train,dat_test,nk_train) 

#testboost = boost_class_all(lst,dat_train,dat_test,nk_train)

#testrf = rf_class_all(lst,dat_train,dat_test,nk_train) # Does not work for 1 variable??
#test_treebag = treebag_class_all(lst,dat_train,dat_test,nk_train)

## Do parallel ##
library(snow)
library(snowfall)
sfInit(parallel=TRUE, cpus=18)

sfLibrary(MASS)
sfLibrary(MCMCpack)
sfLibrary(caret)
sfLibrary(monmlp)
sfExportAll( )

case_0.1 <- sfLapply(1:length(dat_input_0.1),function(x) get_rmse(get_input_dat(dat_input_0.1[[x]])[[1]],
                                                              get_input_dat(dat_input_0.1[[x]])[[2]],
                                                              lapply(test_0.1, '[[',1),nk_train))
#saveRDS(case_0.1, file = "~/Desktop/master_project/base_case/10s_10p_new_way/stacking_with_int/case_0.1.rds")
saveRDS(case_0.1, file = "~/MultiSim/result/10s_10p_new_way/stacking_with_int/case_0.1.rds")

case_0.5 <- sfLapply(1:length(dat_input_0.5),function(x) get_rmse(get_input_dat(dat_input_0.5[[x]])[[1]],
                                                                      get_input_dat(dat_input_0.5[[x]])[[2]],
                                                                      lapply(test_0.5, '[[',1),nk_train))
#saveRDS(case_0.5, file = "~/Desktop/master_project/base_case/10s_10p_new_way/stacking_with_int/case_0.5.rds")
saveRDS(case_0.5, file = "~/MultiSim/result/10s_10p_new_way/stacking_with_int/case_0.5.rds")

case_1 <- sfLapply(1:length(dat_input_1),function(x) get_rmse(get_input_dat(dat_input_1[[x]])[[1]],
                                                                    get_input_dat(dat_input_1[[x]])[[2]],
                                                                    lapply(test_1, '[[',1),nk_train))
#saveRDS(case_1, file = "~/Desktop/master_project/base_case/10s_10p_new_way/stacking_with_int//case_1.rds")
saveRDS(case_1, file = "~/MultiSim/result/10s_10p_new_way/stacking_with_int/case_1.rds")

case_4 <- sfLapply(1:length(dat_input_4),function(x) get_rmse(get_input_dat(dat_input_4[[x]])[[1]],
                                                                    get_input_dat(dat_input_4[[x]])[[2]],
                                                                    lapply(test_4, '[[',1),nk_train))
#saveRDS(case_4, file = "~/Desktop/master_project/base_case/10s_10p_new_way/stacking_with_int/case_4.rds")
saveRDS(case_4, file = "~/MultiSim/result/10s_10p_new_way/stacking_with_int/case_4.rds")
# Output should be a type of list


