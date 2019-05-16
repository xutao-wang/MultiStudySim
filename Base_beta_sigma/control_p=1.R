library(MASS)
library(MCMCpack)
library(caret)
library(monmlp)
### Preliminary Test ###
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
k <- 10
nk <- rep(500,10)
p <- 1

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
test_0 <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0)
dat_input <-  MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta_1.1,sigma_beta_0)

get_input_dat <- function(list){
  maindat <- lapply(list,'[[',1)
  dat_train <- do.call(rbind,maindat)
  return(list(maindat_list=maindat,dat_train=dat_train))
}

lst <- get_input_dat(dat_input)[[1]]
dat_train <- get_input_dat(dat_input)[[2]]
dat_test <- lapply(test_0,'[[',1)[[1]]

fit <- lapply(lst,function(lst) lm(y ~ . + 0, data = lst))
fitM <- lm(y ~ . + 0, data = dat_train)

pred_train <- lapply(fit, function(x) predict(x,dat_train))
pred_test <- lapply(fit, function(x) predict(x,dat_test))
stack_train_x <- as.matrix.data.frame(do.call(cbind,pred_train))  

stack_weight_ls <- function(pred_train,dat_train){
  stack_train_x <- as.matrix.data.frame(do.call(cbind,pred_train))  
  #  stack_test_x <- as.matrix.data.frame(do.call(cbind,pred_test)) 
  nl <- lm(dat_train$y~stack_train_x+0) 
  #  pred <- stack_test_x%*%coef(nl)
  #  pred_norm <- stack_test_x%*%(coef(nl)/sum(coef(nl)))
  return (coef(nl))

}
simp_avg <- function(list1,dat_test){
  result <- rep(0,nrow(dat_test))
  for (i in 1:nrow(dat_test)){
    result[i] <- mean(sapply(list1, function(x) x[i]))
  }
  return(result)
}
# Beta Avg
beta_avg <- mean(sapply(fit,function(x) coef(x)));beta_avg

#Beta_Merge
beta_merge <- coef(fitM);beta_merge

#Beta_Stacking
stacking_train_x <- as.matrix.data.frame(do.call(cbind,pred_train))
nl <- lm(dat_train$y~stack_train_x+0) 
beta_stacking <- coef(nl)[1]*coef(fit[[1]]);beta_stacking

#Beta_Stakcing_with_zero
pred_train_wz <- data.frame(do.call(cbind,pred_train))
index <- lapply(colnames(pred_train_wz),function(x) grep(paste0(x,'.'),rownames(pred_train_wz),fixed = T))
for (i in 1:length(index)){
  pred_train_wz[c(index[[i]]),i] = 0
}
nl_wz <- lm(dat_train$y~as.matrix(pred_train_wz)+0)
kkk <- coef(nl_wz) * 0.9
beta_stacking_wz <- kkk %*% sapply(fit,function(x) coef(x));beta_stacking_wz


