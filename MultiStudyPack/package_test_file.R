k <- 10
nk <- rep(100,k)
p <- 3
m <- 4
mu_x <- matrix(runif(k*p,1,20),nrow = k,ncol = p,byrow = T)
SIG <- MCMCpack::riwish(v = p+1, S = diag(1, p,p))
pii <- matrix(runif(m,0,1),nrow = m,ncol = 1,byrow = T)
mu_beta <- matrix(runif(m*p,1,20),nrow = m,ncol = p,byrow = T)
sigma_beta <- matrix(runif(m*p,1,1.1),nrow = m,ncol = p,byrow = T)

lst_ori <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta)
a <- get_input_data(lst_ori)[[1]]
b <- get_input_data(lst_ori)[[2]]
test <- lapply(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta),'[[',1)
l <- 'glmnet'
nk_train <- nk
result_1 <- method_class_all(a,b,test,'glmnet')

result_2 <- method_class_all(a,b,test,'glmnet',nk_train)


result_4 <- method_class_all(a,b,test,'nnet')


