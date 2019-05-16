stack_weight_int <- function(pred_train,dat_train){
  stack_train_x <- as.matrix.data.frame(cbind(1,-1,do.call(cbind,pred_train)))
  nl <- nnls::nnls(stack_train_x, dat_train$y) 
  return (coef(nl))
}

lst = get_input_dat(MultiStudySim_2(k,nk_100,m,mu_x,sigma_x_0.1,pii,beta_vec))[[1]]
dat_train = get_input_dat(MultiStudySim_2(k,nk_100,m,mu_x,sigma_x_0.1,pii,beta_vec))[[2]]
dat_test = lapply(test_0.1,'[[',1)

ttt = as.matrix.data.frame(do.call(cbind,pred_test[[2]])) %*%nl_int[-c(1:2)]
ttt + (nl_int[1]-nl_int[2])
treebag_class_all <- function(lst,dat_train,dat_test,nk_train){
  con <- caret::trainControl(method = 'cv',number=5,allowParallel=T)
  fit <- lapply(lst,function(lst) caret::train(y~X, data=lst, method = "treebag",trControl=con,preProc = c("center", "scale")))
  
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
  
  return(list(stacking=stacking_rmse_mean,stacking_wz=stacking_wz_rmse_mean,
              simp_avg=simp_avg_rmse_mean,merge=merge_rmse_mean,
              nl_coef=nl,nl_coef_wz=nl_wz))  
}