#' Stacking Regression With Zero Weigthing Method
#'
#' @param pred_train A list that contain multi-study datasets (can be generated from the simulation function)
#' @param dat_train A dataframe that combines all multi-study datasets
#' @param nk_train A vector that contains the study-specific sample size
#' @return The coeficient for each study using stacking regression with zero weighting method
#' @seealso [MultiStudySim()]
#' @seealso [get_input_data()]

#' @examples
#' lst_ori <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta)
#' a <- get_input_data(lst_ori)[[1]]
#' b <- get_input_data(lst_ori)[[2]]
#' fit <- lapply(a, function(a) lm(y~.))
#' b1 <- lapply(fit, function(fit)
#'                      predict(fit,b[,-which(names(b) %in% c('y'))]))
#' n <- nk
#' stack_weight_combine(b1,b,n)


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
