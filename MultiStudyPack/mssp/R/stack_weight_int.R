#' Stacking Regression With Intercept Weigthing Method
#'
#' @param pred_train A list that contain multi-study datasets (can be generated from the simulation function)
#' @param dat_train A dataframe that combines all multi-study datasets
#' @return The coeficient for each study using stacking regression with intercept weighting method
#' @seealso [MultiStudySim()]
#' @seealso [get_input_data()]

#' @examples
#' lst_ori <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta)
#' a <- get_input_data(lst_ori)[[1]]
#' b <- get_input_data(lst_ori)[[2]]
#' fit <- lapply(a, function(a) lm(y~.))
#' b1 <- lapply(fit, function(fit)
#'                      predict(fit,b[,-which(names(b) %in% c('y'))]))
#' stack_weight_int(b1,b)

stack_weight_int <- function(pred_train,dat_train){
  stack_train_x <- as.matrix.data.frame(cbind(1,-1,do.call(cbind,pred_train)))
  nl <- nnls::nnls(stack_train_x, dat_train$y)
  return (coef(nl))
}
