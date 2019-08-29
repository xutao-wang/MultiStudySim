#' Simple Average Weigthing Method
#'
#' @param pred_test A list that cotain the predicted results using the model trained from each single study
#' @examples
#' pred_test <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta)
#' @return A vector of predicted outcome that evaluated taking the mean value of predicted outcome from \code{k} studies


simp_avg <- function(pred_test){
  result <- rep(0,length(as.vector(pred_test[[1]])))
  for (i in 1:length(result)){
    result[i] <- mean(sapply(pred_test, function(x) x[i]))
  }
  return(result)
}
