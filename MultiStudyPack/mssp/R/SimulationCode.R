#' Simulation function for Multi-Study Dataset.
#'
#' @param k A number indicating the number of Studies.
#' @param nk A vector indicating the number of observations within each study
#' @param p A number indicating the number of covariates within each study. Assuming equal covariates for each study right now
#' @param m A number indicating the number of mix component
#' @param mu_x A k by p matrix indicaitng the mean of i_th covairates for i=1,...,p, and in j_th study for j=1,...,k
#' @param SIG A p by p matrix indicating the variance-covariance matrix of mu_x
#' @param pii A m by 1 vector indicating the mixed probability of choosing a specific cluster of beta
#' @param mu_beta A m by p matrix indicating the mean of coefficient of i_th covariates for i=1,...,p, and q_th cluster for q=1,...,m
#' @param sigma_beta A m by p matrix indicating the variance of coefficient of i_th covariates for i=1,...,p, and q_th cluster for q=1,...,m
#' @return Three lists of multi-study datasets. Three sublists within each \code{k} studies.
#' First: Simulated Data
#' Second: beta Coefficients
#' Third: integer indicating the selectd component
#' @examples
#' k <- 10
#' nk <- rep(100,k)
#' p <- 3
#' m <- 4
#' mu_x <- matrix(runif(k*p,1,20),nrow = k,ncol = p,byrow = T)
#' SIG <- MCMCpack::riwish(v = p+1, S = diag(1, p,p))
#' pii <- matrix(runif(m,0,1),nrow = m,ncol = 1,byrow = T)
#' mu_beta <- matrix(runif(m*p,1,20),nrow = m,ncol = p,byrow = T)
#' sigma_beta <- matrix(runif(m*p,1,1.1),nrow = m,ncol = p,byrow = T)
#' MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta)
#' @author
#' Xutao Wang, \email{xutaow@@bu.edu}

MultiStudySim <- function(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta){
  x_vec_list <- lapply(1:k,function(x) MASS::mvrnorm(n = nk[x], mu_x[x,], SIG))
  Ck <- sample(1:m,k,replace = T,prob = pii)
  beta_vec_list <- lapply(1:k,function(z) apply(matrix(c(1:p),nrow = p,ncol = 1),1,function(x)
    rnorm(1,mu_beta[Ck[z],x],sigma_beta[Ck[z],x])))
  output_Y_list <- lapply(1:k,function(x)
    x_vec_list[[x]]%*%beta_vec_list[[x]] + rnorm(nrow(x_vec_list[[x]]),0,1))

  final_data <- lapply(1:k,function(x) list(SimulatedOutput = data.frame(x_vec_list[[x]],y = output_Y_list[[x]],row.names = c(1:nrow(output_Y_list[[x]]))),
                                            BetaValue=beta_vec_list[[x]],Ck = Ck[x]))

  names(final_data) <- paste0('Study',c(1:k))
  return(final_data)
}



