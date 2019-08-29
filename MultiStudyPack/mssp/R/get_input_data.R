#' Modification of the simulated results
#'
#' @param lst The multi-study results can be generated from the simulation function MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta)
#' @seealso [MultiStudySim()]
#' @examples
#' lst <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta)
#' get_input_data(lst)
#' @return
#' Two lists, first list is the multi-study training dataset, second list is the whole training dataset that combines all multi-study datasets
#'

get_input_data <- function(lst){
  maindat <- lapply(lst,'[[',1)
  dat_train <- do.call(rbind,maindat)
  return(list(maindat_list=maindat,dat_train=dat_train))
}
