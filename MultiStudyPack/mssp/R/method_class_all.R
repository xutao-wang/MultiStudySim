#' Model fitting Method
#'
#' @param lst A list that contain trianing multi-study datasets (can be generated from the simulation function)
#' @param dat_train A dataframe that combines all multi-study datasets
#' @param dat_test A list that contains testing multi-study datasets (can be generated from the simulation function)
#' @param learner A string that indicates the model fitting method (build under caret):"glmnet" (for Elastic Net)
#' "gbm" (for gradient boosting), "ranger" (for random forest), "nnet" (for neural network)
#' @param nk_train A vector that only needs when using stacking with zero weighting method
#' @return A list that contain the predicted outcome for each multi-study test dataset across differnt weighting methods using selected model fitting method
#' @seealso [MultiStudySim()]
#' @seealso [get_input_data()]

#' @examples
#' lst_ori <- MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta)
#' a <- get_input_data(lst_ori)[[1]]
#' b <- get_input_data(lst_ori)[[2]]
#' test <- lapply(MultiStudySim(k,nk,p,m,mu_x,SIG,pii,mu_beta,sigma_beta),'[[',1)
#' l <- 'glmnet'
#' nk_train <- nk
#' method_class_all(a,b,test,l)
#' method_class_all(a,b,test,l,nk)
#'
method_class_all <- function(lst,dat_train,dat_test,learner,nk_train){

  con <- caret::trainControl(method = 'cv',number=5,allowParallel = T)
  if (learner == 'glmnet') {

    fit <- lapply(lst,function(lst) caret::train(y ~ ., data=lst, method = "glmnet",trControl=con,preProc = c("center", "scale")))
    pred_train <- lapply(fit, function(fit) predict(fit,dat_train[,-which(names(dat_train) %in% c('y'))]))
    pred_test <- lapply(1:length(dat_test),function(x) predict(fit,dat_test[[x]]))
    pred_avg <- lapply(1:length(pred_test), function(x) simp_avg(pred_test[[x]]))

    nl <- stack_weight_combine(pred_train,dat_train)
    pred_stack <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)

    nl_int <- stack_weight_int(pred_train,dat_train)
    pred_stack_int <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))

    nl_elnet_int <- elnet_weight_int(pred_train, dat_train)
    pred_elent_int <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_elnet_int[-1]+nl_elnet_int[1])

    if(missing(nk_train)){

      predicted <- lapply(1:length(pred_avg), function(x)
        data.frame(cbind(pred_avg[[x]],pred_stack[[x]],pred_stack_int[[x]],pred_elent_int[[x]])))
      column_name <- c('Average','Stacking','Stacking_int','Elnet_int')
      predicted <- lapply(predicted, setNames, column_name)

      return(predicted)

    } else {

      nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
      pred_stack_wz <- lapply(1:length(pred_test), function(x)
        as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)

      predicted <- lapply(1:length(pred_avg), function(x)
        data.frame(cbind(pred_avg[[x]],pred_stack[[x]],pred_stack_int[[x]],pred_stack_wz[[x]],pred_elent_int[[x]])))
      column_name <- c('Average','Stacking','Stacking_int','Stacking_wz','Elnet_int')
      predicted <- lapply(predicted, setNames, column_name)
      return(predicted)

    }


  } else if (learner == 'gbm'){

    fit <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "gbm",trControl=con,preProc = c("center", "scale")))
    pred_train <- lapply(fit, function(fit) predict(fit,dat_train))
    pred_test <- lapply(1:length(dat_test),function(x) predict(fit,dat_test[[x]]  ))

    pred_avg <- lapply(1:length(pred_test), function(x) simp_avg(pred_test[[x]]))

    nl <- stack_weight_combine(pred_train,dat_train)
    pred_stack <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)

    nl_int <- stack_weight_int(pred_train,dat_train)
    pred_stack_int <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))

    nl_elnet_int <- elnet_weight_int(pred_train, dat_train)
    pred_elent_int <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_elnet_int[-1]+nl_elnet_int[1])

    if(missing(nk_train)){

      predicted <- lapply(1:length(pred_avg), function(x)
        data.frame(cbind(pred_avg[[x]],pred_stack[[x]],pred_stack_int[[x]],pred_elent_int[[x]])))
      column_name <- c('Average','Stacking','Stacking_int','Elnet_int')
      predicted <- lapply(predicted, setNames, column_name)

      return(predicted)

    }else {

      nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
      pred_stack_wz <- lapply(1:length(pred_test), function(x)
        as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)

      predicted <- lapply(1:length(pred_avg), function(x)
        data.frame(cbind(pred_avg[[x]],pred_stack[[x]],pred_stack_int[[x]],pred_stack_wz[[x]],pred_elent_int[[x]])))
      column_name <- c('Average','Stacking','Stacking_int','Stacking_wz','Elnet_int')

      predicted <- lapply(predicted, setNames, column_name)
      return(predicted)

    }
  } else if (learner == 'ranger'){

    fit <- lapply(lst,function(lst) caret::train(y~., data=lst, method = "ranger",trControl=con,preProc = c("center", "scale")))
    pred_train <- lapply(fit, function(fit) predict(fit,dat_train))
    pred_test <- lapply(1:length(dat_test),function(x) predict(fit,dat_test[[x]]  ))

    pred_avg <- lapply(1:length(pred_test), function(x) simp_avg(pred_test[[x]]))

    nl <- stack_weight_combine(pred_train,dat_train)
    pred_stack <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)

    nl_int <- stack_weight_int(pred_train,dat_train)
    pred_stack_int <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))

    nl_elnet_int <- elnet_weight_int(pred_train, dat_train)
    pred_elent_int <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_elnet_int[-1]+nl_elnet_int[1])

    if(missing(nk_train)){

      predicted <- lapply(1:length(pred_avg), function(x)
        data.frame(cbind(pred_avg[[x]],pred_stack[[x]],pred_stack_int[[x]],pred_elent_int[[x]])))
      column_name <- c('Average','Stacking','Stacking_int','Elnet_int')
      predicted <- lapply(predicted, setNames, column_name)

      return(predicted)

    }else {

      nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
      pred_stack_wz <- lapply(1:length(pred_test), function(x)
        as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)

      predicted <- lapply(1:length(pred_avg), function(x)
        data.frame(cbind(pred_avg[[x]],pred_stack[[x]],pred_stack_int[[x]],pred_stack_wz[[x]],pred_elent_int[[x]])))
      column_name <- c('Average','Stacking','Stacking_int','Stacking_wz','Elnet_int')

      predicted <- lapply(predicted, setNames, column_name)
      return(predicted)

    }

  } else if (learner == 'nnet'){

    fit <- lapply(lst, function(lst) monmlp::monmlp.fit(x = as.matrix(lst[,names(lst)!='y']), y = as.matrix(lst$y), hidden1 = 2, hidden2 = 2,
                                                        bag = TRUE, iter.max = 500, iter.stopped = 10))
    pred_train <- lapply(fit, function(x) monmlp::monmlp.predict(as.matrix(dat_train[,-which(names(dat_train) %in% c('y'))]),x))
    pred_test <- lapply(1:length(dat_test),function(x)
      lapply(1:length(fit), function(y)  monmlp::monmlp.predict(as.matrix(dat_test[[x]][,names(dat_test[[x]])!='y']),fit[[y]])) )

    pred_avg <- lapply(1:length(pred_test), function(x) simp_avg(pred_test[[x]]))

    nl <- stack_weight_combine(pred_train,dat_train)
    pred_stack <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl)

    nl_int <- stack_weight_int(pred_train,dat_train)
    pred_stack_int <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_int[-c(1:2)]+(nl_int[1]-nl_int[2]))

    nl_elnet_int <- elnet_weight_int(pred_train, dat_train)
    pred_elent_int <- lapply(1:length(pred_test), function(x)
      as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_elnet_int[-1]+nl_elnet_int[1])

    if (missing(nk_train)){

      predicted <- lapply(1:length(pred_avg), function(x)
        data.frame(cbind(pred_avg[[x]],pred_stack[[x]],pred_stack_int[[x]],pred_elent_int[[x]])))
      column_name <- c('Average','Stacking','Stacking_int','Elnet_int')
      predicted <- lapply(predicted, setNames, column_name)

      return(predicted)

    } else {

      nl_wz <- stack_weight_wz(pred_train,dat_train,nk_train)
      pred_stack_wz <- lapply(1:length(pred_test), function(x)
        as.matrix.data.frame(do.call(cbind,pred_test[[x]])) %*%nl_wz)

      predicted <- lapply(1:length(pred_avg), function(x)
        data.frame(cbind(pred_avg[[x]],pred_stack[[x]],pred_stack_int[[x]],pred_stack_wz[[x]],pred_elent_int[[x]])))
      column_name <- c('Average','Stacking','Stacking_int','Stacking_wz','Elnet_int')

      predicted <- lapply(predicted, setNames, column_name)
      return(predicted)

    }


  }


}
