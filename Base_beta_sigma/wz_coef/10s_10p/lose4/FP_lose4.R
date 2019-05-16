case_0 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/case_0.RDS')
case_0.25 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/case_0.25.RDS')
case_0.5 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/case_0.5.RDS')
case_0.75 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/case_0.75.RDS')

case_1 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/case_1.RDS')
case_1.5 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/case_1.5.RDS')
case_2 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/case_2.RDS')
case_2.5 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/case_2.5.RDS')
case_3 <- readRDS('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/case_3.RDS')

library(ggplot2)
library(grid)
library(gridExtra)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  #datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  #ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  #datac$ci <- datac$se * ciMult
  
  return(datac)
}

get_long_data <- function(output,size){
  tt_elnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',1)))[,c(1:5)], weighting, measurement, stacking:merge, factor_key=TRUE),'elnet')
  tt_nnet <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',2)))[,c(1:5)], weighting, measurement, stacking:merge, factor_key=TRUE),'Nnet')
  
  #tt_boost <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',3)))[,c(1:5)], weighting, measurement, stacking:merge, factor_key=TRUE),'boost')
  tt_boost <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',3)))[,c(1:4,6)], weighting, measurement, stacking:stacking_wz, factor_key=TRUE),'boost')
  
  tt_rf <- cbind(tidyr::gather(data.frame(do.call(rbind,lapply(output,'[[',4)))[,c(1:5)], weighting, measurement, stacking:merge, factor_key=TRUE),'rf')
  colnames(tt_elnet)=colnames(tt_nnet)=colnames(tt_boost)=colnames(tt_rf)
  plotting <- data.frame(rbind(tt_elnet,tt_nnet,tt_boost,tt_rf))
  colnames(plotting)[3] <- 'Learners'
  plotting$measurement <- as.numeric(plotting$measurement)
  plotting$measurement <- log(plotting$measurement)
  plotting <- cbind(plotting,size)
  names(plotting)[4] <- 'Size'
  #re <- ggplot(data = plotting, aes(x=Learners, y=measurement,fill=condition)) + geom_boxplot()+ xlab("Learners") + ylab("log(RMSE)")
  return(plotting)
}


dat_0 <- get_long_data(case_0,0)
dat_0.25 <- get_long_data(case_0.25,0.25)
dat_0.5 <- get_long_data(case_0.5,0.5)
dat_0.75 <- get_long_data(case_0.75,0.75)
dat_1 <- get_long_data(case_1,1)
dat_1.5 <- get_long_data(case_1.5,1.5)
dat_2 <- get_long_data(case_2,2)
dat_2.5 <- get_long_data(case_2.5,2.5)
dat_3 <- get_long_data(case_3,3)

dat_total <- rbind(dat_0,dat_0.25,dat_0.5,dat_0.75,dat_1,dat_1.5,dat_2,dat_2.5,dat_3)
write.csv(dat_total,'~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_10p/lose4/dat_total.csv')
get_learner_plotting <- function(dat_total,title,xlab){
  
  dat_elnet <- dat_total[which(dat_total$Learners %in% 'elnet'),names(dat_total)!='Learners']
  dat_elnet_plot <- summarySE(dat_elnet, measurevar="measurement", groupvars=c("weighting","Size"))
  pd <- position_dodge(0)
  p_elnet <- ggplot(dat_elnet_plot, aes(x=Size, y=measurement, color=weighting)) + 
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.1,position = pd) +
    geom_line(position = pd) + xlab(xlab) + ylab("log(RMSE)") +
    geom_point(position = pd) + ggtitle("Elnet") +
    scale_x_continuous(breaks = dat_elnet_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  
  
  dat_nnet <- dat_total[which(dat_total$Learners %in% 'Nnet'),names(dat_total)!='Learners']
  dat_nnet_plot <- summarySE(dat_nnet, measurevar="measurement", groupvars=c("weighting","Size"))
  p_nnet <- ggplot(dat_nnet_plot, aes(x=Size, y=measurement, color=weighting)) + 
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.1,position = pd) +
    geom_line(position = pd) + xlab(xlab) + ylab("log(RMSE)") +
    geom_point(position = pd) + ggtitle("Nnet") +
    scale_x_continuous(breaks = dat_nnet_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  
  
  dat_boost <- dat_total[which(dat_total$Learners %in% 'boost'),names(dat_total)!='Learners']
  dat_boost_plot <- summarySE(dat_boost, measurevar="measurement", groupvars=c("weighting","Size"))
  p_boost <- ggplot(dat_boost_plot, aes(x=Size, y=measurement, color=weighting)) + 
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.1,position = pd) +
    geom_line(position = pd) + xlab(xlab) + ylab("log(RMSE)") +
    geom_point(position = pd) + ggtitle("Boosting") +
    scale_x_continuous(breaks = dat_boost_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  
  
  dat_rf <- dat_total[which(dat_total$Learners %in% 'rf'),names(dat_total)!='Learners']
  dat_rf_plot <- summarySE(dat_rf, measurevar="measurement", groupvars=c("weighting","Size"))
  p_rf <- ggplot(dat_rf_plot, aes(x=Size, y=measurement, color=weighting)) + 
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.1,position = pd) +
    geom_line(position = pd) + xlab(xlab) + ylab("log(RMSE)") +
    geom_point(position = pd) + ggtitle("RF") +
    scale_x_continuous(breaks = dat_rf_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  
  
  library(grid)
  return(grid.arrange(p_elnet,p_nnet,p_boost,p_rf,ncol=2,top=textGrob(title,gp=gpar(fontsize=20,font=3))))
  
}

get_learner_plotting(dat_total,'Beta Value [-3,-2]&[2,3],10S_10P_lose1','Beta Sigma')

get_plot_ratio <- function(dat_total,learner,lose){
  dat_plot_learner <- dat_total[which(dat_total$Learners %in% learner),names(dat_total)!='Learners']
  dat_plot_ratio <- dat_plot_learner[-which(dat_plot_learner$weighting %in% 'stacking'),]
  base <-  dat_plot_learner$measurement[which(dat_plot_learner$weighting %in% 'stacking')]
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'stacking_norm'),]$measurement <- dat_plot_learner[which(dat_plot_learner$weighting %in% 'stacking_norm'),]$measurement-base
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'simp_avg'),]$measurement <- dat_plot_learner[which(dat_plot_learner$weighting %in% 'simp_avg'),]$measurement-base
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'merge'),]$measurement <- dat_plot_learner[which(dat_plot_learner$weighting %in% 'merge'),]$measurement-base
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'stacking_wz'),]$measurement <- dat_plot_learner[which(dat_plot_learner$weighting %in% 'stacking_wz'),]$measurement-base
  
  dat_plot <- summarySE(dat_plot_ratio, measurevar="measurement", groupvars=c("weighting","Size"))
  pd <- position_dodge(0)
  p_ratio <- ggplot(dat_plot, aes(x=Size, y=measurement, color=weighting)) +
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.3,position = pd) +
    geom_line(position = pd) + xlab('Beta Sigma') + ylab("log(RMSE_Ratio)") +
    geom_point(position = pd) + ggtitle(paste(learner,lose)) +
    scale_x_continuous(breaks = dat_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  return(p_ratio)
}
grid.arrange(get_plot_ratio(dat_total,'elnet'),get_plot_ratio(dat_total,'Nnet'),
             get_plot_ratio(dat_total,'boost'),get_plot_ratio(dat_total,'rf'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [2,3] 10S_10P_lose1',gp=gpar(fontsize=20,font=3)))




