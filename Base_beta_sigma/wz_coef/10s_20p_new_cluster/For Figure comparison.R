library(ggplot2)
library(grid)
library(gridExtra)
## For comparison ##
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

get_plot_ratio <- function(dat_total,learner,lose){
  dat_plot_learner <- dat_total[which(dat_total$Learners %in% learner),names(dat_total)!='Learners']
  new_plot <- plyr::ddply(dat_plot_learner,.(weighting,Size))
  dat_plot_ratio <- new_plot[-which(new_plot$weighting %in% 'stacking'),]
  base <-  new_plot$measurement[which(new_plot$weighting %in% 'stacking')]
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'stacking_norm'),]$measurement <- new_plot[which(new_plot$weighting %in% 'stacking_norm'),]$measurement-base
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'simp_avg'),]$measurement <- new_plot[which(new_plot$weighting %in% 'simp_avg'),]$measurement-base
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'merge'),]$measurement <- new_plot[which(new_plot$weighting %in% 'merge'),]$measurement-base
  dat_plot_ratio[which(dat_plot_ratio$weighting %in% 'stacking_wz'),]$measurement <- new_plot[which(new_plot$weighting %in% 'stacking_wz'),]$measurement-base
  
  dat_plot <- summarySE(dat_plot_ratio, measurevar="measurement", groupvars=c("weighting","Size"))
  
  pd <- position_dodge(0)
  p_ratio <- ggplot(dat_plot, aes(x=Size, y=measurement, color=weighting)) +
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.3,position = pd) +
    geom_line(position = pd) + xlab('Beta Sigma') + ylab("log(RMSE_Ratio)") +
    geom_point(position = pd) + ggtitle(paste(learner,lose)) +
    scale_x_continuous(breaks = dat_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  return(p_ratio)
}


dat_lose0 <- read.csv('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_20p_new_cluster/lose0/dat_total.csv')
dat_lose5 <- read.csv('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_20p_new_cluster/lose5/dat_total.csv')
dat_lose15 <- read.csv('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_20p_new_cluster/lose15/dat_total.csv')
dat_lose19 <- read.csv('~/Desktop/master_project/base_case/Base_beta_sigma/wz_coef/10s_20p_new_cluster/lose19/dat_total.csv')

grid.arrange(get_plot_ratio(dat_lose0,'elnet','lose0'),get_plot_ratio(dat_lose5,'elnet','lose5'),
             get_plot_ratio(dat_lose15,'elnet','lose15'),get_plot_ratio(dat_lose19,'elnet','lose19'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P',gp=gpar(fontsize=20,font=3)))
grid.arrange(get_plot_ratio(dat_lose0,'Nnet','lose0'),get_plot_ratio(dat_lose5,'Nnet','lose5'),
             get_plot_ratio(dat_lose15,'Nnet','lose15'),get_plot_ratio(dat_lose19,'Nnet','lose19'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P',gp=gpar(fontsize=20,font=3)))
grid.arrange(get_plot_ratio(dat_lose0,'boost','lose0'),get_plot_ratio(dat_lose5,'boost','lose5'),
             get_plot_ratio(dat_lose15,'boost','lose15'),get_plot_ratio(dat_lose19,'boost','lose19'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P',gp=gpar(fontsize=20,font=3)))
grid.arrange(get_plot_ratio(dat_lose0,'rf','lose0'),get_plot_ratio(dat_lose5,'rf','lose5'),
             get_plot_ratio(dat_lose15,'rf','lose15'),get_plot_ratio(dat_lose19,'rf','lose19'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P',gp=gpar(fontsize=20,font=3)))

### Rearranging According to 0 coef ###
grid.arrange(get_plot_ratio(dat_lose0,'elnet','lose0'),get_plot_ratio(dat_lose5,'Nnet','lose0'),
             get_plot_ratio(dat_lose15,'boost','lose0'),get_plot_ratio(dat_lose19,'rf','lose0'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P',gp=gpar(fontsize=20,font=3)))
grid.arrange(get_plot_ratio(dat_lose0,'elnet','lose5'),get_plot_ratio(dat_lose5,'Nnet','lose5'),
             get_plot_ratio(dat_lose15,'boost','lose5'),get_plot_ratio(dat_lose19,'rf','lose5'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P',gp=gpar(fontsize=20,font=3)))
grid.arrange(get_plot_ratio(dat_lose0,'elnet','lose15'),get_plot_ratio(dat_lose5,'Nnet','lose15'),
             get_plot_ratio(dat_lose15,'boost','lose15'),get_plot_ratio(dat_lose19,'rf','lose15'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P',gp=gpar(fontsize=20,font=3)))
grid.arrange(get_plot_ratio(dat_lose0,'elnet','lose19'),get_plot_ratio(dat_lose5,'Nnet','lose19'),
             get_plot_ratio(dat_lose15,'boost','lose19'),get_plot_ratio(dat_lose19,'rf','lose19'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P',gp=gpar(fontsize=20,font=3)))


get_learner_plotting <- function(dat_total,xlab,learner){
  
  dat_mod <- dat_total[which(dat_total$Learners %in% learner),names(dat_total)!='Learners']
  dat_plot <- summarySE(dat_mod, measurevar="measurement", groupvars=c("weighting","Size"))
  pd <- position_dodge(0)
  p_plot <- ggplot(dat_plot, aes(x=Size, y=measurement, color=weighting)) + 
    geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=0.1,position = pd) +
    geom_line(position = pd) + xlab(xlab) + ylab("log(RMSE)") +
    geom_point(position = pd) + ggtitle(learner) +
    scale_x_continuous(breaks = dat_plot$Size) + theme(axis.text.x = element_text(angle = 0))
  
  return(p_plot)
  
  
}

grid.arrange(get_learner_plotting(dat_lose0,'Beta Sigma lose0','elnet'),get_learner_plotting(dat_lose5,'Beta Sigma lose5','elnet'),
             get_learner_plotting(dat_lose15,'Beta Sigma lose15','elnet'),get_learner_plotting(dat_lose19,'Beta Sigma lose19','elnet'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P Elnet',gp=gpar(fontsize=20,font=3)))
grid.arrange(get_learner_plotting(dat_lose0,'Beta Sigma lose0','Nnet'),get_learner_plotting(dat_lose5,'Beta Sigma lose5','Nnet'),
             get_learner_plotting(dat_lose15,'Beta Sigma lose15','Nnet'),get_learner_plotting(dat_lose19,'Beta Sigma lose19','Nnet'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P Nnet',gp=gpar(fontsize=20,font=3)))
grid.arrange(get_learner_plotting(dat_lose0,'Beta Sigma lose0','boost'),get_learner_plotting(dat_lose5,'Beta Sigma lose5','boost'),
             get_learner_plotting(dat_lose15,'Beta Sigma lose15','boost'),get_learner_plotting(dat_lose19,'Beta Sigma lose19','boost'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P boost',gp=gpar(fontsize=20,font=3)))
grid.arrange(get_learner_plotting(dat_lose0,'Beta Sigma lose0','rf'),get_learner_plotting(dat_lose5,'Beta Sigma lose5','rf'),
             get_learner_plotting(dat_lose15,'Beta Sigma lose15','rf'),get_learner_plotting(dat_lose19,'Beta Sigma lose19','rf'),
             ncol=2,top=textGrob('Beta Value [-3,-2] & [-6,-5] 10S_20P rf',gp=gpar(fontsize=20,font=3)))





