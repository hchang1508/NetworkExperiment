#packages 
library(statnet)

########################################################
#############  set working path  and load data##########
########################################################
#set working path
root = '/home/hc654/NetworkExperiment//'

working_path=paste0(root,'2_simulation_fs/')
data_path=paste0(root,'FinalData/net_complete.RData')


setwd(working_path)
load(data_path)


########################################################

#subject_size=4509
data_path_Y14='/home/hc654/Unified/final_analysis/data'
#data_path_Y14=paste0(data_path_Y14,'/Ds_nat_14_y.csv')
#Y14=as.matrix(read.csv(data_path_Y14,header=FALSE))
#Y14=1-Y14
#y0_opt=Y14[2*(1:subject_size)-1]
#y1_opt=Y14[2*(1:subject_size)]


#data_path_Y2='/home/hc654/Unified/final_analysis/data/'
#data_path_Y2=paste0(data_path_Y2,'/Y_sim2.csv')
#Y_sim2=read.csv(data_path_Y2)
#Y_sim2=Y_sim2[,2:14]

#Import functions
source('./functions/functions_compute_FOSO_inner.R') #functions calculating FOSO probabilities
#source('./functions/functions2.R')
source('./functions/functions_compute_FOSO.R') #functions calculating FOSO probabilitiess
source('./functions/functions_estimators.R') #various estimators
source('./functions/functions_maketable.R')
source('./functions/functions_simulations.R') #scripts for simulations
source('./functions/functions_noharm.R') #no harm estimators
source('./functions/functions_var_estimation.R') #variance bound estimation
source('./functions/functions_opt_logit2.R') #variance bound estimation

