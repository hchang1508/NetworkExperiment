#packages 
library(statnet)

########################################################
#############  set working path  and load data##########
########################################################
#set working path
root = '/home/hc654/NetworkExperiment/'

working_path=paste0(root,'2_simulation_Bernoulli_Ds/')
data_path=paste0(root,'FinalData/net_complete_natvillage.RData')


setwd(working_path)
load(data_path)


########################################################

#subject_size=4509
#data_path_Y14='/home/hc654/Unified/final_analysis/data'
#data_path_Y14=paste0(data_path_Y14,'/Ds_nat_14_y.csv')
#Y14=as.matrix(read.csv(data_path_Y14,header=FALSE))
#Y14=1-Y14
#y0_opt=Y14[2*(1:subject_size)-1]
#y1_opt=Y14[2*(1:subject_size)]

data_path_Y_impute=paste0(root,'FinalData/Y_imputed.csv')
Y_impute=read.csv(data_path_Y_impute)
Y_impute=Y_impute[,2:9]


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

