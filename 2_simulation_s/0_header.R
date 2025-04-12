#packages 
library(statnet)

########################################################
#############  set working path  and load data##########
########################################################
#set working path
root = '/home/hc654/NetworkExperiment/'

working_path=paste0(root,'2_simulation_s/')
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
source('./functions/functions.R')
#source('./functions/functions2.R')
source('./functions/functions3.R')
#source('./functions/functions4.R')
#source('./functions/functions5.R')
#source('./functions/functions6.R')
source('./functions/functions7.R')
#source('./functions/functions8.R')
