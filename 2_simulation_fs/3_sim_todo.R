rm(list=ls())

#########################################################
########Generate Simulation Files########################
#########################################################
setwd('/home/hc654/palmer_scratch/final_analysis_D_fs/')
files=list.files()
temp12=files[grep('Stratified_fs_1_2_FOSO*',files )]
temp13=files[grep('Stratified_fs_1_3_FOSO*',files )]
temp14=files[grep('Stratified_fs_1_4_FOSO*',files )]
temp15=files[grep('Stratified_fs_1_5_FOSO*',files )]
temp16=files[grep('Stratified_fs_1_6_FOSO*',files )]
temp17=files[grep('Stratified_fs_1_7_FOSO*',files )]

#########################################################
########Generate Todo List###############################
#########################################################

source('/home/hc654/NetworkExperiment/2_simulation_fs/0_header.R')

case_to_do12=gsub('Stratified_fs_1_2_FOSO_','',temp12)
case_to_do12=gsub('.Rdata','',case_to_do12)
case_to_do12=setdiff(1:10000,as.numeric(case_to_do12))
write.csv(case_to_do12,'case_to_do_fs12.csv')

case_to_do13=gsub('Stratified_fs_1_3_FOSO_','',temp13)
case_to_do13=gsub('.Rdata','',case_to_do13)
case_to_do13=setdiff(1:10000,as.numeric(case_to_do13))
write.csv(case_to_do13,'case_to_do_fs13.csv')

case_to_do14=gsub('Stratified_fs_1_4_FOSO_','',temp14)
case_to_do14=gsub('.Rdata','',case_to_do14)
case_to_do14=setdiff(1:10000,as.numeric(case_to_do14))
write.csv(case_to_do14,'case_to_do_fs14.csv')

case_to_do15=gsub('Stratified_fs_1_5_FOSO_','',temp15)
case_to_do15=gsub('.Rdata','',case_to_do15)
case_to_do15=setdiff(1:10000,as.numeric(case_to_do15))
write.csv(case_to_do15,'case_to_do_fs15.csv')

case_to_do16=gsub('Stratified_fs_1_6_FOSO_','',temp16)
case_to_do16=gsub('.Rdata','',case_to_do16)
case_to_do16=setdiff(1:10000,as.numeric(case_to_do16))
write.csv(case_to_do16,'case_to_do_fs16.csv')

case_to_do17=gsub('Stratified_fs_1_7_FOSO_','',temp17)
case_to_do17=gsub('.Rdata','',case_to_do17)
case_to_do17=setdiff(1:10000,as.numeric(case_to_do17))
write.csv(case_to_do17,'case_to_do_fs17.csv')


#setwd('/home/hc654/palmer_scratch/final_analysis_Ds_14_AS2/')
#files=list.files()
#temp14=files[grep('*Oct22*',files )]

#temp14_1=temp14[grep('1_14_Sim*',temp14 )]
#case_to_do14_1=gsub('1_14_Sim_','',temp14_1)
#case_to_do14_1=gsub('_Oct22.Rdata','',case_to_do14_1)
#case_to_do14_1=setdiff(1:3000,as.numeric(case_to_do14_1))
#write.csv(case_to_do14_1,'/home/hc654/Unified/final_analysis/Sim_case_to_do_1_14_R.csv')

#temp14_0=temp14[grep('0_14_Sim*',temp14 )]
#case_to_do14_0=gsub('0_14_Sim_','',temp14_0)
#case_to_do14_0=gsub('_Oct22.Rdata','',case_to_do14_0)
#case_to_do14_0=setdiff(1:3000,as.numeric(case_to_do14_0))
#write.csv(case_to_do14_0,'/home/hc654/Unified/final_analysis/Sim_case_to_do_0_14_R.csv')