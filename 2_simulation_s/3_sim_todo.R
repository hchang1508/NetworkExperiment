rm(list=ls())

#########################################################
########Generate Simulation Files########################
#########################################################
setwd('/home/hc654/palmer_scratch/final_analysis_Ds/')
files=list.files()
temp12=files[grep('Stratified_nat_1_2_FOSO*',files )]
temp13=files[grep('Stratified_nat_1_3_FOSO*',files )]
temp14=files[grep('Stratified_nat_1_4_FOSO*',files )]
temp15=files[grep('Stratified_nat_1_5_FOSO*',files )]
temp16=files[grep('Stratified_nat_1_6_FOSO*',files )]
temp17=files[grep('Stratified_nat_1_7_FOSO*',files )]

#########################################################
########Generate Todo List###############################
#########################################################

source('/home/hc654/NetworkExperiment/2_simulation_s/0_header.R')

case_to_do12=gsub('Stratified_nat_1_2_FOSO_','',temp12)
case_to_do12=gsub('.Rdata','',case_to_do12)
case_to_do12=setdiff(1:10000,as.numeric(case_to_do12))
write.csv(case_to_do12,'case_to_do_nat12.csv')

case_to_do13=gsub('Stratified_nat_1_3_FOSO_','',temp13)
case_to_do13=gsub('.Rdata','',case_to_do13)
case_to_do13=setdiff(1:10000,as.numeric(case_to_do13))
write.csv(case_to_do13,'case_to_do_nat13.csv')

case_to_do14=gsub('Stratified_nat_1_4_FOSO_','',temp14)
case_to_do14=gsub('.Rdata','',case_to_do14)
case_to_do14=setdiff(1:10000,as.numeric(case_to_do14))
write.csv(case_to_do14,'case_to_do_nat14.csv')

case_to_do15=gsub('Stratified_nat_1_5_FOSO_','',temp15)
case_to_do15=gsub('.Rdata','',case_to_do15)
case_to_do15=setdiff(1:10000,as.numeric(case_to_do15))
write.csv(case_to_do15,'case_to_do_nat15.csv')

case_to_do16=gsub('Stratified_nat_1_6_FOSO_','',temp16)
case_to_do16=gsub('.Rdata','',case_to_do16)
case_to_do16=setdiff(1:10000,as.numeric(case_to_do16))
write.csv(case_to_do16,'case_to_do_nat16.csv')

case_to_do17=gsub('Stratified_nat_1_7_FOSO_','',temp17)
case_to_do17=gsub('.Rdata','',case_to_do17)
case_to_do17=setdiff(1:10000,as.numeric(case_to_do17))
write.csv(case_to_do17,'case_to_do_nat17.csv')
