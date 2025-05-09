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

temp34=files[grep('Stratified_nat_3_4_FOSO*',files )]
temp35=files[grep('Stratified_nat_3_5_FOSO*',files )]
temp36=files[grep('Stratified_nat_3_6_FOSO*',files )]
temp45=files[grep('Stratified_nat_4_5_FOSO*',files )]
temp46=files[grep('Stratified_nat_4_6_FOSO*',files )]

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

case_to_do34=gsub('Stratified_nat_3_4_FOSO_','',temp34)
case_to_do34=gsub('.Rdata','',case_to_do34)
case_to_do34=setdiff(1:10000,as.numeric(case_to_do34))
write.csv(case_to_do34,'case_to_do_nat34.csv')

case_to_do35=gsub('Stratified_nat_3_5_FOSO_','',temp35)
case_to_do35=gsub('.Rdata','',case_to_do35)
case_to_do35=setdiff(1:10000,as.numeric(case_to_do35))
write.csv(case_to_do35,'case_to_do_nat35.csv')

case_to_do36=gsub('Stratified_nat_3_6_FOSO_','',temp36)
case_to_do36=gsub('.Rdata','',case_to_do36)
case_to_do36=setdiff(1:10000,as.numeric(case_to_do36))
write.csv(case_to_do36,'case_to_do_nat36.csv')

case_to_do45=gsub('Stratified_nat_4_5_FOSO_','',temp45)
case_to_do45=gsub('.Rdata','',case_to_do45)
case_to_do45=setdiff(1:10000,as.numeric(case_to_do45))
write.csv(case_to_do45,'case_to_do_nat45.csv')

case_to_do46=gsub('Stratified_nat_4_6_FOSO_','',temp46)
case_to_do46=gsub('.Rdata','',case_to_do46)
case_to_do46=setdiff(1:10000,as.numeric(case_to_do46))
write.csv(case_to_do46,'case_to_do_nat46.csv')

#########################################################
########Generate Simulation Files########################
#########################################################
setwd('/home/hc654/palmer_scratch/final_analysis_Ds/simulation_output/')
files=list.files()
temp12=files[grep('1_2_*',files )]
temp13=files[grep('1_3*',files )]
temp14=files[grep('1_4*',files )]
temp15=files[grep('1_5*',files )]
temp16=files[grep('1_6*',files )]
temp17=files[grep('1_7*',files )]

#########################################################
########Generate Todo List###############################
#########################################################

source('/home/hc654/NetworkExperiment/2_simulation_s/0_header.R')

case_to_do12=gsub('1_2_Sim_','',temp12)
case_to_do12=gsub('.Rdata','',case_to_do12)
case_to_do12=setdiff(1:10000,as.numeric(case_to_do12))
write.csv(case_to_do12,'Sim_case_to_do_12.csv')

case_to_do13=gsub('1_3_Sim_','',temp13)
case_to_do13=gsub('.Rdata','',case_to_do13)
case_to_do13=setdiff(1:10000,as.numeric(case_to_do13))
write.csv(case_to_do13,'Sim_case_to_do_13.csv')

case_to_do14=gsub('1_4_Sim_','',temp14)
case_to_do14=gsub('.Rdata','',case_to_do14)
case_to_do14=setdiff(1:10000,as.numeric(case_to_do14))
write.csv(case_to_do14,'Sim_case_to_do_14.csv')

case_to_do15=gsub('1_5_Sim_','',temp15)
case_to_do15=gsub('.Rdata','',case_to_do15)
case_to_do15=setdiff(1:10000,as.numeric(case_to_do15))
write.csv(case_to_do15,'Sim_case_to_do_15.csv')

case_to_do16=gsub('1_6_Sim_','',temp16)
case_to_do16=gsub('.Rdata','',case_to_do16)
case_to_do16=setdiff(1:10000,as.numeric(case_to_do16))
write.csv(case_to_do16,'Sim_case_to_do_16.csv')

case_to_do17=gsub('1_7_Sim_','',temp17)
case_to_do17=gsub('.Rdata','',case_to_do17)
case_to_do17=setdiff(1:10000,as.numeric(case_to_do17))
write.csv(case_to_do17,'Sim_case_to_do_17.csv')









