rm(list=ls())

#########################################################
########Generate Simulation Files########################
#########################################################
setwd('/home/hc654/palmer_scratch/final_analysis_Bernoulli/')
files=list.files()
temp12=files[grep('Bernoulli_nat_1_2_FOSO*',files )]
temp13=files[grep('Bernoulli_nat_1_3_FOSO*',files )]
temp14=files[grep('Bernoulli_nat_1_4_FOSO*',files )]
temp15=files[grep('Bernoulli_nat_1_5_FOSO*',files )]
temp16=files[grep('Bernoulli_nat_1_6_FOSO*',files )]
temp17=files[grep('Bernoulli_nat_1_7_FOSO*',files )]
temp45=files[grep('Bernoulli_nat_4_5_FOSO*',files )]

#########################################################
########Generate Todo List###############################
#########################################################

source('/home/hc654/NetworkExperiment/2_simulation_Bernoulli/0_header.R')

case_to_do12=gsub('Bernoulli_nat_1_2_FOSO_','',temp12)
case_to_do12=gsub('.Rdata','',case_to_do12)
case_to_do12=setdiff(1:10000,as.numeric(case_to_do12))
write.csv(case_to_do12,'case_to_do_b12.csv')

case_to_do13=gsub('Bernoulli_nat_1_3_FOSO_','',temp13)
case_to_do13=gsub('.Rdata','',case_to_do13)
case_to_do13=setdiff(1:10000,as.numeric(case_to_do13))
write.csv(case_to_do13,'case_to_do_b13.csv')

case_to_do14=gsub('Bernoulli_nat_1_4_FOSO_','',temp14)
case_to_do14=gsub('.Rdata','',case_to_do14)
case_to_do14=setdiff(1:10000,as.numeric(case_to_do14))
write.csv(case_to_do14,'case_to_do_b14.csv')

case_to_do15=gsub('Bernoulli_nat_1_5_FOSO_','',temp15)
case_to_do15=gsub('.Rdata','',case_to_do15)
case_to_do15=setdiff(1:10000,as.numeric(case_to_do15))
write.csv(case_to_do15,'case_to_do_b15.csv')

case_to_do16=gsub('Bernoulli_nat_1_6_FOSO_','',temp16)
case_to_do16=gsub('.Rdata','',case_to_do16)
case_to_do16=setdiff(1:10000,as.numeric(case_to_do16))
write.csv(case_to_do16,'case_to_do_b16.csv')

case_to_do17=gsub('Bernoulli_nat_1_7_FOSO_','',temp17)
case_to_do17=gsub('.Rdata','',case_to_do17)
case_to_do17=setdiff(1:10000,as.numeric(case_to_do17))
write.csv(case_to_do17,'case_to_do_b17.csv')

case_to_do45=gsub('Bernoulli_nat_4_5_FOSO_','',temp45)
case_to_do45=gsub('.Rdata','',case_to_do17)
case_to_do45=setdiff(1:10000,as.numeric(case_to_do45))
write.csv(case_to_do45,'case_to_do_b45.csv')

#########################################################
########Generate Simulation Files########################
#########################################################
setwd('/home/hc654/palmer_scratch/final_analysis_Bernoulli/simulation_output/')
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

source('/home/hc654/NetworkExperiment/2_simulation_Bernoulli/0_header.R')

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









