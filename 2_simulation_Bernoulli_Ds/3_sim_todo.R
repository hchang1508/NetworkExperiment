rm(list=ls())

#########################################################
########Generate Simulation Files########################
#########################################################
setwd('/home/hc654/palmer_scratch/Bernoulli_compareDs/')
files=list.files()
temp1_34=files[grep('1_Bernoulli_3_4_FOSO*',files )]
temp1_35=files[grep('1_Bernoulli_3_5_FOSO*',files )]
temp1_36=files[grep('1_Bernoulli_3_6_FOSO*',files )]
temp1_45=files[grep('1_Bernoulli_4_5_FOSO*',files )]
temp1_46=files[grep('1_Bernoulli_4_6_FOSO*',files )]

temp2_34=files[grep('2_Bernoulli_3_4_FOSO*',files )]
temp2_35=files[grep('2_Bernoulli_3_5_FOSO*',files )]
temp2_36=files[grep('2_Bernoulli_3_6_FOSO*',files )]
temp2_45=files[grep('2_Bernoulli_4_5_FOSO*',files )]
temp2_46=files[grep('2_Bernoulli_4_6_FOSO*',files )]

temp3_34=files[grep('3_Bernoulli_3_4_FOSO*',files )]
temp3_35=files[grep('3_Bernoulli_3_5_FOSO*',files )]
temp3_36=files[grep('3_Bernoulli_3_6_FOSO*',files )]
temp3_45=files[grep('3_Bernoulli_4_5_FOSO*',files )]
temp3_46=files[grep('3_Bernoulli_4_6_FOSO*',files )]

temp4_34=files[grep('4_Bernoulli_3_4_FOSO*',files )]
temp4_35=files[grep('4_Bernoulli_3_5_FOSO*',files )]
temp4_36=files[grep('4_Bernoulli_3_6_FOSO*',files )]
temp4_45=files[grep('4_Bernoulli_4_5_FOSO*',files )]
temp4_46=files[grep('4_Bernoulli_4_6_FOSO*',files )]

temp5_34=files[grep('5_Bernoulli_3_4_FOSO*',files )]
temp5_35=files[grep('5_Bernoulli_3_5_FOSO*',files )]
temp5_36=files[grep('5_Bernoulli_3_6_FOSO*',files )]
temp5_45=files[grep('5_Bernoulli_4_5_FOSO*',files )]
temp5_46=files[grep('5_Bernoulli_4_6_FOSO*',files )]

temp6_34=files[grep('6_Bernoulli_3_4_FOSO*',files )]
temp6_35=files[grep('6_Bernoulli_3_5_FOSO*',files )]
temp6_36=files[grep('6_Bernoulli_3_6_FOSO*',files )]
temp6_45=files[grep('6_Bernoulli_4_5_FOSO*',files )]
temp6_46=files[grep('6_Bernoulli_4_6_FOSO*',files )]



#########################################################
########Generate Todo List###############################
#########################################################

source('/home/hc654/NetworkExperiment/2_simulation_Bernoulli_Ds/0_header.R')

case_to_do_1_34=gsub('1_Bernoulli_3_4_FOSO_','',temp1_34)
case_to_do_1_34=gsub('.Rdata','',case_to_do_1_34)
case_to_do_1_34=setdiff(1:10000,as.numeric(case_to_do_1_34))
write.csv(case_to_do_1_34,'case_to_do_b1_34.csv')

case_to_do_1_35=gsub('1_Bernoulli_3_5_FOSO_','',temp1_35)
case_to_do_1_35=gsub('.Rdata','',case_to_do_1_35)
case_to_do_1_35=setdiff(1:10000,as.numeric(case_to_do_1_35))
write.csv(case_to_do_1_35,'case_to_do_b1_35.csv')

case_to_do_1_36=gsub('1_Bernoulli_3_6_FOSO_','',temp1_36)
case_to_do_1_36=gsub('.Rdata','',case_to_do_1_36)
case_to_do_1_36=setdiff(1:10000,as.numeric(case_to_do_1_36))
write.csv(case_to_do_1_36,'case_to_do_b1_36.csv')

case_to_do_1_45=gsub('1_Bernoulli_4_5_FOSO_','',temp1_45)
case_to_do_1_45=gsub('.Rdata','',case_to_do_1_45)
case_to_do_1_45=setdiff(1:10000,as.numeric(case_to_do_1_45))
write.csv(case_to_do_1_45,'case_to_do_b1_45.csv')

case_to_do_1_46=gsub('1_Bernoulli_4_6_FOSO_','',temp1_46)
case_to_do_1_46=gsub('.Rdata','',case_to_do_1_46)
case_to_do_1_46=setdiff(1:10000,as.numeric(case_to_do_1_46))
write.csv(case_to_do_1_46,'case_to_do_b1_46.csv')

case_to_do_2_34=gsub('2_Bernoulli_3_4_FOSO_','',temp2_34)
case_to_do_2_34=gsub('.Rdata','',case_to_do_2_34)
case_to_do_2_34=setdiff(1:10000,as.numeric(case_to_do_2_34))
write.csv(case_to_do_2_34,'case_to_do_b2_34.csv')

case_to_do_2_35=gsub('2_Bernoulli_3_5_FOSO_','',temp2_35)
case_to_do_2_35=gsub('.Rdata','',case_to_do_2_35)
case_to_do_2_35=setdiff(1:10000,as.numeric(case_to_do_2_35))
write.csv(case_to_do_2_35,'case_to_do_b2_35.csv')

case_to_do_2_36=gsub('2_Bernoulli_3_6_FOSO_','',temp2_36)
case_to_do_2_36=gsub('.Rdata','',case_to_do_2_36)
case_to_do_2_36=setdiff(1:10000,as.numeric(case_to_do_2_36))
write.csv(case_to_do_2_36,'case_to_do_b2_36.csv')

case_to_do_2_45=gsub('2_Bernoulli_4_5_FOSO_','',temp2_45)
case_to_do_2_45=gsub('.Rdata','',case_to_do_2_45)
case_to_do_2_45=setdiff(1:10000,as.numeric(case_to_do_2_45))
write.csv(case_to_do_2_45,'case_to_do_b2_45.csv')

case_to_do_2_46=gsub('2_Bernoulli_4_6_FOSO_','',temp2_46)
case_to_do_2_46=gsub('.Rdata','',case_to_do_2_46)
case_to_do_2_46=setdiff(1:10000,as.numeric(case_to_do_2_46))
write.csv(case_to_do_2_46,'case_to_do_b2_46.csv')

case_to_do_3_34=gsub('3_Bernoulli_3_4_FOSO_','',temp3_34)
case_to_do_3_34=gsub('.Rdata','',case_to_do_3_34)
case_to_do_3_34=setdiff(1:10000,as.numeric(case_to_do_3_34))
write.csv(case_to_do_3_34,'case_to_do_b3_34.csv')

case_to_do_3_35=gsub('3_Bernoulli_3_5_FOSO_','',temp3_35)
case_to_do_3_35=gsub('.Rdata','',case_to_do_3_35)
case_to_do_3_35=setdiff(1:10000,as.numeric(case_to_do_3_35))
write.csv(case_to_do_3_35,'case_to_do_b3_35.csv')

case_to_do_3_36=gsub('3_Bernoulli_3_6_FOSO_','',temp3_36)
case_to_do_3_36=gsub('.Rdata','',case_to_do_3_36)
case_to_do_3_36=setdiff(1:10000,as.numeric(case_to_do_3_36))
write.csv(case_to_do_3_36,'case_to_do_b3_36.csv')

case_to_do_3_45=gsub('3_Bernoulli_4_5_FOSO_','',temp3_45)
case_to_do_3_45=gsub('.Rdata','',case_to_do_3_45)
case_to_do_3_45=setdiff(1:10000,as.numeric(case_to_do_3_45))
write.csv(case_to_do_3_45,'case_to_do_b3_45.csv')

case_to_do_3_46=gsub('3_Bernoulli_4_6_FOSO_','',temp3_46)
case_to_do_3_46=gsub('.Rdata','',case_to_do_3_46)
case_to_do_3_46=setdiff(1:10000,as.numeric(case_to_do_3_46))
write.csv(case_to_do_3_46,'case_to_do_b3_46.csv')

case_to_do_4_34=gsub('4_Bernoulli_3_4_FOSO_','',temp4_34)
case_to_do_4_34=gsub('.Rdata','',case_to_do_4_34)
case_to_do_4_34=setdiff(1:10000,as.numeric(case_to_do_4_34))
write.csv(case_to_do_4_34,'case_to_do_b4_34.csv')

case_to_do_4_35=gsub('4_Bernoulli_3_5_FOSO_','',temp4_35)
case_to_do_4_35=gsub('.Rdata','',case_to_do_4_35)
case_to_do_4_35=setdiff(1:10000,as.numeric(case_to_do_4_35))
write.csv(case_to_do_4_35,'case_to_do_b4_35.csv')

case_to_do_4_36=gsub('4_Bernoulli_3_6_FOSO_','',temp4_36)
case_to_do_4_36=gsub('.Rdata','',case_to_do_4_36)
case_to_do_4_36=setdiff(1:10000,as.numeric(case_to_do_4_36))
write.csv(case_to_do_4_36,'case_to_do_b4_36.csv')

case_to_do_4_45=gsub('4_Bernoulli_4_5_FOSO_','',temp4_45)
case_to_do_4_45=gsub('.Rdata','',case_to_do_4_45)
case_to_do_4_45=setdiff(1:10000,as.numeric(case_to_do_4_45))
write.csv(case_to_do_4_45,'case_to_do_b4_45.csv')

case_to_do_4_46=gsub('4_Bernoulli_4_6_FOSO_','',temp4_46)
case_to_do_4_46=gsub('.Rdata','',case_to_do_4_46)
case_to_do_4_46=setdiff(1:10000,as.numeric(case_to_do_4_46))
write.csv(case_to_do_4_46,'case_to_do_b4_46.csv')

case_to_do_5_34=gsub('5_Bernoulli_3_4_FOSO_','',temp5_34)
case_to_do_5_34=gsub('.Rdata','',case_to_do_5_34)
case_to_do_5_34=setdiff(1:10000,as.numeric(case_to_do_5_34))
write.csv(case_to_do_5_34,'case_to_do_b5_34.csv')
case_to_do_5_35=gsub('5_Bernoulli_3_5_FOSO_','',temp5_35)
case_to_do_5_35=gsub('.Rdata','',case_to_do_5_35)
case_to_do_5_35=setdiff(1:10000,as.numeric(case_to_do_5_35))
write.csv(case_to_do_5_35,'case_to_do_b5_35.csv')
case_to_do_5_36=gsub('5_Bernoulli_3_6_FOSO_','',temp5_36)
case_to_do_5_36=gsub('.Rdata','',case_to_do_5_36)
case_to_do_5_36=setdiff(1:10000,as.numeric(case_to_do_5_36))
write.csv(case_to_do_5_36,'case_to_do_b5_36.csv')
case_to_do_5_45=gsub('5_Bernoulli_4_5_FOSO_','',temp5_45)
case_to_do_5_45=gsub('.Rdata','',case_to_do_5_45)
case_to_do_5_45=setdiff(1:10000,as.numeric(case_to_do_5_45))
write.csv(case_to_do_5_45,'case_to_do_b5_45.csv')
case_to_do_5_46=gsub('5_Bernoulli_4_6_FOSO_','',temp5_46)
case_to_do_5_46=gsub('.Rdata','',case_to_do_5_46)
case_to_do_5_46=setdiff(1:10000,as.numeric(case_to_do_5_46))
write.csv(case_to_do_5_46,'case_to_do_b5_46.csv')


case_to_do_6_34=gsub('6_Bernoulli_3_4_FOSO_','',temp6_34)
case_to_do_6_34=gsub('.Rdata','',case_to_do_6_34)
case_to_do_6_34=setdiff(1:10000,as.numeric(case_to_do_6_34))
write.csv(case_to_do_6_34,'case_to_do_b6_34.csv')
case_to_do_6_35=gsub('6_Bernoulli_3_5_FOSO_','',temp6_35)
case_to_do_6_35=gsub('.Rdata','',case_to_do_6_35)
case_to_do_6_35=setdiff(1:10000,as.numeric(case_to_do_6_35))
write.csv(case_to_do_6_35,'case_to_do_b6_35.csv')
case_to_do_6_36=gsub('6_Bernoulli_3_6_FOSO_','',temp6_36)
case_to_do_6_36=gsub('.Rdata','',case_to_do_6_36)
case_to_do_6_36=setdiff(1:10000,as.numeric(case_to_do_6_36))
write.csv(case_to_do_6_36,'case_to_do_b5_36.csv')
case_to_do_6_45=gsub('6_Bernoulli_4_5_FOSO_','',temp6_45)
case_to_do_6_45=gsub('.Rdata','',case_to_do_6_45)
case_to_do_6_45=setdiff(1:10000,as.numeric(case_to_do_6_45))
write.csv(case_to_do_6_45,'case_to_do_b6_45.csv')
case_to_do_6_46=gsub('6_Bernoulli_4_6_FOSO_','',temp6_46)
case_to_do_6_46=gsub('.Rdata','',case_to_do_6_46)
case_to_do_6_46=setdiff(1:10000,as.numeric(case_to_do_6_46))
write.csv(case_to_do_6_46,'case_to_do_b6_46.csv')




#########################################################
########Generate Simulation Files########################
#########################################################
setwd('/home/hc654/palmer_scratch/Bernoulli_compareDs/simulation_output/')
files=list.files()
temp4_34=files[grep('4_3_4_*',files )]
temp4_35=files[grep('4_3_5_*',files )]
temp4_36=files[grep('4_3_6_*',files )]
temp4_45=files[grep('4_4_5_*',files )]
temp4_46=files[grep('4_3_6_*',files )]

#########################################################
########Generate Todo List###############################
#########################################################

source('/home/hc654/NetworkExperiment/2_simulation_s/0_header.R')

case_to_do_4_34=gsub('4_3_4_Sim_','',temp4_34)
case_to_do_4_34=gsub('.Rdata','',case_to_do_4_34)
case_to_do_4_34=setdiff(1:10000,as.numeric(case_to_do_4_34))
write.csv(case_to_do_4_34,'Sim_to_do_b4_34.csv')

case_to_do_4_35=gsub('4_3_5_Sim_','',temp4_35)
case_to_do_4_35=gsub('.Rdata','',case_to_do_4_35)
case_to_do_4_35=setdiff(1:10000,as.numeric(case_to_do_4_35))
write.csv(case_to_do_4_35,'Sim_to_do_b4_35.csv')

case_to_do_4_36=gsub('4_3_6_Sim_','',temp4_36)
case_to_do_4_36=gsub('.Rdata','',case_to_do_4_36)
case_to_do_4_36=setdiff(1:10000,as.numeric(case_to_do_4_36))
write.csv(case_to_do_4_36,'Sim_to_do_b4_36.csv')

case_to_do_4_45=gsub('4_4_5_Sim_','',temp4_45)
case_to_do_4_45=gsub('.Rdata','',case_to_do_4_45)
case_to_do_4_45=setdiff(1:10000,as.numeric(case_to_do_4_45))
write.csv(case_to_do_4_45,'Sim_to_do_b4_45.csv')

case_to_do_4_46=gsub('4_4_6_Sim_','',temp4_46)
case_to_do_4_46=gsub('.Rdata','',case_to_do_4_46)
case_to_do_4_46=setdiff(1:10000,as.numeric(case_to_do_4_46))
write.csv(case_to_do_4_46,'Sim_to_do_b4_46.csv')

#########################################################