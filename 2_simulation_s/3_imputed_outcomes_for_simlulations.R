#data imputation for simulations:
rm(list=ls())
library(dplyr)
library(network)
data_save_path='/Users/haoge/Desktop/net_complete_natvillage.RData'
load(data_save_path)

outcome_data_path= '/Users/haoge/Desktop/0422analysis.csv'
outcome_data=read.csv(outcome_data_path)

#outcome data
outcome_data=outcome_data[which(outcome_data$id %in% id_list),c('takeup_survey','id','takeup_survey','male','age','agpop','ricearea_2010','literacy','risk_averse','disaster_prob','delay','intensive','network_obs','network_rate_preintensive','network_rate_presimple')]

#impute number of friends that in first rounds
outcome_data[,'number_friends_first_round_simple']=round(outcome_data$network_obs * outcome_data$network_rate_presimple,0)
outcome_data[,'number_friends_first_round_intensive']=round(outcome_data$network_obs * outcome_data$network_rate_preintensive,0)

table(outcome_data$number_friends_first_round_simple)
table(outcome_data$number_friends_first_round_intensive)

#####################################################################
##########Estimate logit models for each exposure mappings###########
#####################################################################

#first round simple
Y1=outcome_data[which(outcome_data$delay==0 & outcome_data$intensive==0),]
y1_reg=glm(takeup_survey~male+age+agpop+ricearea_2010+literacy+risk_averse+disaster_prob,data=Y1,family = "binomial")
summary(y1_reg)
mean(Y1$takeup_survey)

#first round intensive
Y2=outcome_data[which(outcome_data$delay==0 & outcome_data$intensive==1),]
y2_reg=glm(takeup_survey~male+age+agpop+ricearea_2010+literacy+risk_averse+disaster_prob,data=Y2,family = "binomial")
summary(y2_reg)
mean(Y2$takeup_survey)

#second round simple, with no friends in first round
Y3=outcome_data[which(outcome_data$delay==1 & outcome_data$intensive==0 & outcome_data$number_friends_first_round_intensive==0 & outcome_data$number_friends_first_round_simple==0),]
y3_reg=glm(takeup_survey~male+age+agpop+ricearea_2010+literacy+risk_averse+disaster_prob,data=Y3,family = "binomial")
summary(y3_reg)
mean(Y3$takeup_survey)

#second round simple, with friends in first round simple
Y4=outcome_data[which(outcome_data$delay==1 & outcome_data$intensive==0 & outcome_data$number_friends_first_round_intensive==0 & outcome_data$number_friends_first_round_simple>0),]
y4_reg=glm(takeup_survey~male+age+agpop+ricearea_2010+literacy+risk_averse+disaster_prob,data=Y4,family = "binomial")
summary(y4_reg)
mean(Y4$takeup_survey)

#second round simple, with one friends in first round intensive
Y5=outcome_data[which(outcome_data$delay==1 & outcome_data$intensive==0 & outcome_data$number_friends_first_round_intensive==1 ),]
y5_reg=glm(takeup_survey~male+age+agpop+ricearea_2010+literacy+risk_averse+disaster_prob,data=Y5,family = "binomial")
summary(y5_reg)
mean(Y5$takeup_survey)

#second round simple, with two friends in first round intensive
Y6=outcome_data[which(outcome_data$delay==1 & outcome_data$intensive==0 & outcome_data$number_friends_first_round_intensive==2 ),]
y6_reg=glm(takeup_survey~male+age+agpop+ricearea_2010+literacy+risk_averse+disaster_prob,data=Y6,family = "binomial")
summary(y6_reg)
mean(Y6$takeup_survey)

#second round simple, with two friends in first round intensive
Y7=outcome_data[which(outcome_data$delay==1 & outcome_data$intensive==0 & outcome_data$number_friends_first_round_intensive>2 ),]
y7_reg=glm(takeup_survey~male+age+agpop+ricearea_2010+literacy+risk_averse+disaster_prob,data=Y7,family = "binomial")
summary(y7_reg)
mean(Y7$takeup_survey)



#####################################################################
##########Impute potential outcomes##################################
#####################################################################
#imputed outcomes
Y_imputes=matrix(0,length(id_list),7)

#covariate data
covar=all_info2[which(all_info2$id %in% id_list),c('id','male','age','agpop','ricearea_2010','literacy','risk_averse','disaster_prob')]
#generate latent index
U=runif(length(id_list))

#1st Exposure Mapping: First Round, Simple Info (FR-S)
Y_imputes[,1]= as.numeric( U<= predict(y1_reg,covar,type="response") )
mean(Y_imputes[,1]) 

#2nd Exposure Mapping: First Round, Intensive Info (FR-I)
Y_imputes[,2]= as.numeric( U<= predict(y2_reg,covar,type="response") )  
mean(Y_imputes[,2]) 

#3rd Exposure Mapping: Second Round Simple, No friend in the First Round
Y_imputes[,3]= as.numeric( U<= predict(y3_reg,covar,type="response") )  
mean(Y_imputes[,3])

#4th Exposure Mapping: Second Round Simple, has a friend in FR-S but not FR-I
Y_imputes[,4]= as.numeric( U<= predict(y4_reg,covar,type="response") )  
mean(Y_imputes[,4])

#5th Exposure Maping: Second Round Simple, has only one friend in FR-I
Y_imputes[,5]= as.numeric( U<= predict(y5_reg,covar,type="response") )  
mean(Y_imputes[,5])

#6th Exposure Maping: Second Round Simple, has two friends in FR-I
Y_imputes[,6]= as.numeric( U<= predict(y6_reg,covar,type="response") )  
mean(Y_imputes[,6])

#7th Exposure Maping: Second Round Simple, has more than two friends in FR-I
Y_imputes[,7]= as.numeric( U<= predict(y7_reg,covar,type="response") )  
mean(Y_imputes[,7])

#8th Exposure Mapping: Second Round Intensive, No friend in the First Round

#9th Exposure Mapping: Second Round Intensive, has a friend in FR-S but not FR-I


#10th Exposure Mapping: Second Round Intensive, has only one friend in FR-I


#11th Exposure Maping: Second Round Intensive, has two friends in FR-I


#12th Exposure Maping: Second Round Intensive, has more than two friends in FR-I

Y_imputes=cbind(covar[,c('id')],Y_imputes)

data_save_path='/Users/haoge/Desktop'
data_save_path=paste0(data_save_path,'/Y_imputed.csv')
write.csv(Y_imputes,file=data_save_path)


