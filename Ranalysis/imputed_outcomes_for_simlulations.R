#data imputation for simulations:
rm(list=ls())
library(dplyr)
data_save_path='/Users/haoge/Library/CloudStorage/GoogleDrive-hgchang1508@gmail.com/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/FinalData/'
data_save_path=paste0(data_save_path,'net_complete_natvillage.RData')
load(data_save_path)

subnet=get.inducedSubgraph(net,v=which(net %v%'vertex.names'%in% id_list))
o=as.sociomatrix(subnet)
o=diag(apply(o,1,sum))-o
#imputed outcomes
Y_sim2=matrix(0,length(id_list),12)


#base line take-up rate is 35 percent (P90; also can be computed from the dataset by 
#using people from first group simple session)
#the paper used the following variables for adjustment (Table 2 column 1, coefficients are retrieved from the stata program)
# male age agpop ricearea_2010 literacy risk_averse disaster_prob
# the paper estimated a linear probability model, we convert the coefficients by using the rule
# beta_{LinearRegression}= 4beta_{LogisticRegression}. The 4 is the reciprocal of derivative of exp(x)/(1+exp(x))
# at x=0. 
U=rlogis(length(id_list), location = 0, scale = 1)
#U = o %*% o %*% rnorm(length(id_list)) #preference shock
covar=all_info2[which(all_info2$id %in% id_list),c('id','male','age','agpop','ricearea_2010','literacy','risk_averse','disaster_prob')]
covar=as.matrix(covar)
beta=c(0.0216,0.0044,-0.0089,0.00436,0.0891,0.1091,0.0018) * 4 #from Table 2 column1 in the paper
Utility_index= covar[,2:8] %*% beta


#choose the intercept such that the base line take-up rate is 35 percent on average
intercept_seq=seq(0,5,by=0.02)
base_line= c()
for (intercept in intercept_seq){
  base_line=c(base_line,mean(exp(Utility_index-intercept)/(1+exp(Utility_index-intercept))))
}
#the intercept equals to 2.26
final=intercept_seq[which.min(abs(base_line-0.35))]
Utility_index = Utility_index - final
Utility_index= cbind(covar[,1],Utility_index)


#1st Exposure Mapping: First Round, Simple Info (FR-S)
Y_sim2[,1]= Utility_index[,2] > U
mean(Y_sim2[,1]) 

#2nd Exposure Mapping: First Round, Intensive Info (FR-I)
Y_sim2[,2]= (Utility_index[,2]+0.141/(0.35*0.65)) > U #0.141 from Table 2 Column 1 
mean(Y_sim2[,2]) 

#3rd Exposure Mapping: Second Round Simple, No friend in the First Round
Y_sim2[,3]= Utility_index[,2] > U
mean(Y_sim2[,3])

#4th Exposure Mapping: Second Round Simple, has a friend in FR-S but not FR-I
Y_sim2[,4]=(Utility_index[,2]/(0.35*0.65)) > U #The paper claims no effects
mean(Y_sim2[,4])

#5th Exposure Maping: Second Round Simple, has only one friend in FR-I
Y_sim2[,5]=(Utility_index[,2]+ 0.0970/(0.35*0.65)) > U #Table 2 column 5
mean(Y_sim2[,5])

#6th Exposure Maping: Second Round Simple, has two friends in FR-I
Y_sim2[,6]=(Utility_index[,2]+ 0.177/(0.35*0.65)) > U #from Table 2 Column 5 
mean(Y_sim2[,6])

#7th Exposure Maping: Second Round Simple, has more than two friends in FR-I
Y_sim2[,7]=(Utility_index[,2]+ 0.177/(0.35*0.65)) > U #from Table 2 Column 5 
mean(Y_sim2[,7])

#8th Exposure Mapping: Second Round Intensive, No friend in the First Round
Y_sim2[,8]=  (Utility_index[,2]+0.141/(0.35*0.65)) > U #0.141 from Table 2 Column 1 
mean(Y_sim2[,8])

#9th Exposure Mapping: Second Round Intensive, has a friend in FR-S but not FR-I
Y_sim2[,9]=(Utility_index[,2]+ (0.141+0.001)/(0.35*0.65)) > U #0.141 from Table 2 Column 1 
mean(Y_sim2[,9])

#10th Exposure Mapping: Second Round Intensive, has only one friend in FR-I
Y_sim2[,10]=(Utility_index[,2]+ 0.177/(0.35*0.65)) > U #0.0970 1 from Table 2 Column 3 
mean(Y_sim2[,10])

#11th Exposure Maping: Second Round Intensive, has two friends in FR-I
Y_sim2[,11]=(Utility_index[,2]+ 0.177/(0.35*0.65)) > U #0.177  from Table 2 Column 3 
mean(Y_sim2[,11])

#12th Exposure Maping: Second Round Intensive, has more than two friends in FR-I
Y_sim2[,12]=(Utility_index[,2]+ 0.177/(0.35*0.65)) > U #0.177  from Table 2 Column 3 
mean(Y_sim2[,12])

Y_sim2=cbind(covar[,c('id')],Y_sim2)

data_save_path='/Users/haoge/Desktop'
data_save_path=paste0(data_save_path,'/Y_sim2.csv')
write.csv(Y_sim2,file=data_save_path)

data_save_path='/Users/haoge/Desktop'
data_save_path=paste0(data_save_path,'/Utilitiy_index.csv')
write.csv(Utility_index,data_save_path)
#write.csv(o,file='/Users/haoge/Downloads/sociomatrix.csv')
