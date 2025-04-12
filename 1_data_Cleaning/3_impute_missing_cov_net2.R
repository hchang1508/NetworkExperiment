rm(list=ls())

#####################################################
########Set Directory and Load Data##################
#####################################################
#load files
setwd("/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData")
network_info_path = '/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData/net_complete_natvillage.RData'
load(network_info_path)
target_people=read.csv('target_info.csv')

#####################################################
########Check Missing Patterns#######################
#####################################################
#covar data, dealing with missing values
covars = c('id','male','age', 'agpop', 'ricearea_2010','risk_averse', 'disaster_prob', 'literacy', 'takeup_survey')
#check missing patterns
percent=rep(NA,length(covars))
for (i in 1:length(covars)){
  
  temp=covars[i] 

  percent[i]=sum(is.na(all_info2[which(all_info2$id %in% target_people$id),temp]))/nrow(all_info2)
}
names(percent)=covars

percent 

############################################################################
########Impute Missing Variables with nonmissing data#######################
############################################################################
for (i in 1:length(covars)){
  
  #ith covariates
  temp=covars[i]
  print(temp)
  #rows with missing values
  missing = is.na(all_info2[,temp])
  
  #averages
  avg = mean(all_info2[,temp],na.rm = TRUE) 
  
  if (percent[i]<=0.1){
    #if missing pattern is not severe, replace missing values with averages
    all_info2[missing ,temp] = avg
  }else{
    
    #if missing pattern is severe, missing indicator method
    cov_name = temp
    print(cov_name)
    new_col = paste0(cov_name,'_missing')
    all_info2[,new_col] = missing * avg
    all_info2[missing,temp] = 0
    
  }
  
}

save(net,all_info2,id_list,file=network_info_path)
