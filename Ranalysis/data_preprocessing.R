rm(list=ls())
library(dplyr)
library(statnet)
setwd('/Volumes/GoogleDrive/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/Ranalysis/')


#import data, startwing with 23243 relations
data=read.csv('network.csv')
covar = read.csv('covar.csv')


#covar data, dealing with missing values
covars = c('id','male','age', 'agpop',
           'educ', 'ricearea_2010',
           'rice_inc', 'disaster_yes', 'disaster_loss' ,
           'risk_averse', 'disaster_prob', 'literacy', 'understanding','takeup_survey')
covar=covar[,covars]
sum(complete.cases(covar))#only 2425 complete cases
percent=rep(NA,ncol(covar))
for (i in 1:ncol(covar)){
  percent[i]=sum(is.na(covar[,i]))/nrow(covar)

}
names(percent)=covars
percent #disater_loss missing around 50%

#if a column has less then 10 percent data missing, replace missing data with the mean;
#if more than 10 percent, use missing indicator method;
for (i in 1:ncol(covar)){

  missing = is.na(covar[,i])
  avg = mean(covar[,i],na.rm = TRUE)
  if (percent[i]<=0.1){
        covar[missing ,i] = avg
  }else{
    cov_name = covars[i]
    print(cov_name)
    new_col = paste0(cov_name,'_missing')
    covar[,new_col] = missing * avg
    covar[is.na(covar[,i]),i] = 0
  }

}

#number of unique households: 4984
#number of households' friends: 7248
length(unique(data$network_id)) 
length(unique(data$id))
id_list=unique(data$id)
#4902 households have survey informations
length(unique(covar$id))

#extract ego-centric networks for households whose 
#friendship information are observed; 
# also remove isolated nodes
# 4914 households remain
data=data[which(data$network_id %in% id_list ),]
length(unique(c(data$id)))

#keep households with both survey and network available; 4832 households remain
data=data[which(data$network_id %in% covar$id & data$id %in% covar$id),]
data=data[which(data$id %in% covar$id),]
length(unique(c(data$id,data$network_id)))

#How many edgeds dropped?

#99 is the missing value
data=data[which(data$network_id!=99),]

#compare statistics with Table1, Panel A in the paper:
apply(covar,2,mean,na.rm=TRUE)


#remove self-loop: 4529 households and 17015 rows remain 
data=data[which(data$id != data$network_id ),]
length(unique(data$id))
nrow(data)

#remove repeated friend nomination: 4529 households and 16997 rows remain 
data=distinct(data,data$id,data$network_id,.keep_all=TRUE)
length(unique(data$id))
nrow(data)

#end wih 16997 relations


id_list=unique(data$id)

#each household has on average 3.753 friends
#histogram of network degress
hist(table(data$id))
mean(table(data$id))

#outcomes
outcome=c('id','takeup_survey')
Y=covar[,outcome]

covars = c('id','male','age', 'agpop',
          'educ', 'ricearea_2010',
          'rice_inc', 'disaster_yes', 'disaster_loss' ,'disaster_loss_missing' ,
          'risk_averse', 'disaster_prob', 'literacy', 'understanding')
X=covar[,covars]
assig =c('id','delay','intensive')
assignment = data[,assig]

#create network project and save 
relation_data <- data[,c('id','network_id')]
net <- network(relation_data,matrix.type="edgelist")
data_save_path='/Volumes/GoogleDrive/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/Ranalysis/data/truncated_graph'

data_save_path=paste0(data_save_path,'/net.RData')
save(net,X,Y,id_list,file=data_save_path)

#######################Impute Simulation###################################




