rm(list=ls())
library(dplyr)
library(statnet)
setwd('/Users/haoge/Library/CloudStorage/GoogleDrive-hgchang1508@gmail.com/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/Ranalysis/')


#import data, start with 23243 relations
data=read.csv('network.csv')
covar_raw = read.csv('covar.csv')


#covar data, dealing with missing values
covars = c('id','male','age', 'agpop',
           'educ', 'ricearea_2010',
           'rice_inc', 'disaster_yes', 'disaster_loss' ,
           'risk_averse', 'disaster_prob', 'literacy', 'understanding','takeup_survey')
covar=covar_raw[,covars]

#only 2425 complete cases
sum(complete.cases(covar))
percent=rep(NA,ncol(covar))
for (i in 1:ncol(covar)){
  percent[i]=sum(is.na(covar[,i]))/nrow(covar)
  
}
percent #disater_loss missing around 50%
names(percent)=covars

#if a column has data missing less than 10 percent, impute with averages
#if a column has data missing more than 10 percent, use missing indicator method
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
covar=cbind(covar,covar_raw$village)
colnames(covar)[16]='village'

#number of unique households: 4984
#number of households' friends: 7248
# Some vertices have info missing, especially network informations
length(unique(data$id))
length(unique(data$network_id)) 
id_list=unique(data$id)
#4902 households have survey informations
length(unique(covar$id))


# remove isolated nodes: for those exposure mappings are not identified
# 4661 households remain
# 22812 rows (edges)
data=data[!is.na(data$network_id),]
length(unique(data$id))
total_rows=nrow(data)

#keep households with survey info available; 4587 households remain
data=data[which(data$id %in% covar$id),]
length(unique(c(data$id)))

#99 is the missing value, 4586 household remain
data=data[which(data$network_id!=99),]
length(unique(c(data$id)))
nrow(data)

#remove self-loop: 4586 households and 17015 rows remain 
data=data[which(data$id != data$network_id ),]
length(unique(data$id))
nrow(data)

#remove repeated friend nomination: 4586 households and 16997 rows remain 
data=distinct(data,data$id,data$network_id,.keep_all=TRUE)
length(unique(data$id))
final_rows=nrow(data)

#How many edgeds dropped? Arond 609
total_rows-final_rows

#id_list is the population of interest
id_list=unique(data$id)

#compare statistics with Table1, Panel A in the paper:
apply(covar,2,mean,na.rm=TRUE)

#each household has on average 4.841474 friends
#histogram of network degress
hist(table(data$id))
mean(table(data$id))

#outcomes
outcome=c('id','takeup_survey')
Y=covar[,outcome] #actual outcome

covars = c('id','male','age', 'agpop',
           'educ', 'ricearea_2010',
           'rice_inc', 'disaster_yes', 'disaster_loss' ,'disaster_loss_missing' ,
           'risk_averse', 'disaster_prob', 'literacy', 'understanding','village')
X=covar[,covars]
assig =c('id','delay','intensive')
assignment = data[,assig]


#create network project and save 
relation_data <- data[,c('id','network_id')]
net <- network(relation_data,matrix.type="edgelist")
data_save_path='/Users/haoge/Library/CloudStorage/GoogleDrive-hgchang1508@gmail.com/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/Ranalysis/data/truncated_graph/'

data_save_path=paste0(data_save_path,'/net_complete.RData')
save(net,X,Y,id_list,relation_data,file=data_save_path)



