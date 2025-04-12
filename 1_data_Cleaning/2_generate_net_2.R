rm(list=ls())
library(statnet)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

################################################
####Set working directory########################
################################################
setwd("/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData")

all_people_info_path = '/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData/all_info2_nat.csv'
network_info_path = '/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData/net_complete_natvillage.RData'

################################################################
####Load Preprocessed Dat from script Data_Cleaning.R ##########
################################################################
network_info=read.csv("network_info.csv") #network information
all_info=read.csv("all_info.csv") #all information (including covariates)

#create a net object
relation_data <- network_info[,c('id','network_id')]
net <- network(relation_data,matrix.type="edgelist")

#administrative villages
id_list=net %v% "vertex.names"
villages= all_info[match(id_list,all_info$id),"village"] 
sum(is.na(villages)) #check missing values
net %v% "village" = villages #
#natural village
id_list=net %v% "vertex.names"
nat_villages= all_info[match(id_list,all_info$id),"address"] 
sum(is.na(nat_villages)) 
net %v% "nat_village" = nat_villages

#membership
membership=component.dist(net,connected='weak')$membership
membership_info= cbind(1:length(membership),membership)#info for graphing
net %v% "membership" = membership

#######################Code randomization group##########################
#natural village
nat_village=unique(all_info$address)
all_info['group']=9999
group_index=0
for (i in 1:length(nat_village)){
  print(i)
  
  #nat_village_level_randomization
  nat_village_temp=nat_village[i]
  all_info[which(all_info$address==nat_village_temp),'group']=group_index

  group_index=group_index+1
}
#############################################################################
compare=cbind(get.vertex.attribute(net,"membership"),get.vertex.attribute(net,"vertex.names"))
colnames(compare)=c('membership','id')
all_info2=merge(all_info,compare,by='id')

#merge some clusters
all_info2[which(all_info2$membership==2),'membership']=1
all_info2[which(all_info2$membership==3),'membership']=1
all_info2[which(all_info2$membership==4),'membership']=1

all_info2[which(all_info2$membership==9),'membership']=6
all_info2[which(all_info2$membership==6),'membership']=12

all_info2[which(all_info2$membership==14),'membership']=15
all_info2[which(all_info2$membership==15),'membership']=13

all_info2[which(all_info2$membership==22),'membership']=24
all_info2[which(all_info2$membership==28),'membership']=26
all_info2[which(all_info2$membership==29),'membership']=26
all_info2[which(all_info2$membership==32),'membership']=31

all_info2[which(all_info2$membership==33),'membership']=35
all_info2[which(all_info2$membership==34),'membership']=35

all_info2[which(all_info2$membership==39),'membership']=41

#check membership and group
compare=cbind(all_info2[,'membership'],all_info2[,'group'])
nat_villages= unique(compare[,2])
result=matrix(0,length(nat_villages),3)
result[,1]=nat_villages
for (i in 1:length(nat_villages)){
  temp_village=nat_villages[i]
  membership_temp=compare[which(compare[,2]==temp_village),1]
    
  if (unique(membership_temp)>1){
      print(i)
  }
  result[i,2]=length(membership_temp)
  result[i,3]=unique(membership_temp)
}
colnames(result) = c('nat_village','size','membership')
head(result[order(result[,2],decreasing=FALSE),])

######################################################################
###########For groups smaller than 10, merge with a larger group######
######################################################################

#merge 79 80
all_info2[all_info2$group==80,'group']=79

#merge 47 and 48
all_info2[all_info2$group==47,'group']=48

group_index=unique(all_info2$group)
for (i in group_index){
  
  ind=all_info2$group==i
  nobs=sum(ind)

  if (nobs<10){
    print(i)
    print(nobs)
    
  }
  integer_part= nobs %/% 10
  
  pattern=c(4,3,2,1,4,3,4,3,4,3)
  temp=rep(pattern,integer_part+1)
  temp=temp[1:nobs]
  
  all_info2[ind,'assignment']=temp #order does not matter
  
}

id_list=net %v% "vertex.names"
group= all_info2[match(id_list,all_info2$id),"group"] 
assignment=all_info2[match(id_list,all_info2$id),"assignment"] 
membership=all_info2[match(id_list,all_info2$id),"membership"] 
sum(is.na(villages)) 
net %v% "group" = group
net %v% "assignment"=assignment
net %v% "membership"= membership
id_list=unique(network_info$id)

####save#####
write.csv(all_info2,all_people_info_path)
save(net,all_info2,id_list,file=network_info_path)


#########Old Code########
#number of people in each villages
num_people=table(all_info2[,'address'])
num_people=cbind(as.numeric(rownames(num_people)),as.numeric(num_people))
colnames(num_people)=c('nat_village','num')
###Check two hop villages
nat_village_connect=merge(relation_data,all_info2[,c('id','address')],by=c('id'))
sum(is.na(nat_village_connect$address))
colnames(nat_village_connect)=c('id','network_id','ego_village')
nat_village_connect=merge(nat_village_connect,all_info2[,c('id','address')],by.x=c('network_id'),by.y=c('id'))
colnames(nat_village_connect)=c('id','network_id','ego_village','friend_village')

one_hop = distinct(nat_village_connect,nat_village_connect$ego_village,nat_village_connect$friend_village)
colnames(one_hop)=c('point_from','point_to')
two_hop = merge(one_hop,one_hop, by.x=c('point_to'),by.y=c('point_from'))

two_hop=two_hop[,c(1,3)]
colnames(two_hop)=c('point_from','point_to_two_hop')
two_hop_number_people=merge(two_hop,num_people,by.x=c('point_to_two_hop'), by.y=c('nat_village'))

two_hop_sums=as.data.frame(
two_hop_number_people %>% group_by(point_from) %>% summarise(sum(num)))

