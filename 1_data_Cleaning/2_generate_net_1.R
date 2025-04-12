rm(list=ls())
library(statnet)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

################################################
####Set working directory########################
################################################
setwd("/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData")

all_people_info_path = '/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData/all_info2.csv'
network_info_path = '/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData/net_complete.RData'

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
net %v% "village" = villages #add village information to the network object

#natural village
id_list=net %v% "vertex.names"
nat_villages= all_info[match(id_list,all_info$id),"address"] 
sum(is.na(nat_villages)) 
net %v% "nat_village" = nat_villages

#47 villages
length(unique(net%v%"village"))
#number of households
network.size(net)


#max in and out degree degree
max(degree(net,cmode="outdegree"))
max(degree(net,cmode="indegree"))
max(degree(net,cmode="freeman"))

#inspect graph
components(net,'weak')

#number of households
network.size(net)

#max degree
max(degree(net,gmode="graph"))

##############################################################################
########################component statistics##################################
##############################################################################
#plotting a component
membership=component.dist(net,connected='weak')$membership
membership_info= cbind(1:length(membership),membership)#info for graphing
net %v% "membership" = membership

compare=cbind(get.vertex.attribute(net,"membership"),get.vertex.attribute(net,"nat_village"))
colnames(compare)=c('membership','nat_village')


#######################Code randomization group##########################
#natural village
nat_village=unique(all_info$address)
all_info['group']=9999
all_info['ricearea_2010_per_capita'] = all_info['ricearea_2010'] / all_info['agpop']
group_index=0
for (i in 1:length(nat_village)){
  print(i)
  nat_village_temp=nat_village[i]
  all_info_temp=all_info[which(all_info$address==nat_village_temp),]
  
  print(paste0('Village ',nat_village_temp,' has ',nrow(all_info_temp), ' households.'))
  

  cutoff_agpop_temp=  quantile(all_info_temp$agpop,0.5,na.rm = TRUE)
  cutoff_ricearea_temp = quantile(all_info_temp$ricearea_2010_per_capita,0.5,na.rm = TRUE)
  
  #extract people from four quantile groups
  index1_temp = all_info_temp[ ((all_info_temp$agpop < cutoff_agpop_temp)) &  ((all_info_temp$ricearea_2010_per_capita < cutoff_ricearea_temp)),'id'] #upper left
  index2_temp = all_info_temp[ ((all_info_temp$agpop >= cutoff_agpop_temp)) &  ((all_info_temp$ricearea_2010_per_capita < cutoff_ricearea_temp)),'id'] # lower left
  index3_temp = all_info_temp[ ((all_info_temp$agpop < cutoff_agpop_temp)) &  ((all_info_temp$ricearea_2010_per_capita >= cutoff_ricearea_temp)),'id']  # upper right
  index4_temp = all_info_temp[ ((all_info_temp$agpop >= cutoff_agpop_temp)) &  ((all_info_temp$ricearea_2010_per_capita >= cutoff_ricearea_temp)),'id']  # lower right
  

  
  all_info[all_info$id %in% index1_temp,'group']=group_index + 1
  all_info[all_info$id %in% index2_temp,'group']=group_index + 2
  all_info[all_info$id %in% index3_temp,'group']=group_index + 3
  all_info[all_info$id %in% index4_temp,'group']=group_index + 4
  
  group_index=group_index+4
  
  
}
#############################################################################

compare=cbind(get.vertex.attribute(net,"membership"),get.vertex.attribute(net,"vertex.names")) 
colnames(compare)=c('membership','id') #households belongs to which component
all_info2=merge(all_info,as.data.frame(compare),by='id')

##############################################################################################################################
#######Some groups are too small (less than 4 people) we merge some network components and small groups ######################
##############################################################################################################################
#merge some graph components
all_info2[which(all_info2$membership==2),'membership']=1
all_info2[which(all_info2$membership==3),'membership']=1
all_info2[which(all_info2$membership==4),'membership']=1

all_info2[which(all_info2$membership==9),'membership']=6
all_info2[which(all_info2$membership==6),'membership']=12

all_info2[which(all_info2$membership==14),'membership']=15
all_info2[which(all_info2$membership==15),'membership']=13

all_info2[which(all_info2$membership==17),'membership']=10

all_info2[which(all_info2$membership==22),'membership']=24
all_info2[which(all_info2$membership==28),'membership']=26
all_info2[which(all_info2$membership==29),'membership']=26
all_info2[which(all_info2$membership==32),'membership']=31

all_info2[which(all_info2$membership==33),'membership']=35
all_info2[which(all_info2$membership==34),'membership']=35

all_info2[which(all_info2$membership==39),'membership']=41

membership_vec=sort(unique(all_info2$membership))


##for each component we combind small groups (with size<4) to a larger group.
for (i in membership_vec)  {
  print(i)
  all_info2_temp=all_info2[which(all_info2$membership==i),]
  
  #all groups and their size
  dist_temp=table(all_info2_temp[,'group'])
  
  #group of type 1 (low low)
  dist_temp1=dist_temp[as.numeric(names(dist_temp)) %%4==1]
  dist_temp1=dist_temp1[order(dist_temp1)]
  
  if (!all(dist_temp1>=4)){
  index1A=min((1:length(dist_temp1))[cumsum(dist_temp1)>=4]) #the first entry where sum is larger than 4
  index1B=max((1:length(dist_temp1))[dist_temp1<4]) #groups who have sizes smaller than 4 
  list1=names(dist_temp1)[1:max(index1A,index1B)]
  all_info2[which(all_info2$group %in% list1),'group']=min(as.numeric(list1)) #coarsen to a larger group
  
  }
  
  dist_temp2=dist_temp[as.numeric(names(dist_temp)) %%4==2]
  dist_temp2=dist_temp2[order(dist_temp2)]
  
  if (!all(dist_temp2>=4)){
  index2A=min((1:length(dist_temp2))[cumsum(dist_temp2)>=4]) #the first entry where sum is larger than 4
  index2B=max((1:length(dist_temp2))[dist_temp2<4]) #groups who have sizes smaller than 4 
  list2=names(dist_temp2)[1:max(index2A,index2B)]
  all_info2[which(all_info2$group %in% list2),'group']=min(as.numeric(list2)) #coarsen to a larger group
  

  }
  
  dist_temp3=dist_temp[as.numeric(names(dist_temp)) %%4==3]
  dist_temp3=dist_temp3[order(dist_temp3)]
  
  if (!all(dist_temp3>=4)){
  index3A=min((1:length(dist_temp3))[cumsum(dist_temp3)>=4]) #the first entry where sum is larger than 4
  index3B=max((1:length(dist_temp3))[dist_temp3<4]) #groups who have sizes smaller than 4 
  list3=names(dist_temp3)[1:max(index3A,index3B)]
  all_info2[which(all_info2$group %in% list3),'group']=min(as.numeric(list3)) #coarsen to a larger group
  
  }
  
  dist_temp4=dist_temp[as.numeric(names(dist_temp)) %%4==0]
  dist_temp4=dist_temp4[order(dist_temp4)]
  if (!all(dist_temp4>=4)){
  index4A=min((1:length(dist_temp4))[cumsum(dist_temp4)>=4]) #the first entry where sum is larger than 4
  index4B=max((1:length(dist_temp4))[dist_temp4<4]) #groups who have sizes smaller than 4 
  list4=names(dist_temp4)[1:max(index4A,index4B)]
  all_info2[which(all_info2$group %in% list4),'group']=min(as.numeric(list4)) #coarsen to a larger group
  

}
  
}

#We merge the following group so that the treatment probability for some people in the group is not too small
#Note this is an issue of homophily. Some people nominates all people in the same strata. 

# #Merge group 125 and 127
# all_info2[all_info2['group']==125,]
# all_info2[all_info2['group']==127,]
# 
# all_info2[all_info2['group']==125,'group']=127
# 
# #Merge group 442 and 443
# all_info2[all_info2['group']==442,]
# all_info2[all_info2['group']==443,]
# all_info2[all_info2['group']==442,'group']=443
# 
# #Merge group 648 and 406
# all_info2[all_info2['group']==648,]
# all_info2[all_info2['group']==406,]
# all_info2[all_info2['group']==648,'group']=406
# 
# #Merge group 574 and 575
# all_info2[all_info2['group']==575,]
# all_info2[all_info2['group']==574,]
# all_info2[all_info2['group']==575,'group']=574
# 
# #Merge group 492, 493 and 494
# all_info2[all_info2['group']==492,]
# all_info2[all_info2['group']==493,]
# all_info2[all_info2['group']==494,]
# all_info2[all_info2['group']==492,'group']=494
# all_info2[all_info2['group']==493,'group']=494
# 
# #Merge group 239, 240 and 241
# all_info2[all_info2['group']==239,]
# all_info2[all_info2['group']==240,]
# all_info2[all_info2['group']==241,]
# all_info2[all_info2['group']==239,'group']=241
# all_info2[all_info2['group']==240,'group']=241
# 
# #Merge group 170,171
# all_info2[all_info2['group']==170,]
# all_info2[all_info2['group']==171,]
# all_info2[all_info2['group']==170,'group']=171
# 
# #Merge group 490,491,492
# all_info2[all_info2['group']==489,]
# all_info2[all_info2['group']==490,]
# all_info2[all_info2['group']==492,]
# all_info2[all_info2['group']==489,'group']=492
# all_info2[all_info2['group']==489,'group']=490
# 
# #Merge group 660 and 661
# all_info2[all_info2['group']==660,]
# all_info2[all_info2['group']==661,]
# all_info2[all_info2['group']==660,'group']=661
# 
# #Merge group 134 and 135
# all_info2[all_info2['group']==135,]
# all_info2[all_info2['group']==134,]
# all_info2[all_info2['group']==135,'group']=134
# 
# #Merge group 13 and 14
# all_info2[all_info2['group']==13,]
# all_info2[all_info2['group']==14,]
# all_info2[all_info2['group']==13,'group']=14
# 
# 
# #Merge group 88 86 87
# all_info2[all_info2['group']==88,]
# all_info2[all_info2['group']==86,]
# all_info2[all_info2['group']==87,]
# all_info2[all_info2['group']==88,'group']=86
# all_info2[all_info2['group']==87,'group']=86
# 
# #Merge group 667 and 666
# all_info2[all_info2['group']==667,]
# all_info2[all_info2['group']==666,]
# all_info2[all_info2['group']==667,'group']=666
# 
# #Merge group
# all_info2[all_info2['group']==236,]
# all_info2[all_info2['group']==234,]
# all_info2[all_info2['group']==236,'group']=234

#code assignment patterns 
all_info2['assignment']=99

group_index=unique(all_info2$group)
for (i in group_index){
  
  ind=all_info2$group==i
  nobs=sum(ind)
  
  integer_part= nobs %/% 4

  temp=rep(c(4,3,2,1),integer_part+1)
  temp=temp[1:nobs]
  
  all_info2[ind,'assignment']=temp #order does not matter
  
}

##############################################################################################################################
##################Embed Treatment Information to the Network Object###########################################################
##############################################################################################################################
id_list=net %v% "vertex.names"
group= all_info2[match(id_list,all_info2$id),"group"] 
assignment=all_info2[match(id_list,all_info2$id),"assignment"] 
membership=all_info2[match(id_list,all_info2$id),"membership"] 
sum(is.na(villages)) 
net %v% "group" = group
net %v% "assignment"=assignment
net %v% "membership"= membership
id_list=unique(network_info$id) #people who have network information

##############################################################################################################################
##################Output CSV files############################################################################################
##############################################################################################################################
write.csv(all_info2,all_people_info_path)
save(net,all_info2,id_list,file=network_info_path)
